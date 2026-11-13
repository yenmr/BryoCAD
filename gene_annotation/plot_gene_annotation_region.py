#!/usr/bin/env python3
"""Fig. 3b-style pre-CRF probabilities and GFF3 comparison. See --help.

Requires numpy, xarray, zarr, matplotlib (use the GeneCAD environment).
GFF3 inputs must use the same assembly as the prediction stores.
"""
import argparse
import gzip
from collections import defaultdict
from pathlib import Path
import sys
from urllib.parse import unquote

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

FEATURES = ['intergenic', 'five_prime_utr', 'cds', 'intron', 'three_prime_utr']
NAMES = ['Intergenic', "5′ UTR", 'CDS', 'Intron', "3′ UTR"]
COLORS = ['#ffffff', '#a8cf75', '#78b6ce', '#b29bc8', '#63bdb7']


def load_probabilities(path, chromosome, strand, start, stop):
    """Read only ROI logits; scan coordinates in bounded blocks, including unsorted RC."""
    root = Path(path)
    if (root / '.zgroup').exists() or (root / 'zarr.json').exists():
        stores = [root]
    else:
        stores = sorted(root.glob('predictions.*.zarr')) + sorted(root.glob('segment.*.zarr'))
    if not stores:
        raise ValueError(f'No prediction Zarr stores found in {root}')
    result = np.full((stop - start, 5), np.nan, dtype=np.float32)
    counts = np.zeros(stop - start, dtype=np.int32)
    for store in stores:
        try:
            ds = xr.open_zarr(store, group=strand, consolidated=True)
        except (KeyError, FileNotFoundError):
            continue
        with ds:
            stored_chrom = ds.attrs.get('chromosome_id')
            if stored_chrom and str(stored_chrom) != chromosome:
                raise ValueError(f'{store}: chromosome_id={stored_chrom!r}, expected {chromosome!r}')
            if 'feature_logits' not in ds:
                raise ValueError(f'{store}: missing feature_logits')
            for left in range(0, ds.sizes['sequence'], 1_000_000):
                coords = ds.sequence.isel(sequence=slice(left, left + 1_000_000)).values
                hits = np.flatnonzero((coords >= start) & (coords < stop))
                if not len(hits):
                    continue
                positions = coords[hits].astype(np.int64) - start
                logits = ds.feature_logits.isel(sequence=left + hits).sel(
                    feature=FEATURES).transpose('sequence', 'feature').values
                if np.any(counts[positions]) or len(np.unique(positions)) != len(positions):
                    raise ValueError('Duplicate ROI coordinates across/within stores; choose the correct stores.')
                result[positions] = logits
                counts[positions] += 1
    if np.any(counts != 1):
        missing = np.flatnonzero(counts == 0)
        raise ValueError(f'{len(missing)} ROI bases missing; first missing 0-based coordinate: {start + missing[0]}')
    if not np.isfinite(result).all():
        raise ValueError('Non-finite feature logits')
    result -= result.max(axis=1, keepdims=True)
    np.exp(result, out=result)
    result /= result.sum(axis=1, keepdims=True)
    return result


def attributes(text):
    return {k: v for item in text.split(';') if '=' in item
            for k, v in [item.split('=', 1)]}


def gff_labels(path, chromosome, strand, start, stop, requested_ids):
    """Canonical coding transcript per gene; conflicting overlaps become unknown (-1).

    Infer introns between exons and UTRs from exonic sequence flanking the CDS.
    Noncoding-only models are excluded because the five-state scheme has no ncRNA state.
    """
    transcripts = {}
    children = defaultdict(list)
    aliases = {'5utr': 'five_prime_utr', '3utr': 'three_prime_utr',
               'five_prime_utr': 'five_prime_utr', 'three_prime_utr': 'three_prime_utr'}
    # Detect gzip by its signature, including BGZF and files without a .gz suffix.
    with open(path, 'rb') as probe:
        compressed = probe.read(2) == b'\x1f\x8b'
    opener = gzip.open if compressed else open
    with opener(path, 'rt', encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('##FASTA'):
                break
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'{path}: expected 9-column GFF3')
            chrom, _, kind, lo, hi, _, direction, _, raw = fields
            if chrom != chromosome or direction != strand:
                continue
            attrs = attributes(raw)
            kind = kind.lower()
            lo, hi = int(lo) - 1, int(hi)
            if kind in ('mrna', 'transcript'):
                tid = unquote(attrs.get('ID', ''))
                if tid:
                    transcripts[tid] = unquote(attrs.get('Parent', tid))
            elif kind in ('exon', 'cds', 'intron', 'utr', 'five_prime_utr', 'three_prime_utr', '5utr', '3utr'):
                # Split before URL decoding, since escaped commas may belong to IDs.
                for parent in attrs.get('Parent', '').split(','):
                    if parent:
                        children[unquote(parent)].append((aliases.get(kind, kind), lo, hi))
    coding = {tid: records for tid, records in children.items()
              if any(k == 'cds' for k, _, _ in records)}
    if requested_ids:
        missing = set(requested_ids) - set(coding)
        if missing:
            raise ValueError(f'{path}: coding transcript IDs not found on selected chromosome/strand: {sorted(missing)}')
        chosen = list(dict.fromkeys(requested_ids))
    else:
        by_gene = defaultdict(list)
        for tid in coding:
            by_gene[transcripts.get(tid, tid)].append(tid)
        def rank(tid):
            records = coding[tid]
            return (-sum(b-a for k,a,b in records if k == 'cds'),
                    -sum(b-a for k,a,b in records if k == 'exon'), tid)
        chosen = [sorted(ids, key=rank)[0] for ids in by_gene.values()]
    labels = np.zeros(stop - start, dtype=np.int8)
    used = []
    for tid in chosen:
        records = coding[tid]
        if max(b for _,_,b in records) <= start or min(a for _,a,_ in records) >= stop:
            continue
        used.append(tid)
        local = np.zeros(stop - start, dtype=np.int8)
        cds = [(a,b) for k,a,b in records if k == 'cds']
        cstart, cstop = min(a for a,b in cds), max(b for a,b in cds)
        exons = sorted((a,b) for k,a,b in records if k == 'exon')
        if not exons:
            exons = sorted((a,b) for k,a,b in records if k in ('cds','utr','five_prime_utr','three_prime_utr'))
        merged = []
        for a,b in exons:
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(b, merged[-1][1]))
            else:
                merged.append((a,b))
        def paint(a,b,value):
            a,b = max(a,start), min(b,stop)
            if a < b:
                local[a-start:b-start] = value
        for (_, b), (a, _) in zip(merged, merged[1:]):
            paint(b,a,3)
        for a,b in merged:
            paint(a,min(b,cstart),1 if strand == '+' else 4)
            paint(max(a,cstop),b,4 if strand == '+' else 1)
        for kind in ('intron','five_prime_utr','three_prime_utr','cds'):
            for k,a,b in records:
                if k == kind:
                    paint(a,b,FEATURES.index(kind))
        conflict = (labels != 0) & (local != 0) & (labels != local)
        insert = (labels == 0) & (local != 0)
        labels[insert] = local[insert]
        labels[conflict] = -1
    if not used:
        print(f'WARNING: {path}: no selected coding transcript overlaps ROI; track is intergenic.', file=sys.stderr)
    return labels, used


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--predictions', required=True, help='One prediction Zarr store, or directory containing prediction/segment stores')
    parser.add_argument('--prediction-gff', required=True, help='Predicted GFF3 (e.g. final GeneCAD GFF)')
    parser.add_argument('--annotation-gff', required=True, help='Reference GFF3')
    parser.add_argument('--chromosome', required=True)
    parser.add_argument('--strand', choices=['+','-','positive','negative'], default='+')
    parser.add_argument('--start', type=int, required=True, help='1-based inclusive genomic start')
    parser.add_argument('--end', type=int, required=True, help='1-based inclusive genomic end')
    parser.add_argument('--prediction-transcript', action='append', help='Use specific predicted transcript ID; repeat for multiple models')
    parser.add_argument('--annotation-transcript', action='append', help='Use specific reference transcript ID; repeat for multiple models')
    parser.add_argument('--orientation', choices=['genomic','transcription'], default='transcription', help='Transcription orientation reverses all tracks on the negative strand')
    parser.add_argument('--output', required=True, help='Output figure path (.png, .pdf, .svg)')
    parser.add_argument('--save-npz', help='Optional output of full-resolution probabilities, labels and error mask')
    parser.add_argument('--dpi', type=int, default=250)
    args = parser.parse_args()
    if args.start < 1 or args.end < args.start:
        parser.error('Require 1 <= start <= end')
    start, stop = args.start - 1, args.end
    strand = 'positive' if args.strand in ('+','positive') else 'negative'
    sign = '+' if strand == 'positive' else '-'
    prob = load_probabilities(args.predictions, args.chromosome, strand, start, stop)
    pred, pred_ids = gff_labels(args.prediction_gff, args.chromosome, sign, start, stop, args.prediction_transcript)
    ref, ref_ids = gff_labels(args.annotation_gff, args.chromosome, sign, start, stop, args.annotation_transcript)
    valid = (pred >= 0) & (ref >= 0)
    errors = np.where(valid, (pred != ref).astype(np.int8), -1)
    nerrors, nvalid = int((errors == 1).sum()), int(valid.sum())
    print(f'Predicted transcripts: {pred_ids}\nReference transcripts: {ref_ids}')
    print(f'Mismatched bases: {nerrors}/{nvalid}; ambiguous excluded: {len(valid)-nvalid}')
    if args.save_npz:
        Path(args.save_npz).parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(args.save_npz, positions_0based=np.arange(start,stop),
                            features=np.array(FEATURES), probabilities=prob,
                            prediction_labels=pred, annotation_labels=ref, errors=errors,
                            chromosome=args.chromosome, strand=sign,
                            prediction_transcripts=np.array(pred_ids, dtype=str),
                            annotation_transcripts=np.array(ref_ids, dtype=str))
    reverse = sign == '-' and args.orientation == 'transcription'
    if reverse:
        prob, pred, ref, errors = prob[::-1], pred[::-1], ref[::-1], errors[::-1]
    fig, axes = plt.subplots(4,1,figsize=(14,4.7),sharex=True,
                             gridspec_kw={'height_ratios':[5,0.65,0.65,0.65], 'hspace':0.13})
    im = axes[0].imshow(prob.T, aspect='auto', cmap='Greys', vmin=0,vmax=1,
                        interpolation='nearest', extent=[0,len(prob),4.5,-0.5])
    axes[0].set_yticks(range(5), NAMES)
    axes[0].set_title(f'{args.chromosome}:{args.start:,}–{args.end:,} ({sign} strand; {len(prob):,} bp)')
    cmap = ListedColormap(COLORS)
    cmap.set_bad('#b0b0b0')
    norm = BoundaryNorm(np.arange(-0.5,5.5), 5)
    for ax, labels, title in zip(axes[1:3], [pred,ref], ['Predictions','Annotations']):
        ax.imshow(np.ma.masked_less(labels[None,:],0), aspect='auto', cmap=cmap,
                  norm=norm, interpolation='nearest', extent=[0,len(prob),0,1])
        ax.set_yticks([])
        ax.set_ylabel(title, rotation=0, ha='right', va='center')
    ecmap = ListedColormap(['white','#333333'])
    ecmap.set_bad('#b0b0b0')
    axes[3].imshow(np.ma.masked_less(errors[None,:],0), aspect='auto', cmap=ecmap,
                    vmin=0,vmax=1, interpolation='nearest', extent=[0,len(prob),0,1])
    axes[3].set_yticks([])
    axes[3].set_ylabel('Errors',rotation=0,ha='right',va='center')
    ticks = np.linspace(0,len(prob),6)
    axes[3].set_xticks(ticks, [f'{x:,.0f}' for x in ticks])
    axes[3].set_xlabel('Position in window (bp; 5′ → 3′)' if args.orientation == 'transcription' else 'Position in window (bp; increasing genomic coordinate)')
    fig.subplots_adjust(left=0.12,right=0.83,bottom=0.16,top=0.88)
    cax = fig.add_axes([0.845,0.45,0.012,0.4])
    fig.colorbar(im,cax=cax,label='Probability (pre-CRF)')
    handles = [Patch(facecolor=c,edgecolor='#999999',label=n) for c,n in zip(COLORS,NAMES)]
    handles.append(Patch(facecolor='#b0b0b0',label='Ambiguous / excluded'))
    axes[1].legend(handles=handles,loc='upper left',bbox_to_anchor=(1.01,1),fontsize=8,frameon=False)
    Path(args.output).parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,dpi=args.dpi,bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {args.output}')


if __name__ == '__main__':
    main()
