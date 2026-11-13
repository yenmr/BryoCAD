#!/usr/bin/env python3
"""
Annotation location evaluator. No gene filtering, motifs or BUSCO; FASTA not required.

python evaluate_compatible.py --ref reference.gff3.gz --pred predictions.gff \
    --structure-policy utr --metrics-tsv metrics.tsv --audit-json audit.json

Fixed geometry options when explicit exon+UTR both exist: utr, exon, intersection.
Apply the same option to both inputs. Preserve all genes, transcripts and CDS.
Missing reference information yields NA comparison metrics; prediction counts
and raw feature statistics remain available. Source annotations are never edited.
"""

import argparse
import gzip
import os
import sys
import json
from collections import Counter
from urllib.parse import unquote
from collections import defaultdict


# ─────────────────────────────────────────────────────────────────────────────
# Shared utilities
# ─────────────────────────────────────────────────────────────────────────────


def parse_attributes(attr_str):
    """Parse GFF3 attribute string into a key-value dict."""
    attrs = {}
    for field in attr_str.strip().split(";"):
        if "=" in field:
            k, _, v = field.strip().partition("=")
            attrs[k.strip()] = v.strip()
    return attrs


def parse_gff(path, ftypes=None):
    """Normalize GFF3 feature names and preserve attributes; do not invent exons."""
    aliases = {'gene': 'gene', 'mrna': 'mRNA', 'transcript': 'mRNA',
               'cds': 'CDS', 'exon': 'exon', 'utr': 'UTR',
               'five_prime_utr': 'five_prime_UTR', 'three_prime_utr': 'three_prime_UTR'}
    with open(path, 'rb') as probe:
        compressed = probe.read(2) == b'\x1f\x8b'
    opener = gzip.open if compressed else open
    records = []
    with opener(path, 'rt', encoding='utf-8') as fh:
        for number, line in enumerate(fh, 1):
            if line.startswith('##FASTA'):
                break
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9:
                raise ValueError(f'{path}:{number}: expected 9-column GFF3 (GTF is not supported)')
            ft = aliases.get(parts[2].lower())
            if ft is None:
                continue
            lo, hi = int(parts[3]), int(parts[4])
            if lo < 1 or hi < lo:
                raise ValueError(f'{path}:{number}: invalid interval')
            attrs = {k: unquote(v) for k,v in parse_attributes(parts[8]).items()}
            parents_raw = parse_attributes(parts[8]).get('Parent', '')
            parents = [unquote(x) for x in parents_raw.split(',') if x]
            if ft not in ('gene', 'mRNA') and not parents:
                raise ValueError(f'{path}:{number}: {ft} lacks Parent; cannot attach feature safely')
            records.append(dict(chrom=parts[0], ftype=ft, start=lo, end=hi,
                                strand=parts[6], id=attrs.get('ID',''),
                                parent=parents[0] if len(parents)==1 else '',
                                parents=parents, attrs=attrs, phase=parts[7]))
    if not records:
        raise ValueError(f'{path}: no supported GFF3 features')
    return records


def annotation_models(records, allow_cds_only_exons=False):
    """Order-independent hierarchy: explicit gene Parent, geneID, then audited placeholder."""
    tx = {}
    gene_records = {}
    for r in records:
        if r['ftype'] == 'gene':
            if not r['id'] or r['id'] in gene_records:
                raise ValueError('Missing/duplicate gene ID')
            gene_records[r['id']] = r
        elif r['ftype'] == 'mRNA':
            tid = r['id']
            if not tid or tid in tx:
                raise ValueError(f'Missing/duplicate transcript ID: {tid!r}')
            if r['strand'] not in ('+','-'):
                raise ValueError(f'{tid}: unknown strand; cannot score strand-aware structures')
            parents = r['parents']
            if len(parents) > 1:
                raise ValueError(f'{tid}: multiple gene parents are unsupported')
            geneid = r['attrs'].get('geneID') or r['attrs'].get('gene_id')
            if parents and geneid and parents[0] != geneid:
                raise ValueError(f'{tid}: conflicting Parent and geneID')
            gid = parents[0] if parents else geneid
            source = 'Parent' if parents else ('geneID' if geneid else 'unknown')
            if not gid:
                gid = '__unresolved_gene__:' + tid
            tx[tid] = dict(gene=gid, gene_source=source, chrom=r['chrom'], strand=r['strand'],
                           start=r['start'], end=r['end'], attrs=r['attrs'], features=defaultdict(list), warnings=[])
    for r in records:
        if r['ftype'] in ('gene','mRNA'):
            continue
        for tid in r['parents']:
            if tid not in tx:
                raise ValueError(f"{r['ftype']} Parent={tid!r} has no mRNA/transcript record")
            model = tx[tid]
            if (r['chrom'],r['strand']) != (model['chrom'],model['strand']):
                raise ValueError(f'{tid}: child chromosome/strand differs from transcript')
            if r['start'] < model['start'] or r['end'] > model['end']:
                raise ValueError(f'{tid}: child outside transcript bounds')
            model['features'][r['ftype']].append((r['start'],r['end']))
    for tid,t in tx.items():
        g = gene_records.get(t['gene'])
        if g and ((g['chrom'],g['strand']) != (t['chrom'],t['strand']) or
                  t['start'] < g['start'] or t['end'] > g['end']):
            raise ValueError(f'{tid}: gene bounds/chromosome/strand mismatch')
        cds = sorted(set(t['features'].get('CDS',[])))
        if any(a[1] >= b[0] for a,b in zip(cds,cds[1:])):
            raise ValueError(f'{tid}: overlapping CDS segments within one transcript')
        explicit = t['features'].get('exon',[])
        utr = sum((t['features'].get(k,[]) for k in ('five_prime_UTR','three_prime_UTR','UTR')),[])
        if explicit:
            t['exons'] = tuple(merge_intervals(explicit))
            t['exon_source'] = 'explicit_exon'
            if any(not any(a<=c and d<=b for a,b in t['exons']) for c,d in cds+utr):
                t['warnings'].append('CDS/UTR not contained in explicit exons; explicit feature locations preserved independently')
        elif cds and utr:
            t['exons'] = tuple(merge_intervals(cds+utr))
            t['exon_source'] = 'CDS_plus_UTR_reconstruction'
        elif cds and allow_cds_only_exons:
            t['exons'] = tuple(merge_intervals(cds))
            t['exon_source'] = 'CDS_only_as_annotated'
        else:
            t['exons'] = None
            t['exon_source'] = 'unavailable'
    if not tx:
        raise ValueError('No mRNA/transcript records parsed; unsupported hierarchy')
    return tx


def prepare_models(records, prediction=False, policy="utr"):
    """Apply a fixed geometry policy without gene selection, for all organisms.

    A reference with exon/UTR records supports exonic reconstruction; remaining
    CDS-only transcripts are scored as annotated without UTR. A reference with
    only CDS has no evidence for full-exon/UTR evaluation. Predictions are always
    scored as emitted, so omitted predicted UTRs can incur FN.
    """
    has_structure = any(r['ftype'] in ('exon','five_prime_UTR','three_prime_UTR','UTR') for r in records)
    models=annotation_models(records, allow_cds_only_exons=prediction or has_structure)
    apply_structure_policy(models,policy)
    return models


def subtract_intervals(spans, other):
    result=[]
    other=merge_intervals(other)
    for a,b in merge_intervals(spans):
        for c,d in other:
            if d<a:
                continue
            if c>b:
                break
            if a<c:
                result.append((a,c-1))
            a=max(a,d+1)
            if a>b:
                break
        if a<=b:
            result.append((a,b))
    return result


def apply_structure_policy(models,policy):
    """Resolve geometry consistently without removing genes or transcripts."""
    if policy not in ('utr','exon','intersection'):
        raise ValueError(f'Unknown structure policy: {policy}')
    for tid,t in models.items():
        cds=merge_intervals(t['features'].get('CDS',[]))
        exons=merge_intervals(t['features'].get('exon',[]))
        explicit={'five_prime_utr':list(t['features'].get('five_prime_UTR',[])),
                  'three_prime_utr':list(t['features'].get('three_prime_UTR',[]))}
        has_utr=any(t['features'].get(k) for k in ('five_prime_UTR','three_prime_UTR','UTR'))
        if not cds:
            t['resolved_utr']=explicit
            t['resolution_source']='no_CDS_policy_not_applicable'
            continue
        lo,hi=cds[0][0],cds[-1][1]
        left_type,right_type=('five_prime_utr','three_prime_utr') if t['strand']=='+' else ('three_prime_utr','five_prime_utr')
        for a,b in t['features'].get('UTR',[]):
            if b<lo: explicit[left_type].append((a,b))
            elif a>hi: explicit[right_type].append((a,b))
            else: t['warnings'].append('Generic UTR overlaps coding envelope; cannot assign UTR type')
        inferred={'five_prime_utr':[],'three_prime_utr':[]}
        for a,b in exons:
            if a<lo: inferred[left_type].append((a,min(b,lo-1)))
            if b>hi: inferred[right_type].append((max(a,hi+1),b))
        # CDS takes precedence over UTR at overlapping bases in every policy.
        explicit={k:subtract_intervals(v,cds) for k,v in explicit.items()}
        if exons and has_utr:
            if policy=='utr':
                chosen=explicit
                resolved=merge_intervals(cds+sum(chosen.values(),[]))
            elif policy=='exon':
                chosen=inferred
                resolved=merge_intervals(exons+cds)
            else:
                chosen={k:subtract_intervals(explicit[k],subtract_intervals(explicit[k],inferred[k])) for k in explicit}
                resolved=merge_intervals(cds+sum(chosen.values(),[]))
            if any(merge_intervals(explicit[k])!=merge_intervals(inferred[k]) for k in explicit):
                t['warnings'].append('Explicit UTR differs from exon+CDS inference; resolved using '+policy)
            t['exons']=tuple(resolved)
            t['resolution_source']='exon_and_UTR:'+policy
        elif exons:
            chosen=inferred
            t['exons']=tuple(merge_intervals(exons+cds))
            t['resolution_source']='exon_plus_CDS_only'
        elif has_utr:
            chosen=explicit
            t['exons']=tuple(merge_intervals(cds+sum(chosen.values(),[])))
            t['resolution_source']='CDS_plus_UTR_only'
        else:
            chosen={'five_prime_utr':[],'three_prime_utr':[]}
            t['resolution_source']='CDS_only'
        t['resolved_utr']={k:merge_intervals(v) for k,v in chosen.items()}
    return models


def region_intervals(models):
    """Strand-aware interval unions after fixed-policy coordinate resolution."""
    if any('resolved_utr' not in t for t in models.values()):
        apply_structure_policy(models,'utr')
    regions = {name: defaultdict(list) for name in ('exonic','cds','five_prime_utr','three_prime_utr')}
    for tid,t in models.items():
        key = t['chrom'],t['strand']
        cds = t['features'].get('CDS',[])
        regions['cds'][key].extend(cds)
        if t['exons'] is not None:
            regions['exonic'][key].extend(t['exons'])
        for name,spans in t.get('resolved_utr',{}).items():
            regions[name][key].extend(spans)
    return regions


def interval_base_score(reference, prediction):
    """Union intervals per (sequence ID, strand); overlapping isoforms count once."""
    keys = set(reference)|set(prediction)
    tp = sum(overlap_bases(reference.get(k,[]),prediction.get(k,[])) for k in keys)
    nr = sum(bases_in_merged(v) for v in reference.values())
    np = sum(bases_in_merged(v) for v in prediction.values())
    return prf(tp,np-tp,nr-tp)


def models_to_genes(models, kind):
    genes = {}
    for tid,t in models.items():
        if kind == 'cds':
            intervals = t['features'].get('CDS',[])
            value = frozenset((t['chrom'],a,b,t['strand']) for a,b in intervals)
        elif kind == 'coding_structure':
            intervals = t['features'].get('CDS',[])
            value = tuple(merge_intervals(intervals))
        else:
            value = t['exons']
        if not value:
            continue
        gid = t['gene']
        g = genes.setdefault(gid, dict(chrom=t['chrom'],strand=t['strand'],transcripts={}))
        if (g['chrom'],g['strand']) != (t['chrom'],t['strand']):
            raise ValueError(f'Gene {gid}: inconsistent chromosome/strand')
        g['transcripts'][tid] = value
    return genes


def _build_gene_tx(records, cds_ftype='CDS'):
    return models_to_genes(annotation_models(records),'cds')


def build_ref_cds(records):
    genes = _build_gene_tx(records)
    by_chrom = defaultdict(set)
    for g in genes.values():
        for fs in g['transcripts'].values():
            by_chrom[g['chrom']].add(fs)
    return genes, by_chrom


def build_ref_exon(records):
    return models_to_genes(annotation_models(records),'exon')


def build_pred_exon(records):
    return build_ref_exon(records)



def merge_intervals(ivs):
    """Merge a list of (start, end) tuples."""
    if not ivs:
        return []
    s = sorted(ivs)
    merged = [list(s[0])]
    for a, b in s[1:]:
        if a <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return [tuple(x) for x in merged]


def bases_in_merged(ivs):
    """Return the total number of bases covered by a list of (start, end) intervals."""
    return sum(e - s + 1 for s, e in merge_intervals(ivs))


def overlap_bases(ivs_a, ivs_b):
    """Return the number of bases shared between two lists of (start, end) intervals."""
    ma, mb = merge_intervals(ivs_a), merge_intervals(ivs_b)
    tp = i = j = 0
    while i < len(ma) and j < len(mb):
        os, oe = max(ma[i][0], mb[j][0]), min(ma[i][1], mb[j][1])
        if os <= oe:
            tp += oe - os + 1
        if ma[i][1] < mb[j][1]:
            i += 1
        else:
            j += 1
    return tp


def introns_from_exons(sorted_exons, chrom, strand):
    """Derive intron coordinates from a sorted list of (start, end) exon tuples."""
    return [
        (chrom, sorted_exons[i][1] + 1, sorted_exons[i + 1][0] - 1, strand)
        for i in range(len(sorted_exons) - 1)
        if sorted_exons[i][1] + 1 <= sorted_exons[i + 1][0] - 1
    ]


def prf(tp, fp, fn):
    if min(tp,fp,fn) < 0:
        raise ValueError('Negative confusion count')
    return dict(TP=tp, FP=fp, FN=fn,
                precision=tp/(tp+fp) if tp+fp else None,
                recall=tp/(tp+fn) if tp+fn else None,
                f1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None)


def _fmt(value):
    return 'N/A' if value is None else f'{value:.4f}'




# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1 – CDS-Based Evaluation
# ─────────────────────────────────────────────────────────────────────────────


def _set_score(ref, pred):
    tp = len(ref & pred)
    return prf(tp,len(pred)-tp,len(ref)-tp)


def _entity_score(ref, pred):
    """Maximum one-to-one matching of entities sharing >=1 exact structure."""
    index = defaultdict(set)
    for rid, structures in ref.items():
        for structure in structures:
            index[structure].add(rid)
    edges = {pid: sorted(set().union(*(index[x] for x in structures)))
             for pid,structures in pred.items()}
    # Iterative augmenting paths avoid recursion limits in large duplicate clusters.
    owner, matched_pred = {}, {}
    for root in sorted(pred):
        queue, seen_pred, predecessor = [root], {root}, {}
        free = None
        for pid in queue:
            for rid in edges[pid]:
                if rid in predecessor:
                    continue
                predecessor[rid] = pid
                if rid not in owner:
                    free = rid
                    break
                other = owner[rid]
                if other not in seen_pred:
                    seen_pred.add(other)
                    queue.append(other)
            if free is not None:
                break
        while free is not None:
            pid = predecessor[free]
            previous = matched_pred.get(pid)
            owner[free], matched_pred[pid] = pid, free
            free = previous
    tp = len(owner)
    return prf(tp,len(pred)-tp,len(ref)-tp)


def cds_entities(genes):
    txs, loci, segments = {},{},set()
    for gid,g in genes.items():
        structures = set()
        for tid,fs in g['transcripts'].items():
            if fs:
                signature = tuple(sorted(fs))
                txs[tid] = {signature}
                structures.add(signature)
                segments.update(fs)
        if structures:
            loci[gid] = structures
    return txs,loci,segments


def eval_cds(ref_genes, ref_by_chrom, pred_genes):
    rt,rg,rs = cds_entities(ref_genes)
    pt,pg,ps = cds_entities(pred_genes)
    return dict(locus=_entity_score(rg,pg),tx=_entity_score(rt,pt),cds_ex=_set_score(rs,ps),
                n_ref_loci=len(rg),n_pred_loci=len(pg),n_ref_tx=len(rt),n_pred_tx=len(pt),
                n_ref_cds=len(rs),n_pred_cds=len(ps))


def structural_collections(genes):
    bases = defaultdict(list)
    exons,introns,chains = set(),set(),set()
    loci,txs,multi_loci = {},{},{}
    for gid,g in genes.items():
        ch,st = g['chrom'],g['strand']
        sigs,ic = set(),set()
        for tid,spans in g['transcripts'].items():
            if not spans:
                continue
            signature = (ch,st,tuple(spans))
            sigs.add(signature)
            txs[tid] = {signature}
            for a,b in spans:
                exons.add((ch,a,b,st))
                bases[ch,st].append((a,b))
            gaps = tuple(introns_from_exons(spans,ch,st))
            introns.update(gaps)
            if gaps:
                chains.add(gaps)
                ic.add(gaps)
        if sigs:
            loci[gid] = sigs
        if ic:
            multi_loci[gid] = ic
    return bases,exons,introns,chains,loci,txs,multi_loci


def eval_gffcompare(ref_genes, pred_genes):
    rb,re,ri,rc,rl,rt,rm = structural_collections(ref_genes)
    pb,pe,pi,pc,pl,pt,pm = structural_collections(pred_genes)
    keys = set(rb)|set(pb)
    tp = sum(overlap_bases(rb[k],pb[k]) for k in keys)
    rn = sum(bases_in_merged(x) for x in rb.values())
    pn = sum(bases_in_merged(x) for x in pb.values())
    return dict(base=prf(tp,pn-tp,rn-tp),exon=_set_score(re,pe),intron=_set_score(ri,pi),
                ichain=_set_score(rc,pc),locus_ic=_entity_score(rm,pm),
                full_tx=_entity_score(rt,pt),n_ref_loci=len(rl),n_ref_loci_multi=len(rm),
                n_ref_tx=len(rt),ref_bases=rn,n_pred_loci=len(pl),n_pred_loci_multi=len(pm),
                n_pred_tx=len(pt),pred_bases=pn,n_ref_exons=len(re),n_pred_exons=len(pe),
                n_ref_introns=len(ri),n_pred_introns=len(pi),n_ref_ichains=len(rc),n_pred_ichains=len(pc))




# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3 – Splice Site Analysis
# ─────────────────────────────────────────────────────────────────────────────


def load_fasta(fasta_path):
    """Load a FASTA file into a dict mapping sequence name to sequence string."""
    print(f"[INFO] Loading genome FASTA: {fasta_path}", file=sys.stderr)
    genome = {}
    chrom, seq = None, []
    opener = gzip.open if fasta_path.endswith(".gz") else open
    with opener(fasta_path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if chrom:
                    genome[chrom] = "".join(seq)
                chrom = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if chrom:
        genome[chrom] = "".join(seq)
    print(f"[INFO] Loaded {len(genome)} sequences.", file=sys.stderr)
    return genome


def _collect_tis_tts(cds_genes):
    """
    Extract TIS and TTS from CDS-based genes (output of _build_gene_tx).

    Sites are (chrom, pos, strand) tuples (1-based GFF coordinates).

    TIS – Translation Initiation Site: first coding base on the coding strand.
          For '+': leftmost  start of all CDS intervals in the transcript.
          For '-': rightmost end   of all CDS intervals in the transcript.
    TTS – Translation Termination Site: last coding base on the coding strand.
          For '+': rightmost end   of all CDS intervals.
          For '-': leftmost  start of all CDS intervals.

    Because CDS-based gene dicts store transcripts as frozensets of
    (chrom, start, end, strand) tuples, we handle them directly here.
    """
    tis_set, tts_set = set(), set()
    for gdata in cds_genes.values():
        ch = gdata["chrom"]
        st = gdata["strand"]
        for tx_data in gdata["transcripts"].values():
            if not tx_data:
                continue
            if isinstance(tx_data, frozenset):
                intervals = sorted((s, e) for (_, s, e, _) in tx_data)
            else:
                intervals = sorted(tx_data)
            if not intervals:
                continue
            first_s = intervals[0][0]
            last_e = intervals[-1][1]
            if st == "+":
                tis_set.add((ch, first_s, st))
                tts_set.add((ch, last_e, st))
            else:
                tis_set.add((ch, last_e, st))
                tts_set.add((ch, first_s, st))
    return tis_set, tts_set


def _collect_junctions(exon_genes):
    """
    Extract splice-junction donor and acceptor sites from exon-based genes
    (output of build_ref_exon / build_pred_exon).

    Donor    – last  exon base immediately before each intron.
    Acceptor – first exon base immediately after  each intron.
    """
    donor_set, acceptor_set = set(), set()
    for gdata in exon_genes.values():
        ch = gdata["chrom"]
        st = gdata["strand"]
        for exons in gdata["transcripts"].values():
            if not exons or len(exons) < 2:
                continue
            for i in range(len(exons) - 1):
                exon_end = exons[i][1]
                exon_start = exons[i + 1][0]
                if st == "+":
                    donor_set.add((ch, exon_end, st))
                    acceptor_set.add((ch, exon_start, st))
                else:
                    donor_set.add((ch, exon_start, st))
                    acceptor_set.add((ch, exon_end, st))
    return donor_set, acceptor_set


def eval_sites(ref_cds_genes, pred_cds_genes, ref_exon_genes, pred_exon_genes):
    """
    Compare predicted vs reference genomic sites and return per-type
    precision / recall / F1.

    TIS/TTS  are derived from CDS-based genes (exact translation boundaries).
    Donor/Acceptor are derived from exon-based genes (splice junctions).
    """
    ref_tis, ref_tts = _collect_tis_tts(ref_cds_genes)
    pred_tis, pred_tts = _collect_tis_tts(pred_cds_genes)
    ref_don, ref_acc = _collect_junctions(ref_exon_genes)
    pred_don, pred_acc = _collect_junctions(pred_exon_genes)

    def _prf_sets(pred_s, ref_s):
        tp = len(pred_s & ref_s)
        fp = len(pred_s) - tp
        fn = len(ref_s) - tp
        return prf(tp, fp, fn)

    return {
        "TIS": _prf_sets(pred_tis, ref_tis),
        "TTS": _prf_sets(pred_tts, ref_tts),
        "donor": _prf_sets(pred_don, ref_don),
        "acceptor": _prf_sets(pred_acc, ref_acc),
        "n_ref_tis": len(ref_tis),
        "n_ref_tts": len(ref_tts),
        "n_ref_donor": len(ref_don),
        "n_ref_acceptor": len(ref_acc),
        "n_pred_tis": len(pred_tis),
        "n_pred_tts": len(pred_tts),
        "n_pred_donor": len(pred_don),
        "n_pred_acceptor": len(pred_acc),
    }


def format_site_report(r):
    """Format the site-level error breakdown (Section 5) as a string."""
    lines = [
        "=" * 62,
        "SECTION 5 – Site-Level Error Breakdown",
        "  Each site is a strand-aware (chrom, pos, strand) tuple.",
        "  TIS  = Translation Initiation Site (CDS start on coding strand)",
        "  TTS  = Translation Termination Site (CDS end on coding strand)",
        "  Donor    = last exon base before each intron",
        "  Acceptor = first exon base after each intron",
        "=" * 62,
        "",
        f"  {'Site':<12}  {'Ref':>7}  {'Pred':>7}  {'TP':>7}  {'FP':>7}  {'FN':>7}  "
        f"{'Prec':>7}  {'Rec':>7}  {'F1':>7}",
        "  " + "-" * 76,
    ]
    for site_key, label in [
        ("TIS", "TIS"),
        ("TTS", "TTS"),
        ("donor", "Junc.Donor"),
        ("acceptor", "Junc.Acc."),
    ]:
        m = r[site_key]
        ref = r[f"n_ref_{site_key.lower()}"]
        pred = r[f"n_pred_{site_key.lower()}"]
        lines.append(
            f"  {label:<12}  {ref:>7}  {pred:>7}  {m['TP']:>7}  {m['FP']:>7}  {m['FN']:>7}  "
            f"{_fmt(m['precision']):>7}  {_fmt(m['recall']):>7}  {_fmt(m['f1']):>7}"
        )
    lines.append("")
    return "\n".join(lines)


def _row(title, m):
    """Format a single metric block (TP/FP/FN + P/R/F1) as a list of lines."""
    return [
        f"--- {title} ---",
        f"  TP: {m['TP']}   FP: {m['FP']}   FN: {m['FN']}",
        f"  Precision : {_fmt(m['precision'])}",
        f"  Recall    : {_fmt(m['recall'])}",
        f"  F1        : {_fmt(m['f1'])}",
        "",
    ]


def format_cds_report(r):
    """Format the CDS-based evaluation results (Section 1) as a string."""
    lines = [
        "=" * 62,
        "SECTION 1 – CDS-Based Evaluation (UTRs ignored)",
        "  Gene is correct if its CDS chain matches ANY ref isoform.",
        "=" * 62,
        "",
        f"  Reference : {r['n_ref_loci']} loci | {r['n_ref_tx']} transcripts | {r['n_ref_cds']} unique CDS exons",
        f"  Predicted : {r['n_pred_loci']} loci | {r['n_pred_tx']} transcripts | {r['n_pred_cds']} unique CDS exons",
        "",
    ]
    for title, key in [
        ("CDS-exon-level", "cds_ex"),
        ("Locus-level", "locus"),
        ("Transcript-level", "tx"),
    ]:
        lines += _row(title, r[key])
    return "\n".join(lines)


def format_gffcompare_report(r):
    """Format the gffcompare-style evaluation results (Section 2) as a string."""
    lines = [
        "=" * 62,
        "SECTION 2 – gffcompare-Style Evaluation (full exon incl. UTRs)",
        "  Counts are unique intervals per chromosome.",
        "  Intron chain / Locus[IC]: multi-exon transcripts only.",
        "=" * 62,
        "",
        f"  Reference : {r['n_ref_loci']} loci ({r['n_ref_loci_multi']} multi-exon) | "
        f"{r['n_ref_tx']} transcripts | {r['ref_bases']:,} bases",
        f"              {r['n_ref_exons']} exons | {r['n_ref_introns']} introns | "
        f"{r['n_ref_ichains']} intron chains",
        f"  Predicted : {r['n_pred_loci']} loci ({r['n_pred_loci_multi']} multi-exon) | "
        f"{r['n_pred_tx']} transcripts | {r['pred_bases']:,} bases",
        f"              {r['n_pred_exons']} exons | {r['n_pred_introns']} introns | "
        f"{r['n_pred_ichains']} intron chains",
        "",
    ]
    for title, key in [
        ("Exon level (start/end exact)", "exon"),
        ("Intron level (splice junctions)", "intron"),
        ("Intron chain level (= transcript interior match)", "ichain"),
        ("Locus level [intron-chain based, multi-exon only]", "locus_ic"),
    ]:
        lines += _row(title, r[key])
    return "\n".join(lines)


def main():
    parser=argparse.ArgumentParser(description='Compare annotation locations only: CDS, exons, UTR bases and splice-site coordinates. No FASTA, motif or BUSCO analysis is required.')
    parser.add_argument('--ref',required=True,help='Reference GFF3 or GFF3.gz')
    parser.add_argument('--pred',help='Prediction GFF3 (optional with --audit-only)')
    parser.add_argument('--output',help='Human-readable report; default stdout')
    parser.add_argument('--metrics-tsv',help='Machine-readable location metrics')
    parser.add_argument('--audit-json',help='Feature counts, input warnings and coordinate policy')
    parser.add_argument('--structure-policy',choices=['utr','exon','intersection'],default='utr',help='Fixed geometry rule when both explicit exons and UTRs exist; CDS always retained')
    parser.add_argument('--audit-only',action='store_true')
    parser.add_argument('--seqid-policy',choices=['strict','union'],default='strict',help='Union retains all single-sided sequence features as FP/FN; assumes same assembly')
    parser.add_argument('--seqid',action='append',help='Explicitly restrict both inputs to these seqids; repeat as needed')
    parser.add_argument('--fasta',help='Optional bounds check only; never used for motifs or BUSCO')
    parser.add_argument('--skip-busco',action='store_true',help='Accepted for old commands; BUSCO is always disabled')
    parser.add_argument('--cpu',type=int,default=1,help='Accepted for old commands; location evaluation runs in one process')
    parser.add_argument('--reference-structure',choices=['auto','complete','cds-only'],help=argparse.SUPPRESS)
    args=parser.parse_args()
    if args.reference_structure is not None:
        print('[WARNING] --reference-structure is deprecated and ignored. Use --structure-policy for the fixed coordinate rule.',file=sys.stderr)
    run_compatible_evaluation(args)



def raw_feature_base_counts(records):
    spans=defaultdict(lambda:defaultdict(list))
    for r in records:
        spans[r['ftype']][r['chrom'],r['strand']].append((r['start'],r['end']))
    return {ft:sum(bases_in_merged(v) for v in spans[ft].values()) for ft in ('CDS','exon','five_prime_UTR','three_prime_UTR')}


def input_audit(records, models):
    return dict(feature_counts=dict(Counter(r['ftype'] for r in records)),
                warning_transcripts=sum(bool(t['warnings']) for t in models.values()),
                warning_examples=[dict(transcript=tid,messages=t['warnings']) for tid,t in models.items() if t['warnings']][:30],
                sequence_ids=sorted({r['chrom'] for r in records}),
                transcript_count=len(models),
                coding_transcripts=sum(bool(t['features'].get('CDS')) for t in models.values()),
                gene_link_sources=dict(Counter(t['gene_source'] for t in models.values())),
                exon_sources=dict(Counter(t['exon_source'] for t in models.values())),
                resolution_sources=dict(Counter(t.get('resolution_source','unresolved') for t in models.values())),
                raw_feature_bp=raw_feature_base_counts(records),
                partial_attributes=sum(any(k.lower() in ('partial','start_range','end_range','start_not_found','end_not_found') for k in t['attrs']) for t in models.values()))


def run_compatible_evaluation(args):
    ref_records = parse_gff(args.ref)
    if not args.pred:
        if not args.audit_only:
            raise ValueError('--pred is required unless --audit-only is set')
        if args.seqid:
            requested = set(args.seqid)
            missing = requested - {r['chrom'] for r in ref_records}
            if missing:
                raise ValueError(f'Reference seqid absent: {sorted(missing)}')
            ref_records = [r for r in ref_records if r['chrom'] in requested]
        models = prepare_models(ref_records,policy=args.structure_policy)
        genome = load_fasta(args.fasta) if args.fasta else None
        if genome is not None:
            for r in ref_records:
                if r['chrom'] not in genome or r['end'] > len(genome[r['chrom']]):
                    raise ValueError(f"Reference {r['chrom']}:{r['start']}-{r['end']} outside FASTA")
        region_intervals(models)
        audit = dict(reference=input_audit(ref_records,models), prediction=None,
                     gene_selection='none; all genes and predictions retained',structure_policy=args.structure_policy,
                     fasta_bounds_checked=genome is not None,
                     reference_path=os.path.abspath(args.ref), coordinate_policy='Fixed '+args.structure_policy+' geometry; preserve CDS; do not filter genes')
        report = 'REFERENCE-ONLY AUDIT (no accuracy metrics computed)\n'+json.dumps(audit,indent=2)+'\n'
        if args.audit_json:
            with open(args.audit_json,'w') as f: json.dump(audit,f,indent=2)
        if args.output:
            with open(args.output,'w') as f: f.write(report)
        else:
            print(report)
        return
    pred_records = parse_gff(args.pred)
    requested = set(args.seqid or [])
    if requested:
        for label,records in [('reference',ref_records),('prediction',pred_records)]:
            missing = requested - {r['chrom'] for r in records}
            if missing:
                raise ValueError(f'{label}: requested seqid absent: {sorted(missing)}')
        ref_records = [r for r in ref_records if r['chrom'] in requested]
        pred_records = [r for r in pred_records if r['chrom'] in requested]
    ref_models = prepare_models(ref_records,policy=args.structure_policy)
    # Evaluate all emitted predictions: an omitted predicted UTR is a potential FN,
    # not a reason to discard the prediction or mark the whole comparison NA.
    pred_models = prepare_models(pred_records, prediction=True,policy=args.structure_policy)
    ref_chroms = {r['chrom'] for r in ref_records}
    pred_chroms = {r['chrom'] for r in pred_records}
    # Do not automatically drop reference chromosomes, since that inflates recall.
    if ref_chroms != pred_chroms and args.seqid_policy == 'strict':
        raise ValueError(f'Sequence ID sets differ: reference-only={sorted(ref_chroms-pred_chroms)}, prediction-only={sorted(pred_chroms-ref_chroms)}. Use --seqid explicitly for a subset, and verify assembly names.')
    genome = load_fasta(args.fasta) if args.fasta else None
    if genome is not None:
        for label,records in [('reference',ref_records),('prediction',pred_records)]:
            for r in records:
                if r['chrom'] not in genome or r['end'] > len(genome[r['chrom']]):
                    raise ValueError(f"{label}: {r['chrom']}:{r['start']}-{r['end']} outside FASTA")
    ref_regions,pred_regions = region_intervals(ref_models),region_intervals(pred_models)
    audit = dict(reference=input_audit(ref_records,ref_models),
                 prediction=input_audit(pred_records,pred_models),
                 gene_selection='none; all genes and predictions retained',structure_policy=args.structure_policy,
                 scope=sorted(ref_chroms | pred_chroms), seqid_policy=args.seqid_policy, reference_only_seqids=sorted(ref_chroms-pred_chroms), prediction_only_seqids=sorted(pred_chroms-ref_chroms), fasta_bounds_checked=genome is not None,
                 reference_path=os.path.abspath(args.ref),prediction_path=os.path.abspath(args.pred),
                 coordinate_policy='Fixed '+args.structure_policy+' geometry; preserve CDS; do not filter genes')
    if args.audit_json:
        with open(args.audit_json,'w') as f:
            json.dump(audit,f,indent=2)
    parts = ['SECTION 0 – Input audit and metric availability',json.dumps(audit,indent=2),
             'Missing gene records may be reconstructed from explicit Parent/geneID links.\n'
             'CDS-only reference files do not establish full exon/UTR annotation.\n'
             'CDS interval endpoints are not validated biological start/stop sites; phase is not scored.\n'
             'Explicit exons / CDS+UTR reconstruction provide annotated structures, not proof of completeness.\n'
             'Predicted CDS-only models are scored as emitted: no predicted UTR.\n'
             'One fixed structure policy is applied to both reference and prediction.\n'
             'No genes/transcripts/predictions are excluded for annotation inconsistency.\n'
             'Without both exon and UTR evidence, use available geometry without conflict resolution.\n'
             'Accuracy is agreement with reference annotation, not independent experimental validation.']
    if args.audit_only:
        report='\n\n'.join(parts)+'\n'
        if args.output:
            with open(args.output,'w') as f: f.write(report)
        else: print(report)
        return
    ref_cds = models_to_genes(ref_models,'cds')
    pred_cds = models_to_genes(pred_models,'cds')
    ref_coding = models_to_genes(ref_models,'coding_structure')
    pred_coding = models_to_genes(pred_models,'coding_structure')
    ref_full = models_to_genes(ref_models,'exon')
    pred_full = models_to_genes(pred_models,'exon')
    rows=[]
    def unavailable(section,metric,reason,ref_count=None,pred_count=None):
        rows.append(dict(section=section,metric=metric,status='NA',reason=reason,ref_count=ref_count,pred_count=pred_count))
    def add(section,metric,result):
        rows.append(dict(section=section,metric=metric,status='available',reason='',ref_count=result['TP']+result['FN'],pred_count=result['TP']+result['FP'],**result))
    def known_genes(models,kind):
        return all(t['gene_source']!='unknown' for t in models.values()
                   if (t['features'].get('CDS') if kind=='cds' else t['exons']))
    if ref_cds:
        cds_res=eval_cds(ref_cds,None,pred_cds)
        locus_ok=known_genes(ref_models,'cds') and known_genes(pred_models,'cds')
        if not locus_ok:
            cds_res['locus'] = dict(TP='NA',FP='NA',FN='NA',precision=None,recall=None,f1=None)
        parts.append(format_cds_report(cds_res)+'\nCDS transcript/locus matches use maximum one-to-one matching; duplicate predictions cannot produce negative FN.')
        for key in ('cds_ex','tx','locus'):
            if key=='locus' and not locus_ok:
                unavailable('CDS',key,'Missing gene grouping: no Parent or geneID for at least one coding transcript')
            else: add('CDS',key,cds_res[key])
        coding_res=eval_gffcompare(ref_coding,pred_coding)
        block=['SECTION 1B – Coding-region-only structure (NOT full exon/transcript accuracy)']
        for key,label in [('intron','CDS-supported intron gaps'),('ichain','CDS-supported intron chains')]:
            block+=_row(label,coding_res[key])
            add('coding_structure',key,coding_res[key])
        parts.append('\n'.join(block))
    else:
        reason='Reference has no CDS; coding metrics cannot be evaluated'
        parts.append('SECTION 1 – CDS-Based Evaluation\n  N/A: '+reason)
        for key in ('cds_ex','tx','locus'): unavailable('CDS',key,reason)
        for key in ('intron','ichain'): unavailable('coding_structure',key,reason)
    # A single global metric over a mixed annotation set would silently penalize unavailable UTRs.
    missing_ref=sum(t['exons'] is None for t in ref_models.values())
    missing_pred=sum(t['exons'] is None for t in pred_models.values())
    full_ok=bool(ref_full) and missing_ref==0 and missing_pred==0 
    if full_ok:
        full_res=eval_gffcompare(ref_full,pred_full)
        locus_ok=known_genes(ref_models,'exon') and known_genes(pred_models,'exon')
        if not locus_ok:
            full_res['locus_ic']=dict(TP='NA',FP='NA',FN='NA',precision=None,recall=None,f1=None)
        parts.append(format_gffcompare_report(full_res)+'\n'+'\n'.join(_row('Transcript level (all exon boundaries exact)',full_res['full_tx'])))
        for key in ('exon','intron','ichain','locus_ic','full_tx'):
            if key=='locus_ic' and not locus_ok:
                unavailable('full_structure',key,'Missing gene grouping')
            else: add('full_structure',key,full_res[key])
    else:
        reason=f'Full exonic structure unavailable or disabled (observed-coordinate policy): reference={missing_ref}/{len(ref_models)}, prediction={missing_pred}/{len(pred_models)} transcripts. CDS-only reference files provide no full-exon/UTR annotation evidence.'
        parts.append('SECTION 2 – Full exon/transcript comparison\n  N/A: '+reason)
        predicted_summary=eval_gffcompare({},pred_full)
        parts.append(f"Predicted structures retained: {predicted_summary['n_pred_tx']} transcripts, {predicted_summary['n_pred_exons']} exons, {predicted_summary['n_pred_introns']} introns. Reference comparison is NA.")
        count_keys=dict(exon='n_pred_exons',intron='n_pred_introns',ichain='n_pred_ichains',locus_ic='n_pred_loci_multi',full_tx='n_pred_tx')
        for key in ('exon','intron','ichain','locus_ic','full_tx'):
            unavailable('full_structure',key,reason,pred_count=predicted_summary[count_keys[key]])
    block=['BASE-LEVEL REGION COMPARISONS',
           'Unit: bases. Union by (sequence ID, strand); isoform overlaps counted once.',
           'exonic and UTR bases use the selected structure policy; cds uses original CDS coordinates.',
           'Base accuracy is reported only in base_level: exonic, cds, five_prime_utr, three_prime_utr.']
    for metric in ('exonic','cds','five_prime_utr','three_prime_utr'):
        available = bool(ref_cds) if metric=='cds' else full_ok
        if available:
            result=interval_base_score(ref_regions[metric],pred_regions[metric])
            add('base_level',metric,result)
            block+=_row(metric,result)
        else:
            reason='Reference CDS unavailable' if metric=='cds' else 'Reference lacks full exonic/UTR evidence; CDS-only reference files cannot establish UTR locations'
            pred_count=sum(bases_in_merged(v) for v in pred_regions[metric].values())
            unavailable('base_level',metric,reason,pred_count=pred_count)
            block.append(f'{metric}: predicted bases={pred_count}; reference count/accuracy=N/A – {reason}')
    parts.append('\n'.join(block))
    # Coding endpoints are measurable, but a CDS record alone cannot establish complete TIS/TTS.
    block=['SECTION 5 – Site comparisons',
           'CDS.5prime/CDS.3prime = annotated coding endpoints; NOT validated TIS/TTS.']
    if ref_cds:
        site_res=eval_sites(ref_cds,pred_cds,ref_coding,pred_coding)
        for key,label in [('TIS','CDS.5prime'),('TTS','CDS.3prime'),('donor','Coding.Junc.Donor'),('acceptor','Coding.Junc.Acceptor')]:
            block+=_row(label,site_res[key]); add('sites',label,site_res[key])
    else:
        for label in ('CDS.5prime','CDS.3prime','Coding.Junc.Donor','Coding.Junc.Acceptor'):
            unavailable('sites',label,'Reference has no CDS')
        block.append('Coding endpoint/junction comparisons: N/A (reference has no CDS).')
    if full_ok:
        ds,ac=_collect_junctions(ref_full); pds,pac=_collect_junctions(pred_full)
        for label,rs,ps in [('Full.Junc.Donor',ds,pds),('Full.Junc.Acceptor',ac,pac)]:
            result=_set_score(rs,ps); block+=_row(label,result); add('sites',label,result)
    else:
        pd,pa=_collect_junctions(pred_full)
        for label,ps in [('Full.Junc.Donor',pd),('Full.Junc.Acceptor',pa)]:
            unavailable('sites',label,'Full reference/predicted exon structure unavailable',pred_count=len(ps))
    parts.append('\n'.join(block))
    parts.append('SEQUENCE ANALYSIS: disabled. Location metrics require no FASTA, motif matching or BUSCO.')
    if args.metrics_tsv:
        import csv
        fields=['section','metric','status','reason','TP','FP','FN','precision','recall','f1','ref_count','pred_count']
        with open(args.metrics_tsv,'w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); writer.writeheader()
            for row in rows:
                writer.writerow({k:('NA' if row.get(k) is None else row[k]) for k in fields})
    report='\n\n'.join(parts)+'\n'
    if args.output:
        with open(args.output,'w') as f: f.write(report)
        print(f'[INFO] Report written to {args.output}',file=sys.stderr)
    else: print(report)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as exc:
        print(f'[FATAL] {exc}',file=sys.stderr)
        sys.exit(2)
