#!/usr/bin/env python3
"""Select reproducible whole contigs and matching GFF3 from a genome manifest.

Standard-library only. Selection is annotation-independent: retain sequences at
least --min-contig-bp long, order by a deterministic seeded hash, then take a
prefix reaching --target-bp. The last contig is never cut. Larger targets are
nested prefixes with the same seed, genome ID and minimum size.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from contextlib import contextmanager
from urllib.parse import unquote

VERSION = '1'


def open_text(path):
    with open(path, 'rb') as f:
        compressed = f.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt', encoding='utf-8') if compressed else open(path, encoding='utf-8')


@contextmanager
def gzip_text(path):
    # Suppress timestamps and temporary filenames in gzip headers.
    with open(path, 'wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding='utf-8') as text:
                yield text


def fasta_inventory(path):
    lengths = {}
    current = None
    with open_text(path) as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].split()
                if not name:
                    raise ValueError('Empty FASTA identifier')
                current = name[0]
                if current in lengths:
                    raise ValueError(f'Duplicate FASTA identifier: {current}')
                lengths[current] = 0
            elif line.strip():
                if current is None:
                    raise ValueError('Sequence before first FASTA header')
                lengths[current] += len(''.join(line.split()))
    if not lengths or not sum(lengths.values()):
        raise ValueError('No nonempty FASTA sequences')
    return lengths


def select_contigs(lengths, genome_id, target_bp, min_bp, seed):
    eligible = [name for name, n in lengths.items() if n >= min_bp and n > 0]
    def rank(name):
        key = json.dumps([seed, genome_id, name], ensure_ascii=False, separators=(',', ':'))
        return hashlib.sha256(key.encode()).digest(), name
    eligible.sort(key=rank)
    chosen = []
    total = 0
    for name in eligible:
        chosen.append(name)
        total += lengths[name]
        if total >= target_bp:
            break
    return chosen, eligible


def extract_fasta(source, destination, selected):
    selected = set(selected)
    found = set()
    write = False
    total = 0
    with open_text(source) as src, gzip_text(destination) as dst:
        for line in src:
            if line.startswith('>'):
                name = line[1:].split()[0]
                write = name in selected
                if write:
                    found.add(name)
                    dst.write(line.rstrip('\r\n') + '\n')
            elif write and line.strip():
                seq = ''.join(line.split())
                dst.write(seq + '\n')
                total += len(seq)
    if found != selected:
        raise ValueError(f'Selected FASTA records missing: {sorted(selected-found)}')
    return total


def extract_gff(source, destination, selected, lengths):
    selected = set(selected)
    all_seqids = set()
    annotated = set()
    selected_counts = {}
    selected_ids = set()
    selected_parents = set()
    n_selected = n_total = 0
    with open_text(source) as src, gzip_text(destination) as dst:
        dst.write('##gff-version 3\n')
        for number, line in enumerate(src, 1):
            if line.startswith('##FASTA'):
                break  # Separate FASTA is exported; embedded FASTA is not copied.
            if line.startswith('##gff-version'):
                continue
            if line.startswith('##sequence-region'):
                fields = line.split()
                if len(fields) >= 2 and fields[1] in selected:
                    dst.write(line.rstrip('\r\n') + '\n')
                continue
            if line.startswith('#'):
                dst.write(line.rstrip('\r\n') + '\n')
                continue
            if not line.strip():
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'GFF line {number}: expected 9 columns')
            seqid = fields[0]
            all_seqids.add(seqid)
            n_total += 1
            if seqid not in selected:
                continue
            start, end = int(fields[3]), int(fields[4])
            if not (1 <= start <= end <= lengths[seqid]):
                raise ValueError(f'GFF line {number}: {seqid}:{start}-{end} outside FASTA bounds')
            dst.write(line.rstrip('\r\n') + '\n')
            annotated.add(seqid)
            n_selected += 1
            selected_counts[fields[2]] = selected_counts.get(fields[2], 0) + 1
            at = dict(x.split('=', 1) for x in fields[8].split(';') if '=' in x)
            if 'ID' in at:
                selected_ids.add(unquote(at['ID']))
            selected_parents.update(unquote(x) for x in at.get('Parent', '').split(',') if x)
    if not all_seqids:
        raise ValueError('GFF contains no feature records')
    if not (all_seqids & lengths.keys()):
        raise ValueError('FASTA/GFF share no sequence IDs; check assembly and chromosome naming')
    return dict(gff_total_records=n_total, gff_selected_records=n_selected,
                selected_feature_counts=selected_counts,
                selected_sequences_without_gff=sorted(selected-annotated),
                gff_sequences_absent_from_fasta=sorted(all_seqids-lengths.keys()),
                selected_unresolved_parent_ids=sorted(selected_parents-selected_ids))


def resolve_input(value, manifest_dir):
    path = Path(value).expanduser()
    return (path if path.is_absolute() else manifest_dir/path).resolve()


def identity(path):
    stat = path.stat()
    return dict(path=str(path), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)


def process_genome(row, manifest_dir, output_root, args):
    gid = row['short_name']
    fa = resolve_input(row['fa_path'], manifest_dir)
    gff = resolve_input(row['gff_path'], manifest_dir)
    inputs = dict(fasta=identity(fa), gff=identity(gff))
    parameters = dict(version=VERSION, genome_id=gid, seed=args.seed,
                      target_bp=args.target_bp, min_contig_bp=args.min_contig_bp,
                      fallback_min_contig_bp=args.fallback_min_contig_bp)
    target = output_root/gid
    filenames = dict(fa_path=f'{gid}.selected.fa.gz', gff_path=f'{gid}.selected.gff3.gz',
                     selected_contigs_path='selected_contigs.txt', contig_table_path='selected_contigs.tsv')
    if target.exists():
        if not args.resume:
            raise ValueError(f'Output exists: {target}; use a new outdir or --resume')
        with (target/'selection.json').open() as f:
            saved = json.load(f)
        previous=dict(saved['parameters'])
        previous.setdefault('fallback_min_contig_bp',[])
        current=dict(parameters)
        # A fallback list cannot affect a completed selection that succeeded at
        # the initial threshold. This also permits resuming pre-fallback outputs.
        effective=saved.get('effective_min_contig_bp',previous['min_contig_bp'])
        if effective == args.min_contig_bp:
            previous.pop('fallback_min_contig_bp',None)
            current.pop('fallback_min_contig_bp',None)
        if saved['inputs'] != inputs or previous != current:
            raise ValueError(f'{gid}: resume input/parameter mismatch')
        for name, size in saved['output_sizes'].items():
            if not (target/name).is_file() or (target/name).stat().st_size != size:
                raise ValueError(f'{gid}: missing/changed completed output {name}')
        saved.setdefault('effective_min_contig_bp',effective)
        saved.setdefault('fallback_used',False)
        print(f'  [RESUME] Verified completed output: {gid}',flush=True)
        return saved
    lengths = fasta_inventory(fa)
    attempts=[]
    for effective_min in [args.min_contig_bp]+args.fallback_min_contig_bp:
        chosen, eligible = select_contigs(lengths, gid, args.target_bp, effective_min, args.seed)
        attempts.append(dict(min_contig_bp=effective_min,eligible_contigs=len(eligible),
                             eligible_bp=sum(lengths[name] for name in eligible)))
        if chosen:
            break
    if not chosen:
        raise ValueError(f'{gid}: no sequence meets any minimum threshold {[a["min_contig_bp"] for a in attempts]}; longest={max(lengths.values()):,} bp')
    if effective_min != args.min_contig_bp:
        print(f'  [FALLBACK] {args.min_contig_bp:,} -> {effective_min:,} bp (zero eligible contigs at higher thresholds)',flush=True)
    actual = sum(lengths[name] for name in chosen)
    staging = Path(tempfile.mkdtemp(prefix=f'.{gid}.', dir=output_root))
    try:
        written = extract_fasta(fa, staging/filenames['fa_path'], chosen)
        if written != actual:
            raise ValueError('Extracted FASTA base count differs from inventory')
        gff_stats = extract_gff(gff, staging/filenames['gff_path'], chosen, lengths)
        with (staging/filenames['selected_contigs_path']).open('w') as f:
            f.write('\n'.join(chosen) + '\n')
        with (staging/filenames['contig_table_path']).open('w', newline='') as f:
            writer = csv.writer(f, delimiter='\t', lineterminator='\n')
            writer.writerow(['selection_rank', 'seqid', 'length_bp', 'cumulative_bp'])
            cumulative = 0
            for i, name in enumerate(chosen, 1):
                cumulative += lengths[name]
                writer.writerow([i, name, lengths[name], cumulative])
        stats = dict(parameters=parameters, inputs=inputs,
                     effective_min_contig_bp=effective_min,
                     fallback_used=effective_min!=args.min_contig_bp,threshold_attempts=attempts,
                     output_paths={k:str(target/v) for k,v in filenames.items()},
                     status='target_reached' if actual >= args.target_bp else 'insufficient_eligible_bp',
                     genome_bp=sum(lengths.values()), genome_contigs=len(lengths),
                     eligible_bp=sum(lengths[name] for name in eligible), eligible_contigs=len(eligible),
                     selected_bp=actual, selected_contigs=len(chosen),
                     selected_fraction=actual/sum(lengths.values()),
                     overshoot_bp=max(0, actual-args.target_bp),
                     max_selected_contig_bp=max(lengths[name] for name in chosen), **gff_stats)
        stats['output_sizes'] = {v:(staging/v).stat().st_size for v in filenames.values()}
        # Protect resumability if sources were changed during the two-pass extraction.
        if inputs != dict(fasta=identity(fa), gff=identity(gff)):
            raise ValueError('Input changed during extraction')
        with (staging/'selection.json').open('w') as f:
            json.dump(stats, f, indent=2)
        os.rename(staging, target)
        return stats
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--manifest', required=True, type=Path)
    parser.add_argument('--target-bp', '--size', type=int, default=50_000_000)
    parser.add_argument('--min-contig-bp', '--min-contig-size', type=int, default=100_000)
    parser.add_argument('--fallback-min-contig-bp',nargs='+',type=int,default=[],metavar='BP',
                        help='Descending fallback thresholds; used only if zero contigs meet the previous threshold')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--outdir', required=True, type=Path)
    parser.add_argument('--species', action='append', help='Only these short_name IDs; repeat for a pilot')
    parser.add_argument('--resume', action='store_true', help='Reuse completed outputs with matching source stats and parameters')
    parser.add_argument('--continue-on-error', action='store_true', help='Continue other genomes; exit nonzero if any failed')
    args = parser.parse_args()
    if args.target_bp <= 0 or args.min_contig_bp < 0:
        parser.error('target-bp must be positive; min-contig-bp must be nonnegative')
    thresholds=[args.min_contig_bp]+args.fallback_min_contig_bp
    if any(x<=0 for x in args.fallback_min_contig_bp) or any(a<=b for a,b in zip(thresholds,thresholds[1:])):
        parser.error('Fallback thresholds must be positive and strictly decreasing below min-contig-bp')
    manifest = args.manifest.expanduser().resolve()
    root = args.outdir.expanduser().resolve()
    with manifest.open(newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames or []
        required = {'fa_path','gff_path','short_name'}
        if not required <= set(fields):
            parser.error(f'Manifest requires columns {sorted(required)}')
        rows = list(reader)
    seen = set()
    for row in rows:
        gid = row['short_name']
        if not gid or gid in ('.','..') or '/' in gid or '\\' in gid or gid in seen:
            parser.error(f'Empty, unsafe or duplicate short_name: {gid!r}')
        seen.add(gid)
        if not row.get('fa_path') or not row.get('gff_path'):
            parser.error(f'{gid}: missing fa_path/gff_path')
    if args.species:
        missing = set(args.species)-seen
        if missing:
            parser.error(f'Unknown short_name: {sorted(missing)}')
        rows = [row for row in rows if row['short_name'] in args.species]
    if not rows:
        parser.error('No genomes to process')
    root.mkdir(parents=True, exist_ok=True)
    config = dict(version=VERSION, manifest=str(manifest), target_bp=args.target_bp,
                  min_contig_bp=args.min_contig_bp, seed=args.seed)
    config_path = root/'selection_config.json'
    if config_path.exists():
        with config_path.open() as f:
            if json.load(f) != config:
                parser.error('Output directory belongs to a different manifest/seed/size configuration; use a new outdir')
    else:
        with config_path.open('w') as f:
            json.dump(config,f,indent=2)
    summary_fields = ['short_name','status','target_bp','genome_bp','eligible_bp','selected_bp',
                      'selected_fraction','genome_contigs','eligible_contigs','selected_contigs',
                      'effective_min_contig_bp','fallback_used','overshoot_bp','max_selected_contig_bp','gff_selected_records','sequences_without_gff','error']
    extras = ['original_fa_path','original_gff_path','selected_contigs_path','selected_bp','selected_fraction','selection_status','effective_min_contig_bp','fallback_used']
    output_fields = fields + [x for x in extras if x not in fields]
    summaries, selected_rows, failed = [], [], False
    def save_tables():
        for name, columns, data in [('selection_summary.tsv',summary_fields,summaries),
                                    ('selected_manifest.tsv',output_fields,selected_rows)]:
            # A failed resume must not erase a previously successful manifest.
            if name == 'selected_manifest.tsv' and not data and (root/name).exists():
                continue
            tmp = root/(name+'.tmp')
            with tmp.open('w',newline='') as f:
                writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t', lineterminator='\n', extrasaction='ignore')
                writer.writeheader(); writer.writerows(data)
            os.replace(tmp,root/name)
    for row in rows:
        gid = row['short_name']
        print(f'[SELECT] {gid}', flush=True)
        try:
            stats = process_genome(row, manifest.parent, root, args)
            summaries.append(dict(short_name=gid, target_bp=args.target_bp, error='',
                                  sequences_without_gff=','.join(stats['selected_sequences_without_gff']), **stats))
            selected_rows.append({**row, 'original_fa_path':stats['inputs']['fasta']['path'],
                                  'original_gff_path':stats['inputs']['gff']['path'],
                                  'fa_path':stats['output_paths']['fa_path'], 'gff_path':stats['output_paths']['gff_path'],
                                  'selected_contigs_path':stats['output_paths']['selected_contigs_path'],
                                  'selected_bp':stats['selected_bp'], 'selected_fraction':stats['selected_fraction'],
                                  'effective_min_contig_bp':stats['effective_min_contig_bp'],
                                  'fallback_used':stats['fallback_used'],
                                  'selection_status':stats['status']})
            print(f"  selected {stats['selected_contigs']} contigs, {stats['selected_bp']:,} bp; {stats['status']}", flush=True)
            if stats['selected_sequences_without_gff']:
                print(f"[WARNING] {gid}: selected sequences without GFF features: {stats['selected_sequences_without_gff']}", file=sys.stderr)
        except (ValueError, OSError, KeyError) as exc:
            failed = True
            summaries.append(dict(short_name=gid, status='failed', target_bp=args.target_bp, error=str(exc)))
            print(f'[ERROR] {gid}: {exc}', file=sys.stderr)
            save_tables()
            if not args.continue_on_error:
                return 1
        save_tables()
    print(f'[DONE] {root / "selected_manifest.tsv"}')
    return int(failed)


if __name__ == '__main__':
    sys.exit(main())
