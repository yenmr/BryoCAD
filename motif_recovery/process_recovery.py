#!/usr/bin/env python3
"""BryoCAD recovery summaries, Python >=3.9, standard library only.
Rows are retained as supplied (no deduplication). Invalid metric rows are excluded
and counted. Only reference_is_canonical=True rows enter performance summaries.
"""
import argparse
import csv
import gzip
import json
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path


def reader(path):
    op = gzip.open if str(path).endswith('.gz') else open
    return op(path, 'rt', encoding='utf-8-sig', newline='')


def orthogroup_gene_id(gene, species):
    """Strip only a verified species suffix; allow hyphen/underscore variation.

    Do not alter the gene ID or strip transcript suffixes such as .t1.
    Unrelated pipe suffixes are retained rather than creating a false match.
    """
    gene = gene.strip()
    if '|' not in gene:
        return gene, 'no_suffix'
    bare, suffix = gene.rsplit('|', 1)
    if suffix == species:
        return bare, 'exact_suffix'
    if suffix.replace('-', '_') == species.replace('-', '_'):
        return bare, 'normalized_species_suffix'
    return gene, 'unrecognized_suffix'


def boolean(value):
    value = str(value).strip().lower()
    if value in ('true', '1'): return True
    if value in ('false', '0'): return False
    raise ValueError('Invalid boolean: ' + repr(value))


def table(path, key):
    with reader(path) as handle:
        records = list(csv.DictReader(handle, delimiter='\t'))
    result = {}
    for row in records:
        if not row.get(key) or row[key] in result:
            raise ValueError(f'{path}: missing/duplicate {key}: {row.get(key)}')
        result[row[key]] = row
    return result


def writer(folder, name, fields):
    handle = (folder / name).open('w', newline='')
    obj = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
    obj.writeheader(); handle.flush()
    return handle, obj


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results', required=True, type=Path)
    p.add_argument('--samples', required=True, type=Path)
    p.add_argument('--genomes', required=True, type=Path)
    p.add_argument('--orthogroups', required=True, type=Path)
    p.add_argument('--outdir', required=True, type=Path)
    a = p.parse_args(); a.outdir.mkdir(parents=True, exist_ok=True)
    samples = table(a.samples, 'short_name'); genomes = table(a.genomes, 'short_name')
    for name, row in samples.items():
        for field in ('species', 'division', 'order', 'train'):
            if not row.get(field): raise ValueError(f'{name}: missing {field}')
        row['_train'] = boolean(row['train'])
    seen_orders = {r['order'] for r in samples.values() if r['_train']}
    # Disk-backed index keeps the large long-format orthogroup table out of RAM.
    db = sqlite3.connect(a.outdir / 'orthogroup_index.sqlite')
    db.execute('DROP TABLE IF EXISTS mapping')
    db.execute('CREATE TABLE mapping (species TEXT, gene TEXT, category TEXT, PRIMARY KEY(species,gene))')
    batch = []; category_counts = Counter(); suffix_counts = Counter()
    print('Indexing orthogroups...', file=sys.stderr, flush=True)
    with reader(a.orthogroups) as handle:
        rows = csv.DictReader(handle, delimiter='\t')
        if not {'Species','Gene','OG_category'} <= set(rows.fieldnames or []):
            raise ValueError('Orthogroup file needs Species, Gene, OG_category')
        for row in rows:
            sp = row['Species']
            if sp not in samples: continue
            gene, suffix_mode = orthogroup_gene_id(row['Gene'], sp)
            suffix_counts[(sp, suffix_mode)] += 1
            cat = row['OG_category'].strip() or 'unclassified'
            if cat in ('unmapped','ambiguous'): raise ValueError('Reserved OG_category: '+cat)
            category_counts[cat] += 1
            batch.append((sp, gene, cat))
            if len(batch) >= 50000:
                insert(db, batch); batch.clear()
        insert(db, batch); db.commit()
    sh, sw = writer(a.outdir, 'orthogroup_id_qc.tsv', ['short_name', 'suffix_handling', 'n_rows'])
    for (sp, mode), n in sorted(suffix_counts.items()):
        sw.writerow(dict(short_name=sp, suffix_handling=mode, n_rows=n))
        if mode in ('normalized_species_suffix', 'unrecognized_suffix'):
            print(f'OG ID QC {sp}: {mode}={n}', file=sys.stderr, flush=True)
    sh.close()
    gh, gw = writer(a.outdir, 'gene_category_counts.tsv', ['short_name', 'gene_category', 'n_genes'])
    for sp, cat, n in db.execute('SELECT species,category,COUNT(*) FROM mapping GROUP BY species,category ORDER BY species,category'):
        gw.writerow(dict(short_name=sp, gene_category=cat, n_genes=n))
    gh.close()
    print('Observed OG categories: '+json.dumps(category_counts), file=sys.stderr, flush=True)
    meta_fields = ['short_name','species','division','order','training_group']
    metric_fields = ['motif','n_sites','n_genes','n_exact','exact_recovery_rate','n_predicted_canonical','canonical_recovery_rate']
    outputs = {}
    outputs['main'] = writer(a.outdir,'species_motif_summary.tsv',meta_fields+metric_fields)
    outputs['cat'] = writer(a.outdir,'species_gene_category_summary.tsv',meta_fields+['gene_category']+metric_fields)
    outputs['tts'] = writer(a.outdir,'tts_codon_confusion.tsv',meta_fields+['reference_motif','predicted_motif','n_sites'])
    outputs['qc'] = writer(a.outdir,'analysis_qc.tsv',['short_name','motif','status','file','n_rows','n_noncanonical_reference','n_invalid_rows','n_metric_mismatches','n_duplicate_sites','n_retained','n_unmapped','n_ambiguous','n_mapped_genes','n_retained_genes','gene_mapping_rate','detail'])
    # Separate metadata prevents silent collisions with the source's fields.
    metadata = []
    problems = 0
    try:
        for i, (sp, sample) in enumerate(samples.items(),1):
            group = 'trained' if sample['_train'] else 'unseen_species_seen_order' if sample['order'] in seen_orders else 'unseen_order'
            meta = {k:sample[k] for k in ('species','division','order')}
            meta.update(short_name=sp,training_group=group)
            genome = genomes.get(sp, {})
            metadata.append(dict(meta, **{'genome_'+k:v for k,v in genome.items() if k != 'short_name'}))
            og = dict(db.execute('SELECT gene,category FROM mapping WHERE species=?',(sp,)))
            for motif in ('acceptor','donor','tis','tts'):
                files = sorted((a.results/sp).glob('*_'+motif+'_recovery.tsv*'))
                # Also allow bare motif filenames.
                files = sorted(set(files + list((a.results/sp).glob(motif+'_recovery.tsv*'))))
                files = [f for f in files if f.name.endswith(('.tsv','.tsv.gz'))]
                qc = Counter(); details = []; site_seen = set(); genes = set(); mapped_genes = set()
                total = [0,0,0]; bycat = defaultdict(lambda:[0,0,0]); catgenes = defaultdict(set); confusion = Counter()
                status = 'ok'
                if len(files) != 1:
                    status = 'error'; details.append(f'Expected one file, found {len(files)}')
                else:
                    with reader(files[0]) as handle:
                        rows = csv.DictReader(handle, delimiter='\t')
                        required = {'geneid','chromosome','start','end','direction','reference_motif','predicted_motif','reference_is_canonical','exact_recovery','predicted_is_canonical'}
                        missing = required-set(rows.fieldnames or [])
                        if missing:
                            status='error'; details.append('Missing columns: '+','.join(sorted(missing)))
                        else:
                            for line, row in enumerate(rows,2):
                                qc['n_rows'] += 1
                                try:
                                    if None in row or any(row.get(k) is None for k in required): raise ValueError('Malformed TSV row')
                                    canonical = boolean(row['reference_is_canonical'])
                                    if not canonical:
                                        qc['n_noncanonical_reference'] += 1
                                        continue
                                    ref = row['reference_motif'].strip().upper(); pred = row['predicted_motif'].strip().upper()
                                    exact = boolean(row['exact_recovery']); pc = boolean(row['predicted_is_canonical'])
                                    expected_len = 3 if motif in ('tis','tts') else 2
                                    if len(ref)!=expected_len or len(pred)!=expected_len or set(ref+pred)-set('ACGT'): raise ValueError('Invalid motif length/bases')
                                    if not row['geneid']: raise ValueError('Empty geneid')
                                    start,end = int(row['start']),int(row['end'])
                                    if start<1 or end<start or row['direction'] not in ('+','-'): raise ValueError('Invalid coordinate/strand')
                                    if exact != (ref==pred):
                                        qc['n_metric_mismatches'] += 1
                                        raise ValueError('exact_recovery disagrees with motif strings')
                                    if motif=='tts' and (canonical != (ref in ('TAA','TAG','TGA')) or pc != (pred in ('TAA','TAG','TGA'))):
                                        qc['n_metric_mismatches'] += 1
                                        raise ValueError('TTS canonical flag disagrees with stop codon set')
                                except (ValueError,TypeError) as exc:
                                    qc['n_invalid_rows'] += 1
                                    if len(details)<5: details.append(f'line {line}: {exc}')
                                    continue
                                key = (row['chromosome'],start,end,row['direction'])
                                if key in site_seen: qc['n_duplicate_sites'] += 1
                                site_seen.add(key)
                                gene = row['geneid']; cat = og.get(gene,'unmapped')
                                genes.add(gene); catgenes[cat].add(gene)
                                if cat not in ('unmapped','ambiguous'): mapped_genes.add(gene)
                                if cat in ('unmapped','ambiguous'): qc['n_'+cat] += 1
                                for stats in (total,bycat[cat]):
                                    stats[0]+=1;stats[1]+=int(exact);stats[2]+=int(pc)
                                if motif=='tts': confusion[(ref,pred)]+=1
                if total[0]:
                    outputs['main'][1].writerow(dict(meta,**metrics(motif,total,len(genes))))
                    for cat, stats in sorted(bycat.items()):
                        outputs['cat'][1].writerow(dict(meta,gene_category=cat,**metrics(motif,stats,len(catgenes[cat]))))
                    for (ref,pred),n in sorted(confusion.items()):
                        outputs['tts'][1].writerow(dict(meta,reference_motif=ref,predicted_motif=pred,n_sites=n))
                elif status!='error':
                    status='error';details.append('No valid canonical-reference sites')
                if status=='ok' and any(qc[k] for k in ('n_invalid_rows','n_duplicate_sites','n_ambiguous')): status='warning'
                if status!='ok': problems+=1
                qrow = dict(qc,short_name=sp,motif=motif,status=status,file=str(files[0]) if len(files)==1 else '',n_retained=total[0],n_retained_genes=len(genes),n_mapped_genes=len(mapped_genes),gene_mapping_rate=len(mapped_genes)/len(genes) if genes else 'NA',detail='; '.join(details))
                outputs['qc'][1].writerow({k:qrow.get(k,0) for k in outputs['qc'][1].fieldnames})
                for h,_ in outputs.values(): h.flush()
                print(f'[{i}/{len(samples)}] {sp} {motif}: {status}; retained={total[0]}; unmapped={qc["n_unmapped"]}; '+ '; '.join(details),file=sys.stderr,flush=True)
    finally:
        for h,_ in outputs.values(): h.close()
        db.close()
    fields = meta_fields+sorted(set().union(*(r.keys() for r in metadata))-set(meta_fields))
    h,w=writer(a.outdir,'species_metadata.tsv',fields)
    for r in metadata:w.writerow(r)
    h.close()
    (a.outdir/'run_config.json').write_text(json.dumps({'arguments':{k:str(v) for k,v in vars(a).items()},'canonical_filter':'reference_is_canonical','deduplicated':False,'OG_category_counts':category_counts,'problem_species_motif_pairs':problems},indent=2))
    print(f'Done. {problems} species/motif pairs require QC review.',file=sys.stderr)


def insert(db, batch):
    db.executemany('INSERT INTO mapping VALUES (?,?,?) ON CONFLICT(species,gene) DO UPDATE SET category=CASE WHEN category=excluded.category THEN category ELSE "ambiguous" END',batch)


def metrics(motif,stats,genes):
    n,exact,canonical=stats
    return dict(motif=motif,n_sites=n,n_genes=genes,n_exact=exact,exact_recovery_rate=exact/n,n_predicted_canonical=canonical if motif=='tts' else 'NA',canonical_recovery_rate=canonical/n if motif=='tts' else 'NA')


if __name__=='__main__': main()
