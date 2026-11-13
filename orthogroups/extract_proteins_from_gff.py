#!/usr/bin/env python3
"""Select the longest coding transcript per gene, then translate genomic CDS.

GFF3 + FASTA only; standard genetic code; Python standard library only.
Length means spliced exon length (CDS+UTR union if exons are absent).
Internal CDS phase bases complete split codons and must not be removed.
"""
import argparse
import csv
import gzip
import itertools
import json
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

CODONS = dict(zip((''.join(c) for c in itertools.product('TCAG', repeat=3)),
    'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
COMPLEMENT = str.maketrans('ACGTRYMKBDHVNSW', 'TGCAYRKMVHDBNSW')

def open_text(path):
    return gzip.open(path, 'rt', encoding='utf-8', errors='replace') if str(path).endswith('.gz') else open(path, encoding='utf-8', errors='replace')

def read_fasta(path):
    seqs={}; name=None; parts=[]
    with open_text(path) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if name is not None: seqs[name]=''.join(parts).upper()
                name=line[1:].split()[0]; parts=[]
                if name in seqs: raise ValueError(f'Duplicate FASTA ID: {name}')
            else:
                if name is None: raise ValueError('Sequence before first FASTA header')
                parts.append(''.join(line.split()))
    if name is not None: seqs[name]=''.join(parts).upper()
    if not seqs: raise ValueError('Empty FASTA')
    return seqs

def attrs(text):
    result={}
    if '=' not in text:
        if text=='.':return result
        raise ValueError('Expected GFF3 key=value attributes; GTF is not supported')
    for part in text.split(';'):
        if not part: continue
        # Some annotations use unescaped semicolons in free-text Note values.
        # These fragments carry no ID/Parent information.
        if '=' not in part:continue
        k,v=part.split('=',1); result[k]=[unquote(x) for x in v.split(',')]
    return result

def union_length(intervals):
    end=-1; length=0
    for a,b in sorted(set(intervals)):
        length+=max(0,b-max(a,end));end=max(end,b)
    return length

def read_gff(path):
    nodes={}; child=defaultdict(lambda:defaultdict(list)); genes=set()
    feature_types={'cds','exon','five_prime_utr','three_prime_utr','utr'}
    with open_text(path) as f:
        for line_no,line in enumerate(f,1):
            if line.startswith('##FASTA'): break
            if not line.strip() or line.startswith('#'): continue
            # Preserve occasional literal tabs inside free-text descriptions.
            cols=line.rstrip('\n').split('\t',8)
            if len(cols)!=9: raise ValueError(f'GFF3 line {line_no}: expected 9 columns')
            seqid,_,kind,start,end,_,strand,phase,attributes=cols
            start,end=int(start)-1,int(end)
            if start<0 or end<=start: raise ValueError(f'GFF3 line {line_no}: invalid coordinates')
            a=attrs(attributes);fid=a.get('ID',[None])[0];parents=a.get('Parent',[]);kind=kind.lower()
            if fid and kind in {'gene','pseudogene','transposable_element_gene','mrna','transcript','pseudogenic_transcript'}:
                if fid in nodes: raise ValueError(f'Duplicate gene/transcript ID: {fid}')
                nodes[fid]=(kind,parents,seqid,strand)
                if kind in {'gene','pseudogene','transposable_element_gene'}: genes.add(fid)
            if kind in feature_types:
                if not parents: raise ValueError(f'GFF3 line {line_no}: {kind} lacks Parent')
                for parent in parents:child[parent][kind].append((seqid,start,end,strand,phase))
    transcripts=defaultdict(list)
    for tid,features in child.items():
        if not features['cds']: continue
        if tid not in nodes: raise ValueError(f'CDS parent {tid} lacks a transcript feature')
        kind,parents,seqid,strand=nodes[tid]
        if kind in {'gene','pseudogene','transposable_element_gene'}: raise ValueError(f'CDS directly under gene {tid}; explicit transcripts are required')
        if len(parents)!=1 or parents[0] not in genes: raise ValueError(f'Transcript {tid} requires one annotated gene Parent')
        all_features=[v for vs in features.values() for v in vs]
        if strand not in {'+','-'} or any(v[0]!=seqid or v[3]!=strand for v in all_features):
            raise ValueError(f'Inconsistent contig/strand in {tid}')
        exons=features['exon'];basis='exon_union'
        if not exons:exons=all_features;basis='CDS_UTR_union'
        transcripts[parents[0]].append(dict(id=tid,seqid=seqid,strand=strand,features=features,
            transcript_length=union_length((v[1],v[2]) for v in exons),
            cds_length=union_length((v[1],v[2]) for v in features['cds']),length_basis=basis))
    return genes,transcripts

def translate_cds(tx,genome):
    seqid=tx['seqid'];strand=tx['strand']
    if seqid not in genome: raise ValueError('missing_FASTA_contig')
    seq=genome[seqid]
    parts=sorted(set(tx['features']['cds']),key=lambda v:v[1],reverse=strand=='-')
    genomic=sorted(parts,key=lambda v:v[1])
    if any(b[1]<a[2] for a,b in zip(genomic,genomic[1:])): raise ValueError('overlapping_CDS')
    chunks=[];warnings=[];length=0;first_phase=0
    for i,(_,start,end,_,phase) in enumerate(parts):
        if end>len(seq): raise ValueError('CDS_outside_FASTA')
        if phase not in {'0','1','2'}: raise ValueError('missing_or_invalid_CDS_phase')
        phase=int(phase);piece=seq[start:end]
        if strand=='-': piece=piece.translate(COMPLEMENT)[::-1]
        if i==0:
            first_phase=phase;piece=piece[phase:]
            if phase:warnings.append('partial_5prime')
        elif phase!=(-length)%3: raise ValueError('inconsistent_CDS_phase')
        chunks.append(piece);length+=len(piece)
    dna=''.join(chunks);remainder=len(dna)%3
    if remainder:warnings.append('partial_3prime');dna=dna[:-remainder]
    if not dna: raise ValueError('empty_CDS')
    if set(dna)-set('ACGTRYMKBDHVNSW'): raise ValueError('invalid_DNA_characters')
    aa=''.join(CODONS.get(dna[i:i+3],'X') for i in range(0,len(dna),3))
    if aa.endswith('*'):aa=aa[:-1]
    if '*' in aa:raise ValueError('internal_stop_codon')
    return aa,warnings,first_phase,remainder

def extract(args):
    genome=read_fasta(args.fasta);genes,txs=read_gff(args.gff)
    output=Path(args.output);output.parent.mkdir(parents=True,exist_ok=True)
    report=Path(args.report or str(output)+'.selection.tsv');report.parent.mkdir(parents=True,exist_ok=True)
    rows=[];records=[];output_ids=set()
    for gene in sorted(genes):
        row=dict(gene_id=gene,transcript_id='',protein_id='',transcript_length='',cds_length='',length_basis='',protein_length='',status='',detail='')
        candidates=txs.get(gene,[])
        if not candidates:row['status']='no_coding_transcript';rows.append(row);continue
        key='transcript_length' if args.length_metric=='transcript' else 'cds_length'
        tx=sorted(candidates,key=lambda t:(-t[key],-t['cds_length'],t['id']))[0]
        row.update(transcript_id=tx['id'],transcript_length=tx['transcript_length'],cds_length=tx['cds_length'],length_basis=tx['length_basis'])
        protein_id=gene+args.id_suffix
        if '\ufffd' in protein_id or any(c.isspace() for c in protein_id) or protein_id in output_ids:raise ValueError('Invalid/duplicate output identifier')
        output_ids.add(protein_id);row['protein_id']=protein_id
        try:
            aa,warnings,_,_=translate_cds(tx,genome);row['protein_length']=len(aa);row['detail']=';'.join(warnings)
            nterm=len(aa)-len(aa.lstrip('X'))
            if len(aa)<args.min_len:row['status']='short_protein'
            elif args.nterm_x>0 and nterm>=args.nterm_x:row['status']='Nterminal_X_run'
            elif aa.count('X')/max(1,len(aa))>args.max_xfrac:row['status']='high_X_fraction'
            else:row['status']='kept';records.append((protein_id,aa))
        except ValueError as e:row['status']='translation_rejected';row['detail']=str(e)
        rows.append(row)
    with output.open('w') as f:
        for name,aa in records:
            f.write('>'+name+'\n')
            for i in range(0,len(aa),80):f.write(aa[i:i+80]+'\n')
    fields=['gene_id','transcript_id','protein_id','transcript_length','cds_length','length_basis','protein_length','status','detail']
    with report.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
    counts=defaultdict(int)
    for row in rows:counts[row['status']]+=1
    print(json.dumps({'genes':len(genes),'proteins':len(records),'status_counts':dict(counts),'genetic_code':1,'length_metric':args.length_metric},indent=2))
    return len(records)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta',required=True,help='Genome FASTA or FASTA.gz')
    p.add_argument('--gff',required=True,help='GFF3 or GFF3.gz; explicit gene/transcript/CDS hierarchy')
    p.add_argument('-o','--output',required=True)
    p.add_argument('--report',help='Selection/QC TSV; default OUTPUT.selection.tsv')
    p.add_argument('--length-metric',choices=['transcript','cds'],default='transcript')
    p.add_argument('--id-suffix',default='',help='Optional suffix, e.g. |Mpo, matching OG identifiers')
    p.add_argument('--min_len',type=int,default=50)
    p.add_argument('--max_xfrac',type=float,default=.30)
    p.add_argument('--nterm_x',type=int,default=20,help='Reject this many initial X; 0 disables')
    a=p.parse_args()
    if a.min_len<0 or a.nterm_x<0 or not 0<=a.max_xfrac<=1:p.error('Invalid QC thresholds')
    if extract(a)==0:raise SystemExit('No proteins exported; inspect selection report')

if __name__=='__main__':main()
