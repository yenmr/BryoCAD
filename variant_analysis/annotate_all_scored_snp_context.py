#!/usr/bin/env python3
import argparse, gzip, re
from collections import defaultdict
import pandas as pd
import pysam
from intervaltree import IntervalTree
from Bio import SeqIO
from tqdm import tqdm

def open_text(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "rt")

def clean_gene_id(x):
    s = str(x).strip()
    if s == "" or s.lower() in {"nan", "none", "."}:
        return ""
    return s.split()[0].split("|")[0]

def parse_attrs(s):
    d = {}
    for x in s.strip().split(";"):
        if not x:
            continue
        if "=" in x:
            k, v = x.split("=", 1)
            d[k] = v
    return d

def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()

CODON = {
"TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
"TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
"CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
"CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
"ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
"AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
"GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
"GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}

def load_genome(fasta):
    genome = {}
    with open_text(fasta) as f:
        for rec in SeqIO.parse(f, "fasta"):
            genome[rec.id] = str(rec.seq).upper()
    return genome

def load_og(path, species):
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    df = df[df["Species"] == species].copy()
    df["gene_clean"] = df["Gene"].map(clean_gene_id)
    return dict(zip(df["gene_clean"], df["OG_category"]))

def build_gff(gff, genome):
    trees = defaultdict(IntervalTree)
    tx_to_gene = {}
    cds_by_tx = defaultdict(list)

    with open_text(gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue

            chrom, _, feat, start, end, _, strand, _, attrs = p
            if chrom not in genome:
                continue

            start0, end0 = int(start)-1, int(end)
            a = parse_attrs(attrs)
            fid = a.get("ID", "")
            parent = a.get("Parent", "")

            if feat == "gene":
                gid = clean_gene_id(fid or a.get("Name", ""))
                trees[chrom].addi(start0, end0, {"type":"gene","gene_id":gid,"tx_id":"","strand":strand})

                if strand == "+":
                    ps, pe = max(0, start0-1000), start0
                else:
                    ps, pe = end0, min(len(genome[chrom]), end0+1000)
                if pe > ps:
                    trees[chrom].addi(ps, pe, {"type":"promoter_1kb","gene_id":gid,"tx_id":"","strand":strand})

            elif feat in {"mRNA", "transcript"}:
                if fid and parent:
                    tx_to_gene[fid] = clean_gene_id(parent.split(",")[0])

            elif feat == "CDS":
                for tx in parent.split(","):
                    gid = clean_gene_id(tx_to_gene.get(tx, tx))
                    trees[chrom].addi(start0, end0, {"type":"CDS","gene_id":gid,"tx_id":tx,"strand":strand})
                    cds_by_tx[tx].append((chrom,start0,end0,strand,gid))

            elif feat in {"five_prime_UTR", "three_prime_UTR", "UTR"}:
                for tx in parent.split(","):
                    gid = clean_gene_id(tx_to_gene.get(tx, tx))
                    trees[chrom].addi(start0, end0, {"type":"UTR","gene_id":gid,"tx_id":tx,"strand":strand})

            elif feat == "exon":
                for tx in parent.split(","):
                    gid = clean_gene_id(tx_to_gene.get(tx, tx))
                    chrom_len = len(genome[chrom])
                    if strand == "+":
                        intervals = [
                            (max(0,start0-2), start0, "splice_acceptor"),
                            (end0, min(chrom_len,end0+2), "splice_donor")
                        ]
                    else:
                        intervals = [
                            (max(0,start0-2), start0, "splice_donor"),
                            (end0, min(chrom_len,end0+2), "splice_acceptor")
                        ]
                    for s,e,t in intervals:
                        if e > s:
                            trees[chrom].addi(s, e, {"type":t,"gene_id":gid,"tx_id":tx,"strand":strand})
    return trees, cds_by_tx

def cds_effect(chrom,pos0,ref,alt,genome,cds_by_tx,tx):
    cds = cds_by_tx.get(tx, [])
    if not cds:
        return "CDS_unknown","",""
    strand = cds[0][3]
    gid = cds[0][4]
    ordered = sorted(cds, key=lambda x:x[1], reverse=(strand=="-"))

    seqs, gpos = [], []
    for c,s,e,st,g in ordered:
        seq = genome[c][s:e]
        pos = list(range(s,e))
        if strand == "-":
            seq = revcomp(seq)
            pos = pos[::-1]
        seqs.append(seq); gpos.extend(pos)

    cds_seq = "".join(seqs)
    if pos0 not in gpos:
        return "CDS_unknown",gid,tx

    idx = gpos.index(pos0)
    codon = cds_seq[(idx//3)*3:(idx//3)*3+3]
    if len(codon) != 3 or any(b not in "ACGT" for b in codon):
        return "CDS_unknown",gid,tx

    alt_base = revcomp(alt) if strand == "-" else alt
    alt_codon = list(codon)
    alt_codon[idx % 3] = alt_base
    alt_codon = "".join(alt_codon)

    aa1, aa2 = CODON.get(codon,"X"), CODON.get(alt_codon,"X")

    if idx < 3 and aa1 == "M" and aa2 != "M":
        ctx = "start_lost_CDS"
    elif aa1 == aa2:
        ctx = "synonymous_CDS"
    elif aa2 == "*":
        ctx = "stop_gained_CDS"
    elif aa1 == "*":
        ctx = "stop_lost_CDS"
    else:
        ctx = "missense_CDS"

    return ctx,gid,tx

def choose_context(overlaps, chrom,pos0,ref,alt,genome,cds_by_tx):
    if not overlaps:
        return "distal_intergenic","",""

    priority = {"splice_donor":0,"splice_acceptor":1,"CDS":2,"UTR":3,"promoter_1kb":4,"gene":5}
    best = sorted(overlaps, key=lambda iv: priority.get(iv.data["type"],99))[0].data

    if best["type"] == "CDS":
        return cds_effect(chrom,pos0,ref,alt,genome,cds_by_tx,best["tx_id"])
    if best["type"] == "gene":
        return "intron_or_gene_body", best["gene_id"], best["tx_id"]
    return best["type"], best["gene_id"], best["tx_id"]

def haploid_af(rec):
    alt_count, n_called = 0, 0
    for s in rec.samples.values():
        gt = s.get("GT")
        if gt is None:
            continue
        for a in gt:
            if a is None:
                continue
            n_called += 1
            if a == 1:
                alt_count += 1
    if n_called == 0:
        return None,None,None,None
    af = alt_count / n_called
    maf = min(af, 1-af)
    return alt_count,n_called,af,maf

def get_score(info):
    m = re.search(r"plantCAD_zero_shot=([^;]+)", str(info))
    return None if m is None else float(m.group(1).split(",")[0])

def fig4d_group(ctx):
    if ctx in {"missense_CDS","stop_gained_CDS","stop_lost_CDS","start_lost_CDS"}:
        return "Missense"
    if ctx == "synonymous_CDS":
        return "Synonymous"
    return "Noncoding"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--og_table", required=True)
    ap.add_argument("--og_species", default="M_polymorpha_v61")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    genome = load_genome(args.fasta)
    og = load_og(args.og_table, args.og_species)
    trees, cds_by_tx = build_gff(args.gff, genome)
    vcf = pysam.VariantFile(args.vcf)

    rows = []
    for rec in tqdm(vcf.fetch(), desc="Annotating scored SNPs"):
        if rec.alts is None or len(rec.alts) != 1:
            continue
        ref, alt = rec.ref.upper(), rec.alts[0].upper()
        if ref not in "ACGT" or alt not in "ACGT":
            continue

        alt_count,n_called,af,maf = haploid_af(rec)
        if maf is None:
            continue

        chrom, pos0 = rec.chrom, rec.pos - 1
        overlaps = trees[chrom].at(pos0) if chrom in trees else []
        ctx,gid,tx = choose_context(overlaps, chrom,pos0,ref,alt,genome,cds_by_tx)
        gclean = clean_gene_id(gid)
        gclass = og.get(gclean, "none" if ctx == "distal_intergenic" else "unknown")
        score = rec.info.get("plantCAD_zero_shot", None)
        if isinstance(score, tuple):
            score = score[0]

        rows.append({
            "chrom": chrom, "pos": rec.pos, "ref": ref, "alt": alt,
            "alt_count": alt_count, "n_called": n_called,
            "alt_frequency": af, "maf": maf,
            "zero_shot_score": score,
            "deleterious_score": -float(score),
            "context": ctx,
            "fig4d_group": fig4d_group(ctx),
            "gene_id": gid,
            "gene_class": gclass
        })

    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote: {args.out}")
    print(df["context"].value_counts())
    print(df["fig4d_group"].value_counts())

if __name__ == "__main__":
    main()
