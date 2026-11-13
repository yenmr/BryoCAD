#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, gzip, random
from collections import defaultdict
import pandas as pd
import pysam
from Bio import SeqIO
from intervaltree import IntervalTree
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
        elif " " in x:
            k, v = x.split(" ", 1)
            d[k] = v.strip('"')
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

def load_core_genes(og_table, species):
    df = pd.read_csv(og_table, sep="\t", dtype=str, low_memory=False)
    df = df[(df["Species"] == species) & (df["OG_category"] == "core")].copy()
    genes = set(df["Gene"].map(clean_gene_id))
    print(f"[INFO] core genes: {len(genes):,}")
    return genes

def build_features(gff, genome, core_genes):
    tx_to_gene = {}
    genes = {}
    cds_by_tx = defaultdict(list)
    trees = defaultdict(IntervalTree)

    raw = []

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

            start0, end0 = int(start) - 1, int(end)
            a = parse_attrs(attrs)
            fid = a.get("ID", "")
            parent = a.get("Parent", "")

            if feat == "gene":
                gid = clean_gene_id(fid or a.get("Name", ""))
                if gid in core_genes:
                    genes[gid] = (chrom, start0, end0, strand)
                    raw.append((chrom, start0, end0, "gene", gid, "", strand))

            elif feat in {"mRNA", "transcript"}:
                if fid and parent:
                    tx_to_gene[fid] = clean_gene_id(parent.split(",")[0])

            elif feat == "CDS":
                for tx in parent.split(","):
                    gid = clean_gene_id(tx_to_gene.get(tx, tx))
                    if gid in core_genes:
                        raw.append((chrom, start0, end0, "CDS", gid, tx, strand))
                        cds_by_tx[tx].append((chrom, start0, end0, strand, gid))

            elif feat in {"five_prime_UTR", "three_prime_UTR", "UTR"}:
                for tx in parent.split(","):
                    gid = clean_gene_id(tx_to_gene.get(tx, tx))
                    if gid in core_genes:
                        raw.append((chrom, start0, end0, "UTR", gid, tx, strand))

            elif feat == "exon":
                for tx in parent.split(","):
                    gid = clean_gene_id(tx_to_gene.get(tx, tx))
                    if gid not in core_genes:
                        continue

                    chrom_len = len(genome[chrom])

                    splice_intervals = []

                    if strand == "+":
                        splice_intervals.append((max(0, start0 - 2), start0, "splice_acceptor"))
                        splice_intervals.append((end0, min(chrom_len, end0 + 2), "splice_donor"))
                    else:
                        splice_intervals.append((max(0, start0 - 2), start0, "splice_donor"))
                        splice_intervals.append((end0, min(chrom_len, end0 + 2), "splice_acceptor"))

                    for s, e, splice_type in splice_intervals:
                        if e <= s:
                            continue
                        raw.append((chrom, s, e, splice_type, gid, tx, strand))

    for gid, (chrom, s, e, strand) in genes.items():
        if strand == "+":
            ps, pe = max(0, s - 1000), s
        else:
            ps, pe = e, min(len(genome[chrom]), e + 1000)
        if pe > ps:
            raw.append((chrom, ps, pe, "promoter_1kb", gid, "", strand))

    for chrom, s, e, typ, gid, tx, strand in raw:
        trees[chrom].addi(s, e, {
            "type": typ,
            "gene_id": gid,
            "tx_id": tx,
            "strand": strand
        })

    print(f"[INFO] raw feature intervals: {len(raw):,}")
    return trees, cds_by_tx

def choose_best(overlaps):
    priority = {
        "splice_donor": 0,
        "splice_acceptor": 1,
        "CDS": 2,
        "UTR": 3,
        "promoter_1kb": 4,
        "gene": 5,
    }
    if not overlaps:
        return None
    return sorted(overlaps, key=lambda iv: priority.get(iv.data["type"], 99))[0].data

def cds_effect(chrom, pos0, ref, alt, genome, cds_by_tx, tx):
    cds = cds_by_tx.get(tx, [])
    if not cds:
        return "CDS_unknown", "", tx, None, None

    strand = cds[0][3]
    gid = cds[0][4]
    ordered = sorted(cds, key=lambda x: x[1], reverse=(strand == "-"))

    cds_seq = []
    gpos = []

    for c, s, e, st, g in ordered:
        seq = genome[c][s:e]
        pos = list(range(s, e))
        if strand == "-":
            seq = revcomp(seq)
            pos = pos[::-1]
        cds_seq.append(seq)
        gpos.extend(pos)

    cds_seq = "".join(cds_seq)

    if pos0 not in gpos:
        return "CDS_unknown", gid, tx, None, None

    idx = gpos.index(pos0)
    is_start_codon = idx < 3
    codon_start = (idx // 3) * 3
    codon = cds_seq[codon_start:codon_start + 3]

    if len(codon) != 3 or any(b not in "ACGT" for b in codon):
        return "CDS_unknown", gid, tx, None, None

    alt_base = revcomp(alt) if strand == "-" else alt
    alt_codon = list(codon)
    alt_codon[idx % 3] = alt_base
    alt_codon = "".join(alt_codon)

    aa_ref = CODON.get(codon, "X")
    aa_alt = CODON.get(alt_codon, "X")

    if is_start_codon and aa_ref == "M" and aa_alt != "M":
        context = "start_lost_CDS"
    elif aa_ref == aa_alt:
        context = "synonymous_CDS"
    elif aa_alt == "*":
        context = "stop_gained_CDS"
    elif aa_ref == "*":
        context = "stop_lost_CDS"
    else:
        context = "missense_CDS"

    return context, gid, tx, idx + 1, (idx % 3) + 1

def reservoir_add(bucket, seen, key, item, max_n, rng):
    seen[key] += 1
    if len(bucket[key]) < max_n:
        bucket[key].append(item)
    else:
        j = rng.randrange(seen[key])
        if j < max_n:
            bucket[key][j] = item

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--og_table", required=True)
    ap.add_argument("--og_species", default="M_polymorpha_v61")
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--max_per_category", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    nts = ["A", "C", "G", "T"]

    genome = load_genome(args.fasta)
    core_genes = load_core_genes(args.og_table, args.og_species)
    trees, cds_by_tx = build_features(args.gff, genome, core_genes)

    target_contexts = {
        "stop_gained_CDS",
        "stop_lost_CDS",
        "start_lost_CDS",
        "splice_donor",
        "splice_acceptor",
        "missense_CDS",
        "synonymous_CDS",
        "UTR",
        "promoter_1kb",
        "intron_or_gene_body",
    }

    bucket = defaultdict(list)
    seen = defaultdict(int)

    # Build unique positions from intervals, but process each position once.
    positions_by_chrom = defaultdict(set)
    for chrom, tree in trees.items():
        for iv in tree:
            if iv.data["type"] in {"CDS", "UTR", "promoter_1kb", "splice_20bp", "gene"}:
                positions_by_chrom[chrom].update(range(iv.begin, iv.end))

    total_positions = sum(len(x) for x in positions_by_chrom.values())
    print(f"[INFO] unique core-feature positions to scan: {total_positions:,}")

    for chrom in tqdm(sorted(positions_by_chrom.keys()), desc="Scanning chromosomes"):
        seq = genome[chrom]
        for pos0 in tqdm(sorted(positions_by_chrom[chrom]), desc=f"{chrom}", leave=False):
            ref = seq[pos0]
            if ref not in nts:
                continue

            best = choose_best(trees[chrom].at(pos0))
            if best is None:
                continue

            for alt in nts:
                if alt == ref:
                    continue

                if best["type"] == "CDS":
                    context, gid, tx, cds_pos, codon_pos = cds_effect(
                        chrom, pos0, ref, alt, genome, cds_by_tx, best["tx_id"]
                    )
                elif best["type"] == "gene":
                    context, gid, tx, cds_pos, codon_pos = (
                        "intron_or_gene_body", best["gene_id"], best["tx_id"], None, None
                    )
                else:
                    context, gid, tx, cds_pos, codon_pos = (
                        best["type"], best["gene_id"], best["tx_id"], None, None
                    )

                if context not in target_contexts:
                    continue

                row = {
                    "chrom": chrom,
                    "pos": pos0 + 1,
                    "id": ".",
                    "ref": ref,
                    "alt": alt,
                    "context": context,
                    "gene_id": gid,
                    "transcript_id": tx,
                    "cds_position": cds_pos,
                    "codon_position": codon_pos,
                    "variant_id": f"{chrom}:{pos0+1}:{ref}:{alt}"
                }

                reservoir_add(bucket, seen, context, row, args.max_per_category, rng)

    rows = []
    for k in sorted(bucket.keys()):
        rows.extend(bucket[k])

    df = pd.DataFrame(rows)
    contig_order = {c: i for i, c in enumerate(genome.keys())}
    df["contig_order"] = df["chrom"].map(contig_order).fillna(10**9).astype(int)
    df = df.sort_values(["contig_order", "pos", "ref", "alt"]).drop(columns=["contig_order"])

    out_tsv = args.out_prefix + ".annotation.tsv"
    out_vcf = args.out_prefix + ".vcf.gz"

    df.to_csv(out_tsv, sep="\t", index=False)

    header = pysam.VariantHeader()
    header.add_meta("fileformat", value="VCFv4.2")
    for chrom, seq in genome.items():
        header.contigs.add(chrom, length=len(seq))
    header.add_meta("INFO", items=[
        ("ID", "SYN_CONTEXT"),
        ("Number", "1"),
        ("Type", "String"),
        ("Description", "Synthetic mutation context")
    ])

    writer = pysam.VariantFile(out_vcf, "wz", header=header)
    for _, r in tqdm(df.iterrows(), total=len(df), desc="Writing VCF"):
        rec = writer.new_record(
            contig=str(r["chrom"]),
            start=int(r["pos"]) - 1,
            stop=int(r["pos"]),
            alleles=(str(r["ref"]), str(r["alt"]))
        )
        rec.info["SYN_CONTEXT"] = str(r["context"])
        writer.write(rec)
    writer.close()
    pysam.tabix_index(out_vcf, preset="vcf", force=True)

    print(f"Wrote: {out_vcf}")
    print(f"Wrote: {out_tsv}")
    print(df["context"].value_counts())

if __name__ == "__main__":
    main()
