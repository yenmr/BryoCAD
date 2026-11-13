#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import random
from pathlib import Path

import pandas as pd
from tqdm import tqdm


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def clean_gene_id(x):
    """
    Clean gene ID for matching GFF gene IDs and OG table gene IDs.

    Example:
    AagrBONN_evm.TU.Sc2ySwM_117.1025|A_agrestis_v10
    -> AagrBONN_evm.TU.Sc2ySwM_117.1025
    """
    return str(x).strip().split()[0].split("|")[0]


def read_fasta(path):
    genome = {}
    chrom = None
    seqs = []

    with open_text(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if chrom is not None:
                    genome[chrom] = "".join(seqs).upper()
                chrom = line[1:].split()[0]
                seqs = []
            else:
                seqs.append(line)

        if chrom is not None:
            genome[chrom] = "".join(seqs).upper()

    return genome


def parse_attrs(attr_string):
    out = {}
    for item in attr_string.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        elif " " in item:
            k, v = item.split(" ", 1)
            out[k] = v.strip('"')
    return out


def read_gff_cds(path):
    """
    Parse GFF and return CDS rows with both transcript_id and gene_id.

    Logic:
    1. Read gene features:
       gene ID = ID
    2. Read mRNA/transcript features:
       transcript ID = ID
       gene ID = Parent
    3. Read CDS features:
       transcript ID = Parent
       gene ID = transcript_to_gene[Parent]
    """

    gene_ids = set()
    transcript_to_gene = {}
    cds_raw = []

    with open_text(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attrs = parts
            attr = parse_attrs(attrs)

            if feature == "gene":
                gene_id = attr.get("ID")
                if gene_id:
                    gene_ids.add(gene_id)

            elif feature in {"mRNA", "transcript"}:
                transcript_id = attr.get("ID")
                parent_gene = attr.get("Parent")

                if transcript_id and parent_gene:
                    # Parent can sometimes be comma-separated
                    parent_gene = parent_gene.split(",")[0]
                    transcript_to_gene[transcript_id] = parent_gene

            elif feature == "CDS":
                parent = attr.get("Parent")
                cds_id = attr.get("ID", "NA")

                if parent is None:
                    continue

                # CDS Parent may contain multiple transcripts
                for transcript_id in parent.split(","):
                    cds_raw.append({
                        "chrom": chrom,
                        "start0": int(start) - 1,
                        "end0": int(end),
                        "strand": strand,
                        "phase": phase,
                        "transcript_id": transcript_id,
                        "cds_id": cds_id,
                    })

    rows = []

    for r in cds_raw:
        transcript_id = r["transcript_id"]

        gene_id = transcript_to_gene.get(transcript_id)

        # fallback: sometimes CDS Parent directly points to gene
        if gene_id is None:
            if transcript_id in gene_ids:
                gene_id = transcript_id
            else:
                gene_id = transcript_id

        rows.append({
            "chrom": r["chrom"],
            "start0": r["start0"],
            "end0": r["end0"],
            "strand": r["strand"],
            "phase": r["phase"],
            "transcript_id": transcript_id,
            "gene_id": gene_id,
            "cds_id": r["cds_id"],
        })

    df = pd.DataFrame(rows)

    if len(df) == 0:
        raise ValueError(f"No CDS entries parsed from {path}")

    return df


def genomic_interval_to_oriented(start0, end0, chrom_len, strand):
    if strand == "+":
        return start0, end0
    else:
        return chrom_len - end0, chrom_len - start0


def oriented_interval_to_genomic(start_ori, end_ori, chrom_len, strand):
    if strand == "+":
        return start_ori, end_ori
    else:
        return chrom_len - end_ori, chrom_len - start_ori


def oriented_pos_to_genomic(pos_ori, chrom_len, strand):
    if strand == "+":
        return pos_ori
    else:
        return chrom_len - 1 - pos_ori


def fetch_oriented(genome, chrom, start_ori, end_ori, strand):
    chrom_len = len(genome[chrom])

    if start_ori < 0 or end_ori > chrom_len or start_ori >= end_ori:
        return None

    if strand == "+":
        return genome[chrom][start_ori:end_ori].upper()
    else:
        g_start = chrom_len - end_ori
        g_end = chrom_len - start_ori
        return revcomp(genome[chrom][g_start:g_end])


def make_ref_mut_sequences(
    genome,
    chrom,
    strand,
    anchor_ori,
    del_start_ori,
    del_end_ori,
    seq_length,
    mode,
):
    """
    mode:
      - deletion_after_anchor: deletion starts at anchor, used for CDS deletion
      - deletion_before_anchor: deletion ends at anchor, used for TIS upstream deletion
    """

    half = seq_length // 2

    if mode == "deletion_after_anchor":
        # ref: centered at deletion start
        ref_start = anchor_ori - half
        ref_end = anchor_ori + half

        # mut: left half before deletion + right half after deletion
        mut_left_start = anchor_ori - half
        mut_left_end = anchor_ori
        mut_right_start = del_end_ori
        mut_right_end = del_end_ori + half

    elif mode == "deletion_before_anchor":
        # ref: centered at TIS anchor
        ref_start = anchor_ori - half
        ref_end = anchor_ori + half

        # mut: take extra upstream sequence, delete region immediately before anchor
        mut_left_start = anchor_ori - half - (del_end_ori - del_start_ori)
        mut_left_end = del_start_ori
        mut_right_start = anchor_ori
        mut_right_end = anchor_ori + half

    else:
        raise ValueError(f"Unknown mode: {mode}")

    ref_seq = fetch_oriented(genome, chrom, ref_start, ref_end, strand)
    mut_left = fetch_oriented(genome, chrom, mut_left_start, mut_left_end, strand)
    mut_right = fetch_oriented(genome, chrom, mut_right_start, mut_right_end, strand)

    if ref_seq is None or mut_left is None or mut_right is None:
        return None, None

    mut_seq = mut_left + mut_right

    if len(ref_seq) != seq_length or len(mut_seq) != seq_length:
        return None, None

    return ref_seq, mut_seq


def n_ratio(seq):
    return seq.count("N") / len(seq)


def build_transcripts(cds_df, genome):
    transcripts = []

    for tid, sub in cds_df.groupby("transcript_id"):
        chrom = sub["chrom"].iloc[0]
        strand = sub["strand"].iloc[0]
        gene_id = sub["gene_id"].iloc[0]

        if chrom not in genome:
            continue

        chrom_len = len(genome[chrom])

        intervals = []
        for _, r in sub.iterrows():
            s_ori, e_ori = genomic_interval_to_oriented(
                r["start0"],
                r["end0"],
                chrom_len,
                strand,
            )
            intervals.append((s_ori, e_ori))

        intervals = sorted(intervals)

        if len(intervals) == 0:
            continue

        tis_ori = intervals[0][0]

        transcripts.append({
            "transcript_id": tid,
            "gene_id": gene_id,
            "chrom": chrom,
            "strand": strand,
            "cds_intervals_ori": intervals,
            "tis_ori": tis_ori,
        })

    return transcripts


def make_tis_upstream_candidate(
    genome,
    species,
    transcript,
    seq_length,
    deletion_length,
    max_n_ratio,
):
    chrom = transcript["chrom"]
    strand = transcript["strand"]
    chrom_len = len(genome[chrom])

    anchor_ori = transcript["tis_ori"]
    del_start_ori = anchor_ori - deletion_length
    del_end_ori = anchor_ori

    ref_seq, mut_seq = make_ref_mut_sequences(
        genome=genome,
        chrom=chrom,
        strand=strand,
        anchor_ori=anchor_ori,
        del_start_ori=del_start_ori,
        del_end_ori=del_end_ori,
        seq_length=seq_length,
        mode="deletion_before_anchor",
    )

    if ref_seq is None:
        return None

    if n_ratio(ref_seq) > max_n_ratio or n_ratio(mut_seq) > max_n_ratio:
        return None

    g_start0, g_end0 = oriented_interval_to_genomic(
        del_start_ori,
        del_end_ori,
        chrom_len,
        strand,
    )

    anchor_g0 = oriented_pos_to_genomic(anchor_ori, chrom_len, strand)

    return {
        "sv_id": None,
        "species": species,
        "sv_type": "TIS_upstream_deletion",
        "chromosome": chrom,
        "strand": strand,
        "gene_id": transcript["gene_id"],
        "transcript_id": transcript["transcript_id"],
        "anchor_type": "TIS_ATG",
        "anchor_genomic_1based": anchor_g0 + 1,
        "anchor_oriented_0based": anchor_ori,
        "deletion_length": deletion_length,
        "deletion_genomic_start_1based": g_start0 + 1,
        "deletion_genomic_end_1based": g_end0,
        "deletion_oriented_start_0based": del_start_ori,
        "deletion_oriented_end_0based": del_end_ori,
        "is_cds": False,
        "is_in_frame": None,
        "ref_sequence": ref_seq,
        "mut_sequence": mut_seq,
        "ref_n_ratio": n_ratio(ref_seq),
        "mut_n_ratio": n_ratio(mut_seq),
    }


def choose_cds_deletion_from_transcript(
    genome,
    species,
    transcript,
    seq_length,
    deletion_length,
    cds_flank,
    max_n_ratio,
    rng,
    sv_type,
    require_codon_boundary=True,
):
    chrom = transcript["chrom"]
    strand = transcript["strand"]
    chrom_len = len(genome[chrom])

    valid_sites = []

    for exon_idx, (exon_s, exon_e) in enumerate(transcript["cds_intervals_ori"]):
        exon_len = exon_e - exon_s

        if exon_len < deletion_length + 2 * cds_flank:
            continue

        min_offset = cds_flank
        max_offset = exon_len - deletion_length - cds_flank

        for offset in range(min_offset, max_offset + 1):
            if require_codon_boundary and offset % 3 != 0:
                continue

            del_start_ori = exon_s + offset
            del_end_ori = del_start_ori + deletion_length
            valid_sites.append((exon_idx, del_start_ori, del_end_ori, offset))

    if not valid_sites:
        return None

    exon_idx, del_start_ori, del_end_ori, offset = rng.choice(valid_sites)

    anchor_ori = del_start_ori

    ref_seq, mut_seq = make_ref_mut_sequences(
        genome=genome,
        chrom=chrom,
        strand=strand,
        anchor_ori=anchor_ori,
        del_start_ori=del_start_ori,
        del_end_ori=del_end_ori,
        seq_length=seq_length,
        mode="deletion_after_anchor",
    )

    if ref_seq is None:
        return None

    if n_ratio(ref_seq) > max_n_ratio or n_ratio(mut_seq) > max_n_ratio:
        return None

    g_start0, g_end0 = oriented_interval_to_genomic(
        del_start_ori,
        del_end_ori,
        chrom_len,
        strand,
    )

    anchor_g0 = oriented_pos_to_genomic(anchor_ori, chrom_len, strand)

    return {
        "sv_id": None,
        "species": species,
        "sv_type": sv_type,
        "chromosome": chrom,
        "strand": strand,
        "gene_id": transcript["gene_id"],
        "transcript_id": transcript["transcript_id"],
        "anchor_type": "CDS_deletion_start",
        "anchor_genomic_1based": anchor_g0 + 1,
        "anchor_oriented_0based": anchor_ori,
        "deletion_length": deletion_length,
        "deletion_genomic_start_1based": g_start0 + 1,
        "deletion_genomic_end_1based": g_end0,
        "deletion_oriented_start_0based": del_start_ori,
        "deletion_oriented_end_0based": del_end_ori,
        "cds_exon_index": exon_idx,
        "cds_offset_in_exon": offset,
        "is_cds": True,
        "is_in_frame": deletion_length % 3 == 0,
        "ref_sequence": ref_seq,
        "mut_sequence": mut_seq,
        "ref_n_ratio": n_ratio(ref_seq),
        "mut_n_ratio": n_ratio(mut_seq),
    }


def choose_deep_intron_deletion_from_transcript(
    genome,
    species,
    transcript,
    seq_length,
    deletion_length,
    intron_flank,
    max_n_ratio,
    rng,
):
    chrom = transcript["chrom"]
    strand = transcript["strand"]
    chrom_len = len(genome[chrom])

    intervals = sorted(transcript["cds_intervals_ori"])

    if len(intervals) < 2:
        return None

    valid_sites = []

    for intron_idx in range(len(intervals) - 1):
        left_exon_s, left_exon_e = intervals[intron_idx]
        right_exon_s, right_exon_e = intervals[intron_idx + 1]

        intron_s = left_exon_e
        intron_e = right_exon_s
        intron_len = intron_e - intron_s

        if intron_len < deletion_length + 2 * intron_flank:
            continue

        min_start = intron_s + intron_flank
        max_start = intron_e - intron_flank - deletion_length

        if max_start <= min_start:
            continue

        for del_start_ori in range(min_start, max_start + 1):
            del_end_ori = del_start_ori + deletion_length
            valid_sites.append((intron_idx, del_start_ori, del_end_ori))

    if not valid_sites:
        return None

    intron_idx, del_start_ori, del_end_ori = rng.choice(valid_sites)

    anchor_ori = del_start_ori

    ref_seq, mut_seq = make_ref_mut_sequences(
        genome=genome,
        chrom=chrom,
        strand=strand,
        anchor_ori=anchor_ori,
        del_start_ori=del_start_ori,
        del_end_ori=del_end_ori,
        seq_length=seq_length,
        mode="deletion_after_anchor",
    )

    if ref_seq is None:
        return None

    if n_ratio(ref_seq) > max_n_ratio or n_ratio(mut_seq) > max_n_ratio:
        return None

    g_start0, g_end0 = oriented_interval_to_genomic(
        del_start_ori,
        del_end_ori,
        chrom_len,
        strand,
    )

    anchor_g0 = oriented_pos_to_genomic(anchor_ori, chrom_len, strand)

    return {
        "sv_id": None,
        "species": species,
        "sv_type": "Deep_intron_deletion",
        "chromosome": chrom,
        "strand": strand,
        "gene_id": transcript["gene_id"],
        "transcript_id": transcript["transcript_id"],
        "anchor_type": "deep_intron_deletion_start",
        "anchor_genomic_1based": anchor_g0 + 1,
        "anchor_oriented_0based": anchor_ori,
        "deletion_length": deletion_length,
        "deletion_genomic_start_1based": g_start0 + 1,
        "deletion_genomic_end_1based": g_end0,
        "deletion_oriented_start_0based": del_start_ori,
        "deletion_oriented_end_0based": del_end_ori,
        "intron_index": intron_idx,
        "is_cds": False,
        "is_in_frame": None,
        "ref_sequence": ref_seq,
        "mut_sequence": mut_seq,
        "ref_n_ratio": n_ratio(ref_seq),
        "mut_n_ratio": n_ratio(mut_seq),
    }


def sample_rows(rows, max_n, rng):
    if max_n is None or max_n <= 0:
        return rows
    if len(rows) <= max_n:
        return rows
    return rng.sample(rows, max_n)


def process_one_genome(
    species,
    fasta_path,
    gff_path,
    args,
    rng,
    allowed_genes=None,
):
    print(f"\nProcessing {species}")
    print(f"  FASTA: {fasta_path}")
    print(f"  GFF  : {gff_path}")

    genome = read_fasta(fasta_path)
    cds_df = read_gff_cds(gff_path)

    transcripts = build_transcripts(cds_df, genome)

    print(f"  transcripts with CDS: {len(transcripts):,}")

    tis_rows = []
    cds_inframe_rows = []
    cds_frameshift_rows = []
    deep_intron_rows = []

    for tx in tqdm(transcripts):
        gene_clean = clean_gene_id(tx["gene_id"])
        if allowed_genes is not None:
            if (species, gene_clean) not in allowed_genes:
                continue

        row = make_tis_upstream_candidate(
            genome=genome,
            species=species,
            transcript=tx,
            seq_length=args.seq_length,
            deletion_length=args.tis_upstream_del_len,
            max_n_ratio=args.max_n_ratio,
        )
        if row is not None:
            tis_rows.append(row)

        row = choose_cds_deletion_from_transcript(
            genome=genome,
            species=species,
            transcript=tx,
            seq_length=args.seq_length,
            deletion_length=args.cds_inframe_del_len,
            cds_flank=args.cds_flank,
            max_n_ratio=args.max_n_ratio,
            rng=rng,
            sv_type="CDS_inframe_deletion",
            require_codon_boundary=True,
        )
        if row is not None:
            cds_inframe_rows.append(row)

        row = choose_cds_deletion_from_transcript(
            genome=genome,
            species=species,
            transcript=tx,
            seq_length=args.seq_length,
            deletion_length=args.cds_frameshift_del_len,
            cds_flank=args.cds_flank,
            max_n_ratio=args.max_n_ratio,
            rng=rng,
            sv_type="CDS_frameshift_deletion",
            require_codon_boundary=True,
        )
        if row is not None:
            cds_frameshift_rows.append(row)

        row = choose_deep_intron_deletion_from_transcript(
            genome=genome,
            species=species,
            transcript=tx,
            seq_length=args.seq_length,
            deletion_length=args.deep_intron_del_len,
            intron_flank=args.intron_flank,
            max_n_ratio=args.max_n_ratio,
            rng=rng,
        )
        if row is not None:
            deep_intron_rows.append(row)

    tis_rows = sample_rows(tis_rows, args.max_per_type_per_genome, rng)
    cds_inframe_rows = sample_rows(cds_inframe_rows, args.max_per_type_per_genome, rng)
    cds_frameshift_rows = sample_rows(cds_frameshift_rows, args.max_per_type_per_genome, rng)
    deep_intron_rows = sample_rows(deep_intron_rows, args.max_per_type_per_genome, rng)

    rows = tis_rows + cds_inframe_rows + cds_frameshift_rows + deep_intron_rows

    print(f"  TIS upstream deletion: {len(tis_rows):,}")
    print(f"  CDS in-frame deletion: {len(cds_inframe_rows):,}")
    print(f"  CDS frameshift deletion: {len(cds_frameshift_rows):,}")
    print(f"  Deep intron deletion: {len(deep_intron_rows):,}")

    return rows


def read_manifest(path):
    df = pd.read_csv(path, sep="\t")

    required = {"species", "fasta", "gff"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"Manifest missing columns: {missing}")

    return df


def load_allowed_core_genes(core_gene_table, core_only=False, shared_core_only=False):
    if core_gene_table is None:
        return None

    df = pd.read_csv(core_gene_table, sep="\t")

    required = {"geneid", "species", "OG_category"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"core_gene_table missing columns: {missing}")

    df["geneid_clean"] = df["geneid"].map(clean_gene_id)
    df["OG_category"] = df["OG_category"].astype(str)

    if shared_core_only:
        if "shared_core" not in df.columns:
            raise ValueError("--shared_core_only requires shared_core column.")

        df["shared_core_bool"] = (
            df["shared_core"]
            .astype(str)
            .str.lower()
            .isin(["true", "1", "yes"])
        )

        df = df[
            (df["OG_category"] == "core") &
            (df["shared_core_bool"])
        ].copy()

    elif core_only:
        df = df[df["OG_category"] == "core"].copy()

    allowed = set(zip(df["species"].astype(str), df["geneid_clean"].astype(str)))

    print(f"[INFO] allowed genes from core table: {len(allowed):,}")

    return allowed


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--manifest", required=True,
                    help="TSV with columns: species, fasta, gff")
    ap.add_argument("--out", default="synthetic_deletion_input.tsv")

    ap.add_argument("--seq_length", type=int, default=8192)
    ap.add_argument("--max_n_ratio", type=float, default=0.05)

    ap.add_argument("--tis_upstream_del_len", type=int, default=100)
    ap.add_argument("--cds_inframe_del_len", type=int, default=9)
    ap.add_argument("--cds_frameshift_del_len", type=int, default=10)
    ap.add_argument("--deep_intron_del_len", type=int, default=10)

    ap.add_argument("--intron_flank", type=int, default=100,
                    help="Require deep intron deletion to be this many bp away from CDS/exon boundary.")
    
    ap.add_argument("--cds_flank", type=int, default=10,
                    help="Require deletion to be at least this many bp away from CDS exon boundary.")

    ap.add_argument("--core_gene_table", default=None,
                    help="TSV with geneid, species, OG_category, shared_core columns.")

    ap.add_argument("--core_only", action="store_true",
                    help="Only use genes with OG_category == core.")

    ap.add_argument("--max_per_type_per_genome", type=int, default=1000,
                    help="Maximum rows per SV type per genome. 0 = keep all.")

    ap.add_argument("--seed", type=int, default=42)

    args = ap.parse_args()

    if args.seq_length % 2 != 0:
        raise ValueError("--seq_length must be even.")

    rng = random.Random(args.seed)

    manifest = read_manifest(args.manifest)

    allowed_genes = load_allowed_core_genes(
        core_gene_table=args.core_gene_table,
        core_only=args.core_only,
    )

    all_rows = []

    for _, r in manifest.iterrows():
        rows = process_one_genome(
            species=r["species"],
            fasta_path=r["fasta"],
            gff_path=r["gff"],
            args=args,
            rng=rng,
            allowed_genes=allowed_genes,
        )
        all_rows.extend(rows)



    if not all_rows:
        raise ValueError("No rows generated.")

    for i, row in enumerate(all_rows, start=1):
        row["sv_id"] = f"synSV_{i:07d}"

    df = pd.DataFrame(all_rows)

    out = Path(args.out)
    df.to_csv(out, sep="\t", index=False)

    print(f"\nWrote: {out}")
    print(f"Rows: {len(df):,}")

    print("\nCounts:")
    print(df.groupby(["species", "sv_type"]).size())


if __name__ == "__main__":
    main()
