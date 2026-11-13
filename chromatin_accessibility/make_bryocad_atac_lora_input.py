#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build BryoCAD ATAC-seq LoRA binary-classification input.

Positive:
    A fixed-length sequence centered on each ATAC-seq peak summit.

Negative:
    Fixed-length whole-genome windows that do not overlap any ATAC-seq peak
    by even 1 bp.

Filtering:
    Remove any positive or negative sequence containing more than --max_n_bp Ns.

Splitting:
    Chromosome-level train/validation split to prevent nearby genomic regions
    from leaking between training and validation sets.

Negative sampling:
    --neg_per_pos > 0:
        Randomly sample this many eligible negative windows per retained positive
        within each split.
    --neg_per_pos 0:
        Keep all eligible whole-genome negative windows.
"""

import argparse
import gzip
import random
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def read_fasta(path: str) -> Dict[str, str]:
    genome: Dict[str, str] = {}
    chrom: Optional[str] = None
    seqs: List[str] = []

    with open_text(path) as handle:
        for line in handle:
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

    if not genome:
        raise ValueError(f"No FASTA records found in {path}")

    return genome


def gc_content(seq: str) -> float:
    seq = seq.upper()
    a = seq.count("A")
    c = seq.count("C")
    g = seq.count("G")
    t = seq.count("T")
    valid = a + c + g + t

    if valid == 0:
        return np.nan

    return (g + c) / valid


def count_n(seq: str) -> int:
    return seq.upper().count("N")


def extract_interval(
    genome: Dict[str, str],
    chrom: str,
    start: int,
    end: int,
) -> Optional[str]:
    if chrom not in genome:
        return None
    if start < 0 or end > len(genome[chrom]) or start >= end:
        return None

    seq = genome[chrom][start:end].upper()
    if len(seq) != end - start:
        return None

    return seq


def extract_centered_window(
    genome: Dict[str, str],
    chrom: str,
    center: int,
    window: int,
) -> Tuple[Optional[str], Optional[int], Optional[int]]:
    half = window // 2
    left = center - half
    right = left + window
    seq = extract_interval(genome, chrom, left, right)

    if seq is None or len(seq) != window:
        return None, None, None

    return seq, left, right


def read_narrowpeak(path: str) -> pd.DataFrame:
    cols = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "signalValue",
        "pValue",
        "qValue",
        "summit_offset",
    ]

    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=cols,
        usecols=range(10),
        comment="#",
    )

    for col in ["start", "end", "summit_offset"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["chrom", "start", "end", "summit_offset"]).copy()
    df[["start", "end", "summit_offset"]] = df[
        ["start", "end", "summit_offset"]
    ].astype(int)

    df = df[df["end"] > df["start"]].copy()
    df["summit"] = df["start"] + df["summit_offset"]
    df["peak_id"] = [f"peak_{i}" for i in range(len(df))]

    return df


def merge_intervals(
    starts: Sequence[int],
    ends: Sequence[int],
) -> Tuple[np.ndarray, np.ndarray]:
    if len(starts) == 0:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)

    starts_array = np.asarray(starts, dtype=np.int64)
    ends_array = np.asarray(ends, dtype=np.int64)
    order = np.argsort(starts_array)

    starts_array = starts_array[order]
    ends_array = ends_array[order]

    merged_starts: List[int] = []
    merged_ends: List[int] = []

    cur_start = int(starts_array[0])
    cur_end = int(ends_array[0])

    for start, end in zip(starts_array[1:], ends_array[1:]):
        start = int(start)
        end = int(end)

        if start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            merged_starts.append(cur_start)
            merged_ends.append(cur_end)
            cur_start = start
            cur_end = end

    merged_starts.append(cur_start)
    merged_ends.append(cur_end)

    return (
        np.asarray(merged_starts, dtype=np.int64),
        np.asarray(merged_ends, dtype=np.int64),
    )


def make_peak_interval_index(
    peaks: pd.DataFrame,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    interval_starts: Dict[str, np.ndarray] = {}
    interval_ends: Dict[str, np.ndarray] = {}

    for chrom, sub in peaks.groupby("chrom", sort=False):
        starts, ends = merge_intervals(
            sub["start"].to_numpy(),
            sub["end"].to_numpy(),
        )
        interval_starts[chrom] = starts
        interval_ends[chrom] = ends

    return interval_starts, interval_ends


def interval_overlaps_peak(
    interval_starts: Dict[str, np.ndarray],
    interval_ends: Dict[str, np.ndarray],
    chrom: str,
    start: int,
    end: int,
) -> bool:
    """Return True when [start, end) overlaps any merged peak interval."""
    if chrom not in interval_starts:
        return False

    starts = interval_starts[chrom]
    ends = interval_ends[chrom]

    # Candidate peak immediately before `end`.
    idx = np.searchsorted(starts, end, side="left") - 1

    # BED intervals overlap when peak_start < end and peak_end > start.
    return idx >= 0 and int(ends[idx]) > start


def split_chromosomes(
    peaks: pd.DataFrame,
    val_fraction: float,
    seed: int,
    n_trials: int = 5000,
) -> Tuple[Set[str], Set[str]]:
    rng = random.Random(seed)
    counts = peaks["chrom"].value_counts().to_dict()
    chroms = list(counts)

    if len(chroms) < 2:
        raise ValueError(
            "Chromosome-level train/validation splitting requires at least two "
            "chromosomes containing peaks."
        )

    total = sum(counts.values())
    target = total * val_fraction

    best_set: Optional[Set[str]] = None
    best_diff = float("inf")

    for _ in range(n_trials):
        shuffled = chroms.copy()
        rng.shuffle(shuffled)

        val: Set[str] = set()
        current = 0

        for chrom in shuffled:
            if current < target:
                val.add(chrom)
                current += counts[chrom]

        if len(val) == len(chroms):
            val.remove(shuffled[-1])
            current -= counts[shuffled[-1]]

        diff = abs(current - target)

        if diff < best_diff:
            best_diff = diff
            best_set = set(val)

    if not best_set:
        raise RuntimeError("Unable to create a non-empty validation chromosome set.")

    train = set(chroms) - best_set
    return train, best_set


def build_positive_rows(
    genome: Dict[str, str],
    peaks: pd.DataFrame,
    val_chroms: Set[str],
    seq_length: int,
    max_n_bp: int,
) -> Tuple[List[dict], dict]:
    rows: List[dict] = []
    stats = {
        "input_peaks": len(peaks),
        "kept_positive": 0,
        "skipped_positive_boundary": 0,
        "skipped_positive_n": 0,
    }

    for peak in tqdm(
        peaks.itertuples(index=False),
        total=len(peaks),
        desc="Positive windows",
    ):
        chrom = peak.chrom
        summit = int(peak.summit)

        seq, start, end = extract_centered_window(
            genome=genome,
            chrom=chrom,
            center=summit,
            window=seq_length,
        )

        if seq is None:
            stats["skipped_positive_boundary"] += 1
            continue

        n_bp = count_n(seq)
        if n_bp > max_n_bp:
            stats["skipped_positive_n"] += 1
            continue

        rows.append(
            {
                "sequence": seq,
                "label": 1,
                "chromosome": chrom,
                "center": summit,
                "start": start,
                "end": end,
                "source": "peak_summit",
                "peak_id": peak.peak_id,
                "summit": summit,
                "gc_content": gc_content(seq),
                "n_bp": n_bp,
                "split": "validation" if chrom in val_chroms else "train",
            }
        )
        stats["kept_positive"] += 1

    return rows, stats


def generate_negative_candidates(
    genome: Dict[str, str],
    interval_starts: Dict[str, np.ndarray],
    interval_ends: Dict[str, np.ndarray],
    train_chroms: Set[str],
    val_chroms: Set[str],
    seq_length: int,
    window_stride: int,
    max_n_bp: int,
) -> Tuple[Dict[str, List[dict]], dict]:
    candidates: Dict[str, List[dict]] = {"train": [], "validation": []}
    stats = {
        "genome_windows_examined": 0,
        "negative_windows_overlapping_peak": 0,
        "negative_windows_excess_n": 0,
        "eligible_negative_windows": 0,
    }

    selected_chroms = train_chroms | val_chroms

    for chrom in tqdm(sorted(selected_chroms), desc="Whole-genome negatives"):
        chrom_seq = genome[chrom]
        chrom_len = len(chrom_seq)
        split = "validation" if chrom in val_chroms else "train"

        if chrom_len < seq_length:
            continue

        # Fixed genomic windows: [start, start + seq_length)
        for start in range(0, chrom_len - seq_length + 1, window_stride):
            end = start + seq_length
            stats["genome_windows_examined"] += 1

            if interval_overlaps_peak(
                interval_starts=interval_starts,
                interval_ends=interval_ends,
                chrom=chrom,
                start=start,
                end=end,
            ):
                stats["negative_windows_overlapping_peak"] += 1
                continue

            seq = chrom_seq[start:end]
            n_bp = count_n(seq)

            if n_bp > max_n_bp:
                stats["negative_windows_excess_n"] += 1
                continue

            candidates[split].append(
                {
                    "sequence": seq,
                    "label": 0,
                    "chromosome": chrom,
                    "center": start + seq_length // 2,
                    "start": start,
                    "end": end,
                    "source": "whole_genome_nonpeak",
                    "peak_id": "",
                    "summit": "",
                    "gc_content": gc_content(seq),
                    "n_bp": n_bp,
                    "split": split,
                }
            )
            stats["eligible_negative_windows"] += 1

    return candidates, stats


def select_negative_rows(
    positive_rows: List[dict],
    negative_candidates: Dict[str, List[dict]],
    neg_per_pos: int,
    seed: int,
) -> Tuple[List[dict], dict]:
    rng = random.Random(seed)
    selected: List[dict] = []
    stats = {
        "kept_negative": 0,
        "available_train_negative": len(negative_candidates["train"]),
        "available_validation_negative": len(negative_candidates["validation"]),
    }

    positive_counts = {
        "train": sum(row["split"] == "train" for row in positive_rows),
        "validation": sum(row["split"] == "validation" for row in positive_rows),
    }

    for split in ["train", "validation"]:
        pool = negative_candidates[split]

        if neg_per_pos == 0:
            chosen = pool
        else:
            requested = positive_counts[split] * neg_per_pos
            if len(pool) < requested:
                raise ValueError(
                    f"Not enough eligible {split} negatives: requested {requested}, "
                    f"available {len(pool)}. Reduce --neg_per_pos, reduce "
                    f"--window_stride, increase --max_n_bp, or use --neg_per_pos 0."
                )
            chosen = rng.sample(pool, requested)

        selected.extend(chosen)

    stats["kept_negative"] = len(selected)
    return selected, stats


def write_outputs(
    rows: List[dict],
    outdir: Path,
    summary_metrics: dict,
    args: argparse.Namespace,
) -> None:
    if not rows:
        raise ValueError("No samples generated.")

    df = pd.DataFrame(rows)
    df["label"] = df["label"].astype(int)

    # Reproducible row shuffling within each split.
    train_df = (
        df[df["split"] == "train"]
        .sample(frac=1, random_state=args.seed)
        .reset_index(drop=True)
    )
    val_df = (
        df[df["split"] == "validation"]
        .sample(frac=1, random_state=args.seed)
        .reset_index(drop=True)
    )

    if train_df.empty or val_df.empty:
        raise ValueError("Train or validation dataset is empty.")

    outdir.mkdir(parents=True, exist_ok=True)

    train_tsv = outdir / "Train.tsv"
    val_tsv = outdir / "Validation.tsv"
    summary_out = outdir / "dataset_summary.tsv"

    train_df.to_csv(train_tsv, sep="\t", index=False)
    val_df.to_csv(val_tsv, sep="\t", index=False)

    if args.write_parquet:
        train_df.to_parquet(
            outdir / "Train.parquet",
            index=False,
            compression="zstd",
        )
        val_df.to_parquet(
            outdir / "Validation.parquet",
            index=False,
            compression="zstd",
        )

    final_metrics = dict(summary_metrics)
    final_metrics.update(
        {
            "train_rows": len(train_df),
            "validation_rows": len(val_df),
            "train_positive": int((train_df["label"] == 1).sum()),
            "train_negative": int((train_df["label"] == 0).sum()),
            "validation_positive": int((val_df["label"] == 1).sum()),
            "validation_negative": int((val_df["label"] == 0).sum()),
            "seq_length": args.seq_length,
            "window_stride": args.window_stride,
            "neg_per_pos": args.neg_per_pos,
            "max_n_bp": args.max_n_bp,
            "val_fraction": args.val_fraction,
            "seed": args.seed,
        }
    )

    summary = pd.DataFrame(
        [{"metric": key, "value": value} for key, value in final_metrics.items()]
    )
    summary.to_csv(summary_out, sep="\t", index=False)

    print(f"\nWrote: {train_tsv}")
    print(f"Wrote: {val_tsv}")
    print(f"Wrote: {summary_out}")

    if args.write_parquet:
        print(f"Wrote: {outdir / 'Train.parquet'}")
        print(f"Wrote: {outdir / 'Validation.parquet'}")

    print("\nLabel counts:")
    print(df.groupby(["split", "label"]).size())

    print("\nGC summary:")
    print(df.groupby(["split", "label"])["gc_content"].describe())

    print("\nN-bp summary:")
    print(df.groupby(["split", "label"])["n_bp"].describe())


def build_dataset(args: argparse.Namespace) -> None:
    print("Reading genome...")
    genome = read_fasta(args.fasta)

    print("Reading narrowPeak...")
    peaks = read_narrowpeak(args.narrowpeak)
    peaks = peaks[peaks["chrom"].isin(genome)].copy()

    if peaks.empty:
        raise ValueError(
            "No peaks remain after matching narrowPeak chromosome names "
            "against the FASTA."
        )

    train_chroms, val_chroms = split_chromosomes(
        peaks=peaks,
        val_fraction=args.val_fraction,
        seed=args.seed,
    )

    print(f"Train chromosomes: {sorted(train_chroms)}")
    print(f"Validation chromosomes: {sorted(val_chroms)}")

    interval_starts, interval_ends = make_peak_interval_index(peaks)

    positive_rows, positive_stats = build_positive_rows(
        genome=genome,
        peaks=peaks,
        val_chroms=val_chroms,
        seq_length=args.seq_length,
        max_n_bp=args.max_n_bp,
    )

    if not positive_rows:
        raise ValueError("No positive samples remain after boundary/N filtering.")

    negative_candidates, negative_candidate_stats = generate_negative_candidates(
        genome=genome,
        interval_starts=interval_starts,
        interval_ends=interval_ends,
        train_chroms=train_chroms,
        val_chroms=val_chroms,
        seq_length=args.seq_length,
        window_stride=args.window_stride,
        max_n_bp=args.max_n_bp,
    )

    negative_rows, negative_selection_stats = select_negative_rows(
        positive_rows=positive_rows,
        negative_candidates=negative_candidates,
        neg_per_pos=args.neg_per_pos,
        seed=args.seed,
    )

    all_rows = positive_rows + negative_rows
    summary_metrics = {}
    summary_metrics.update(positive_stats)
    summary_metrics.update(negative_candidate_stats)
    summary_metrics.update(negative_selection_stats)

    write_outputs(
        rows=all_rows,
        outdir=Path(args.outdir),
        summary_metrics=summary_metrics,
        args=args,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Create ATAC LoRA inputs using peak-summit positives and "
            "whole-genome non-peak windows as negatives."
        )
    )

    parser.add_argument("--fasta", required=True, help="Genome FASTA or FASTA.gz")
    parser.add_argument(
        "--narrowpeak",
        required=True,
        help="MACS2 narrowPeak file",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory",
    )

    parser.add_argument(
        "--seq_length",
        type=int,
        default=512,
        help="Length of positive and negative sequences. Default: 512",
    )
    parser.add_argument(
        "--window_stride",
        type=int,
        default=None,
        help=(
            "Stride between whole-genome negative windows. "
            "Default: seq_length, producing non-overlapping genomic tiles."
        ),
    )
    parser.add_argument(
        "--neg_per_pos",
        type=int,
        default=0,
        help=(
            "Number of eligible negatives sampled per positive in each split. "
            "Use 0 to keep every eligible whole-genome negative window. Default: 0"
        ),
    )
    parser.add_argument(
        "--max_n_bp",
        type=int,
        default=0,
        help=(
            "Maximum number of N bases allowed in any sequence. "
            "For example, --max_n_bp 10 removes windows with >10 Ns. Default: 0"
        ),
    )

    parser.add_argument("--val_fraction", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--write_parquet", action="store_true")

    args = parser.parse_args()

    if args.seq_length <= 0 or args.seq_length % 2 != 0:
        raise ValueError("--seq_length must be a positive even integer.")

    if args.window_stride is None:
        args.window_stride = args.seq_length

    if args.window_stride <= 0:
        raise ValueError("--window_stride must be > 0.")

    if args.neg_per_pos < 0:
        raise ValueError("--neg_per_pos must be >= 0.")

    if args.max_n_bp < 0:
        raise ValueError("--max_n_bp must be >= 0.")

    if not 0 < args.val_fraction < 1:
        raise ValueError("--val_fraction must be between 0 and 1.")

    build_dataset(args)


if __name__ == "__main__":
    main()
