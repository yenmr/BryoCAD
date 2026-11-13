#!/usr/bin/env python3

import argparse
import os
from collections import defaultdict
from contextlib import redirect_stdout
from pathlib import Path

import pandas as pd
from tqdm import tqdm
from evo2 import Evo2


TASK_RULES = {
    "donor": {"motifs": {"GT"}, "length": 2},
    "acceptor": {"motifs": {"AG"}, "length": 2},
    "tis": {"motifs": {"ATG"}, "length": 3},
    "tts": {"motifs": {"TAA", "TAG", "TGA"}, "length": 3},
}
BASES = set("ACGT")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def prepare_input_row(row, task):
    """Read the annotated motif directly; canonical status is not a filter."""
    sequence = str(row["sequence"]).strip().upper()
    motif_start = int(row["motif_start"])
    motif_end = int(row["motif_end"])
    motif_len = TASK_RULES[task]["length"]
    if motif_end - motif_start != motif_len:
        raise ValueError(
            f"motif range length is {motif_end - motif_start}, expected {motif_len}"
        )
    if motif_start < 0 or motif_end > len(sequence):
        raise ValueError("motif range is outside sequence")
    sequence_motif = sequence[motif_start:motif_end]
    annotated_motif = str(row["motif"]).strip().upper()
    if sequence_motif != annotated_motif:
        raise ValueError(
            f"motif column {annotated_motif!r} != sequence range {sequence_motif!r}"
        )
    return sequence, motif_start, motif_end, annotated_motif


def greedy_generate(model, prompts, n_tokens):
    """Generate deterministic top-1 continuations with Evo2's official API."""
    # Vortex prints every full prompt/output/score by default. Suppress that
    # output because 4-kb genomic prompts can otherwise fill the Slurm log.
    with open(os.devnull, "w") as devnull, redirect_stdout(devnull):
        output = model.generate(
            prompt_seqs=prompts,
            n_tokens=n_tokens,
            temperature=1.0,
            top_k=1,
        )
    if not hasattr(output, "sequences"):
        raise RuntimeError("Evo2 generate() result has no sequences attribute")
    if len(output.sequences) != len(prompts):
        raise RuntimeError(
            f"Generated {len(output.sequences)} sequences for {len(prompts)} prompts"
        )

    generated = []
    for prompt, sequence in zip(prompts, output.sequences):
        sequence = str(sequence).upper()
        # Vortex versions may return either continuation alone or
        # prompt + continuation. The final n_tokens are generated bases.
        continuation = sequence[-n_tokens:]
        if len(continuation) != n_tokens or any(
            base not in BASES for base in continuation
        ):
            raise RuntimeError(
                "Invalid Evo2 continuation: "
                f"prompt_length={len(prompt)}, output={sequence!r}"
            )
        generated.append(continuation)
    return generated


def score_context_batch(model, contexts, reference_motifs):
    """Compute greedy motif recovery and per-position agreement."""
    motif_len = len(reference_motifs[0])
    context_len = len(contexts[0])
    if any(len(context) != context_len for context in contexts):
        raise ValueError("Contexts in a batch have different lengths")
    if any(len(motif) != motif_len for motif in reference_motifs):
        raise ValueError("Reference motifs in a batch have different lengths")

    # Generate once. Token accuracy is the fraction of generated motif
    # positions matching the reference; motif accuracy is evaluated later.
    predicted_motifs = greedy_generate(model, contexts, motif_len)
    token_accuracy = [
        sum(pred == ref for pred, ref in zip(predicted, reference)) / motif_len
        for predicted, reference in zip(predicted_motifs, reference_motifs)
    ]
    return token_accuracy, predicted_motifs


def score_batch(model, rows, task):
    motif_len = TASK_RULES[task]["length"]
    motif_start = rows[0]["motif_start"]
    if any(row["motif_start"] != motif_start for row in rows):
        raise ValueError("A batch contains different motif start positions")

    fwd_contexts = [row["sequence"][:motif_start] for row in rows]
    fwd_references = [row["reference_motif"] for row in rows]
    fwd_token_accuracy, fwd_predictions = score_context_batch(
        model, fwd_contexts, fwd_references
    )

    # Original downstream becomes the upstream prompt in RC orientation.
    rc_contexts = [
        reverse_complement(row["sequence"][motif_start + motif_len:])
        for row in rows
    ]
    rc_references = [
        reverse_complement(row["reference_motif"])
        for row in rows
    ]
    rc_token_accuracy, rc_generated = score_context_batch(
        model, rc_contexts, rc_references
    )
    rc_predictions = [reverse_complement(motif) for motif in rc_generated]

    valid_motifs = TASK_RULES[task]["motifs"]
    results = []
    for i, row in enumerate(rows):
        result = dict(row["metadata"])
        result.update(
            {
                "_input_order": row["_input_order"],
                "reference_motif": row["reference_motif"],
                "motif_start_0based": row["motif_start"],
                "motif_end_0based_exclusive": row["motif_start"] + motif_len,
                "reference_is_canonical": row["reference_motif"] in valid_motifs,
                "fwd_predicted_motif": fwd_predictions[i],
                "fwd_token_accuracy": fwd_token_accuracy[i],
                "fwd_exact_recovery": fwd_predictions[i] == row["reference_motif"],
                "fwd_predicted_is_canonical": fwd_predictions[i] in valid_motifs,
                "rc_reference_motif": rc_references[i],
                "rc_generated_motif": rc_generated[i],
                "rc_predicted_motif": rc_predictions[i],
                "rc_token_accuracy": rc_token_accuracy[i],
                "rc_exact_recovery": rc_predictions[i] == row["reference_motif"],
                "rc_predicted_is_canonical": rc_predictions[i] in valid_motifs,
            }
        )
        results.append(result)
    return results


def print_summary(result_df):
    print("\nSummary")
    print(f"N = {len(result_df)}")
    print(
        "Evo2-fwd mean token accuracy = "
        f"{result_df['fwd_token_accuracy'].mean():.6f}"
    )
    print(
        "Evo2-fwd exact recovery = "
        f"{result_df['fwd_exact_recovery'].mean():.6f}"
    )
    print(
        "Evo2-fwd predicted canonical fraction = "
        f"{result_df['fwd_predicted_is_canonical'].mean():.6f}"
    )
    print(
        "Evo2-rc mean token accuracy = "
        f"{result_df['rc_token_accuracy'].mean():.6f}"
    )
    print(
        "Evo2-rc exact recovery = "
        f"{result_df['rc_exact_recovery'].mean():.6f}"
    )
    print(
        "Evo2-rc predicted canonical fraction = "
        f"{result_df['rc_predicted_is_canonical'].mean():.6f}"
    )
    print("\nRecovery by reference class")
    for reference_is_canonical, group in result_df.groupby(
        "reference_is_canonical", sort=False
    ):
        label = "canonical" if reference_is_canonical else "noncanonical"
        print(
            f"{label}: N={len(group)}, "
            f"fwd_exact={group['fwd_exact_recovery'].mean():.6f}, "
            f"rc_exact={group['rc_exact_recovery'].mean():.6f}, "
            f"fwd_predicted_canonical="
            f"{group['fwd_predicted_is_canonical'].mean():.6f}, "
            f"rc_predicted_canonical="
            f"{group['rc_predicted_is_canonical'].mean():.6f}"
        )


def main():
    parser = argparse.ArgumentParser(
        description="Evo2 official greedy forward/RC motif recovery"
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument(
        "--task", required=True, type=str.lower, choices=sorted(TASK_RULES)
    )
    parser.add_argument("--model", default="evo2_1b_base")
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--batch-size", type=int, default=15)
    parser.add_argument(
        "--max-samples",
        type=int,
        default=None,
        help="Optional first-N-row limit for validation runs",
    )
    args = parser.parse_args()

    if args.batch_size < 1:
        raise ValueError("--batch-size must be at least 1")

    input_path = Path(args.input)
    output_path = Path(args.output)
    df = pd.read_csv(input_path, sep="\t")
    missing = {"sequence", "geneid", "motif_start", "motif_end", "motif"} - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    if args.max_samples is not None:
        if args.max_samples < 1:
            raise ValueError("--max-samples must be at least 1")
        df = df.iloc[:args.max_samples].reset_index(drop=True)

    prepared = []
    failures = []
    metadata_columns = [column for column in df.columns if column != "sequence"]
    for input_order, (index, row) in enumerate(df.iterrows()):
        try:
            sequence, motif_start, motif_end, reference_motif = prepare_input_row(
                row, args.task
            )
            prepared.append(
                {
                    "_input_order": input_order,
                    "sequence": sequence,
                    "motif_start": motif_start,
                    "motif_end": motif_end,
                    "reference_motif": reference_motif,
                    "metadata": {column: row[column] for column in metadata_columns},
                }
            )
        except Exception as exc:
            failures.append(
                {"row": index, "geneid": row["geneid"], "error": str(exc)}
            )

    print(f"Input rows: {len(df)}")
    print(f"Task: {args.task}")
    print(f"Model: {args.model}")
    print(f"Valid rows: {len(prepared)}")
    print(f"Input-validation failures: {len(failures)}")
    if failures:
        print("First input-validation failures:")
        for failure in failures[:5]:
            print(failure)
    if not prepared:
        raise RuntimeError("No valid sequences to score")

    grouped = defaultdict(list)
    for row in prepared:
        grouped[row["motif_start"]].append(row)

    print("Loading Evo2 model")
    model = Evo2(args.model)

    all_results = []
    total_batches = sum(
        (len(rows) + args.batch_size - 1) // args.batch_size
        for rows in grouped.values()
    )
    progress = tqdm(total=total_batches, desc="Scoring batches")
    for motif_start in sorted(grouped):
        rows = grouped[motif_start]
        for start in range(0, len(rows), args.batch_size):
            all_results.extend(
                score_batch(model, rows[start:start + args.batch_size], args.task)
            )
            progress.update(1)
    progress.close()

    result_df = (
        pd.DataFrame(all_results)
        .sort_values("_input_order")
        .drop(columns="_input_order")
        .reset_index(drop=True)
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(output_path, sep="\t", index=False)
    print(f"Output written to: {output_path}")
    print_summary(result_df)


if __name__ == "__main__":
    main()
