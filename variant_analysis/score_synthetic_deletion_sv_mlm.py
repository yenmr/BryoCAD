#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm
from transformers import AutoModelForMaskedLM, AutoTokenizer


NUCLEOTIDES = ["A", "C", "G", "T"]


def load_model(model_dir, device):
    print(f"[INFO] Loading model from: {model_dir}")

    model = AutoModelForMaskedLM.from_pretrained(
        model_dir,
        trust_remote_code=True,
    )
    tokenizer = AutoTokenizer.from_pretrained(
        model_dir,
        trust_remote_code=True,
    )

    model.to(device)
    model.eval()
    return model, tokenizer


def get_vocab_indices(tokenizer):
    vocab = tokenizer.get_vocab()

    nts_lower = ["a", "c", "g", "t"]
    nts_upper = ["A", "C", "G", "T"]

    if all(x in vocab for x in nts_lower):
        nts = nts_lower
    elif all(x in vocab for x in nts_upper):
        nts = nts_upper
    else:
        raise ValueError("Tokenizer vocab does not contain A/C/G/T tokens.")

    vocab_indices = [vocab[x] for x in nts]
    nt_to_col = {x.upper(): i for i, x in enumerate(nts)}

    return vocab_indices, nt_to_col


def build_junction_indices(row, seq_length, flank_size):
    half = seq_length // 2
    sv_type = row["sv_type"]
    deletion_length = int(row["deletion_length"])

    if sv_type == "TIS_upstream_deletion":
        # Ref:
        # deletion is [half - deletion_length, half)
        # TIS/ATG anchor is at half
        #
        # Mut:
        # deleted upstream region is removed
        # new junction is at half
        ref_indices = (
            list(range(half - deletion_length - flank_size, half - deletion_length)) +
            list(range(half, half + flank_size))
        )

        mut_indices = (
            list(range(half - flank_size, half)) +
            list(range(half, half + flank_size))
        )

    else:
        # CDS deletion starts at center.
        # Ref compares left flank before deletion + right flank after deletion.
        # Mut compares the new junction at center.
        ref_indices = (
            list(range(half - flank_size, half)) +
            list(range(half + deletion_length, half + deletion_length + flank_size))
        )

        mut_indices = (
            list(range(half - flank_size, half)) +
            list(range(half, half + flank_size))
        )

    return ref_indices, mut_indices


def extract_indexed_sequence(seq, indices):
    seq = str(seq).upper()
    return "".join(seq[i] if 0 <= i < len(seq) else "N" for i in indices)


def score_sequence_pll(
    sequence,
    positions,
    model,
    tokenizer,
    device,
    mask_batch_size=16,
):
    """
    Calculate mean masked pseudo-log-likelihood across selected positions.

    For each evaluated position:
      1. Mask only that position.
      2. Predict the original nucleotide.
      3. Record its log-probability.
    """
    sequence = str(sequence).upper()

    encoding = tokenizer(
        sequence,
        truncation=False,
        padding=False,
        return_tensors="pt",
        return_attention_mask=False,
        return_token_type_ids=False,
        add_special_tokens=False,
    )

    original_ids = encoding["input_ids"][0]

    if original_ids.shape[0] != len(sequence):
        raise ValueError(
            f"Tokenizer length mismatch: sequence={len(sequence)}, "
            f"tokens={original_ids.shape[0]}"
        )

    valid_positions = []
    target_token_ids = []

    for pos in positions:
        if pos < 0 or pos >= len(sequence):
            continue

        nt = sequence[pos]
        if nt not in {"A", "C", "G", "T"}:
            continue

        valid_positions.append(pos)
        target_token_ids.append(int(original_ids[pos]))

    if not valid_positions:
        return np.nan

    position_scores = []

    for start in range(0, len(valid_positions), mask_batch_size):
        batch_positions = valid_positions[start:start + mask_batch_size]
        batch_targets = target_token_ids[start:start + mask_batch_size]

        # One copy of the sequence per masked position.
        masked_ids = original_ids.unsqueeze(0).repeat(
            len(batch_positions), 1
        )

        row_indices = torch.arange(len(batch_positions))
        col_indices = torch.tensor(batch_positions, dtype=torch.long)

        masked_ids[row_indices, col_indices] = tokenizer.mask_token_id
        masked_ids = masked_ids.to(device)

        with torch.inference_mode():
            logits = model(input_ids=masked_ids).logits.float()

        row_indices_device = torch.arange(
            len(batch_positions),
            device=device,
        )
        col_indices_device = torch.tensor(
            batch_positions,
            device=device,
            dtype=torch.long,
        )
        target_ids_device = torch.tensor(
            batch_targets,
            device=device,
            dtype=torch.long,
        )

        selected_logits = logits[
            row_indices_device,
            col_indices_device,
            :,
        ]

        log_probs = torch.log_softmax(selected_logits, dim=-1)

        target_log_probs = log_probs[
            row_indices_device,
            target_ids_device,
        ]

        position_scores.extend(
            target_log_probs.detach().cpu().tolist()
        )

        del masked_ids, logits, selected_logits, log_probs

    return float(np.mean(position_scores))


def score_all_sequences(sequences, indices_list, model, tokenizer, device, batch_size, desc):
    all_scores = []

    for i in tqdm(range(0, len(sequences), batch_size), desc=desc):
        batch_seqs = sequences[i:i + batch_size]
        batch_indices = indices_list[i:i + batch_size]

        scores = score_batch_sequences(
            sequences=batch_seqs,
            indices_list=batch_indices,
            model=model,
            tokenizer=tokenizer,
            device=device,
        )

        all_scores.extend(scores)

    return np.array(all_scores, dtype=float)


def score_all_sequences_pll(
    sequences,
    indices_list,
    model,
    tokenizer,
    device,
    mask_batch_size,
    desc,
):
    scores = []

    iterator = zip(sequences, indices_list)

    for sequence, positions in tqdm(
        iterator,
        total=len(sequences),
        desc=desc,
    ):
        score = score_sequence_pll(
            sequence=sequence,
            positions=positions,
            model=model,
            tokenizer=tokenizer,
            device=device,
            mask_batch_size=mask_batch_size,
        )
        scores.append(score)

    return np.asarray(scores, dtype=float)


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--input", required=True)
    ap.add_argument("--output", default="synthetic_deletion_scored.tsv")
    ap.add_argument("--report", default="synthetic_deletion_score_position_report.tsv")
    ap.add_argument("--model", required=True)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--batch_size", type=int, default=4)
    ap.add_argument("--seq_length", type=int, default=8192)
    ap.add_argument("--flank_size", type=int, default=10)

    args = ap.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    required = {"ref_sequence", "mut_sequence", "sv_type", "deletion_length"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input missing columns: {missing}")

    ref_indices = []
    mut_indices = []
    ref_scored_seq = []
    mut_scored_seq = []

    for _, row in df.iterrows():
        r_idx, m_idx = build_junction_indices(
            row,
            seq_length=args.seq_length,
            flank_size=args.flank_size,
        )

        ref_indices.append(r_idx)
        mut_indices.append(m_idx)

        ref_s = extract_indexed_sequence(row["ref_sequence"], r_idx)
        mut_s = extract_indexed_sequence(row["mut_sequence"], m_idx)

        ref_scored_seq.append(ref_s)
        mut_scored_seq.append(mut_s)

    same_scored_seq = [
        r == m for r, m in zip(ref_scored_seq, mut_scored_seq)
    ]

    report_df = df.drop(columns=["ref_sequence", "mut_sequence"]).copy()
    report_df["ref_score_indices"] = [",".join(map(str, x)) for x in ref_indices]
    report_df["mut_score_indices"] = [",".join(map(str, x)) for x in mut_indices]
    report_df["ref_scored_sequence"] = ref_scored_seq
    report_df["mut_scored_sequence"] = mut_scored_seq
    report_df["ref_mut_scored_sequence_identical"] = same_scored_seq

    n_bad = int((~pd.Series(same_scored_seq)).sum())
    print(f"[CHECK] ref/mut scored sequence identical: {len(df) - n_bad}/{len(df)}")
    print(f"[CHECK] mismatched rows: {n_bad}")

    if n_bad > 0:
        print("[WARNING] Some ref/mut scored sequences are not identical.")
        print(report_df.loc[
            ~report_df["ref_mut_scored_sequence_identical"],
            ["sv_id", "sv_type", "deletion_length", "ref_scored_sequence", "mut_scored_sequence"]
        ].head(10))

    report_df.to_csv(args.report, sep="\t", index=False)
    print(f"[INFO] Wrote report: {args.report}")

    model, tokenizer = load_model(args.model, args.device)

    ref_scores = score_all_sequences_pll(
        sequences=df["ref_sequence"].astype(str).tolist(),
        indices_list=ref_indices,
        model=model,
        tokenizer=tokenizer,
        device=args.device,
        mask_batch_size=args.batch_size,
        desc="Reference PLL",
    )

    mut_scores = score_all_sequences_pll(
        sequences=df["mut_sequence"].astype(str).tolist(),
        indices_list=mut_indices,
        model=model,
        tokenizer=tokenizer,
        device=args.device,
        mask_batch_size=args.batch_size,
        desc="Mutated PLL",
    )

    out = df.drop(columns=["ref_sequence", "mut_sequence"]).copy()

    out["ref_junction_logprob"] = ref_scores
    out["mut_junction_logprob"] = mut_scores
    out["delta_pll_mut_minus_ref"] = mut_scores - ref_scores
    out["flank_size"] = args.flank_size
    out["scored_positions"] = args.flank_size * 2
    out["ref_mut_scored_sequence_identical"] = same_scored_seq

    out.to_csv(args.output, sep="\t", index=False)

    print(f"\nWrote: {args.output}")
    print()
    print(out.groupby("sv_type")["score_drop_ref_minus_mut"].describe())


if __name__ == "__main__":
    main()