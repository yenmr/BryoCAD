#!/usr/bin/env python3
"""Score PlantCAD2/BryoCAD zero-shot motif recovery from local TSV input."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset
from tqdm import tqdm
from transformers import AutoModelForMaskedLM, AutoTokenizer


LOGGER = logging.getLogger("plantcad2_motif_recovery")
NUCLEOTIDES = ("A", "C", "G", "T")
TASK_RULES = {
    "donor": {"motifs": {"GT"}, "length": 2},
    "acceptor": {"motifs": {"AG"}, "length": 2},
    "tis": {"motifs": {"ATG"}, "length": 3},
    "tts": {"motifs": {"TAA", "TAG", "TGA"}, "length": 3},
}


def resolve_dtype(dtype_name: str, device: str) -> torch.dtype:
    """Resolve explicit/automatic precision for the selected device."""
    choices = {
        "float32": torch.float32,
        "float16": torch.float16,
        "bfloat16": torch.bfloat16,
    }
    if dtype_name != "auto":
        return choices[dtype_name]
    if not device.startswith("cuda") or not torch.cuda.is_available():
        return torch.float32
    device_obj = torch.device(device)
    major, _ = torch.cuda.get_device_capability(device_obj)
    if major >= 8 and torch.cuda.is_bf16_supported():
        return torch.bfloat16
    # Half-precision capability alone does not establish Triton kernel compatibility.
    # Keep pre-Ampere inference in FP32 unless the user explicitly selects FP16.
    return torch.float32


def require_device(device: str) -> torch.device:
    dev = torch.device(device)
    if dev.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError(f"CUDA device requested but CUDA is unavailable: {device}")
    return dev


def load_model(model_name: str, device: str, dtype_name: str):
    dev = require_device(device)
    dtype = resolve_dtype(dtype_name, device)
    if dev.type == "cpu" and dtype in {torch.float16, torch.bfloat16}:
        raise ValueError("float16/bfloat16 CPU inference is not supported; use --dtype float32")
    LOGGER.info("Loading %s on %s with dtype=%s", model_name, dev, dtype)
    kwargs = {"trust_remote_code": True, "torch_dtype": dtype}
    model = AutoModelForMaskedLM.from_pretrained(model_name, **kwargs)
    tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
    if tokenizer.mask_token_id is None:
        raise ValueError("Tokenizer has no mask token")
    model.to(device=dev, dtype=dtype)
    model.eval()
    return model, tokenizer, dev, dtype


def prepare_dataframe(df: pd.DataFrame, task: str) -> pd.DataFrame:
    required = {"sequence", "geneid", "motif_start", "motif_end", "motif"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing input columns: {sorted(missing)}")
    out = df.copy()
    out["sequence"] = out["sequence"].astype(str).str.strip().str.upper()
    out["motif"] = out["motif"].astype(str).str.strip().str.upper()
    out["motif_start"] = pd.to_numeric(out["motif_start"], errors="raise").astype(int)
    out["motif_end"] = pd.to_numeric(out["motif_end"], errors="raise").astype(int)
    expected_len = TASK_RULES[task]["length"]
    errors = []
    for idx, row in out.iterrows():
        start, end, sequence = row["motif_start"], row["motif_end"], row["sequence"]
        if start < 0 or end > len(sequence) or end <= start:
            errors.append(f"row {idx}: motif range {start}:{end} outside sequence length {len(sequence)}")
        elif end - start != expected_len:
            errors.append(f"row {idx}: motif length {end - start}, expected {expected_len}")
        elif sequence[start:end] != row["motif"]:
            errors.append(
                f"row {idx}: motif={row['motif']!r}, sequence range={sequence[start:end]!r}"
            )
        if len(errors) >= 10:
            break
    if errors:
        raise ValueError("Invalid recovery input:\n" + "\n".join(errors))
    return out


class MotifMaskDataset(Dataset):
    def __init__(self, df: pd.DataFrame, tokenizer):
        self.df = df.reset_index(drop=True)
        self.tokenizer = tokenizer

    def __len__(self):
        return len(self.df)

    def __getitem__(self, index: int):
        row = self.df.iloc[index]
        sequence = row["sequence"]
        start, end = int(row["motif_start"]), int(row["motif_end"])
        encoded = self.tokenizer(
            sequence,
            return_tensors="pt",
            padding=False,
            truncation=False,
            add_special_tokens=False,
            return_attention_mask=False,
            return_token_type_ids=False,
        )["input_ids"]
        if encoded.shape[1] != len(sequence):
            raise ValueError(
                "PlantCAD2 scoring requires one token per nucleotide: "
                f"sequence length={len(sequence)}, token length={encoded.shape[1]}"
            )
        encoded[0, start:end] = self.tokenizer.mask_token_id
        return encoded.squeeze(0)


def nucleotide_token_ids(tokenizer) -> list[int]:
    vocab = tokenizer.get_vocab()
    ids = []
    for base in NUCLEOTIDES:
        if base.lower() in vocab:
            ids.append(vocab[base.lower()])
        elif base in vocab:
            ids.append(vocab[base])
        else:
            raise KeyError(f"Tokenizer vocabulary has no token for {base}/{base.lower()}")
    return ids


def run_inference(model, tokenizer, loader, device, motif_len: int) -> np.ndarray:
    """Return A/C/G/T probabilities shaped [number of inputs, motif_len, 4]."""
    base_ids = torch.tensor(nucleotide_token_ids(tokenizer), device=device)
    batches = []
    for input_ids in tqdm(loader, desc="PlantCAD2 masked inference"):
        input_ids = input_ids.to(device)
        mask = input_ids.eq(tokenizer.mask_token_id)
        counts = mask.sum(dim=1)
        if not torch.all(counts == motif_len):
            raise ValueError(f"Expected {motif_len} mask tokens per input; got {counts.tolist()}")
        with torch.inference_mode():
            logits = model(input_ids=input_ids).logits
        selected = logits[mask].reshape(input_ids.shape[0], motif_len, -1)
        selected = selected.index_select(2, base_ids)
        batches.append(torch.softmax(selected.float(), dim=-1).cpu().numpy())
    if not batches:
        return np.empty((0, motif_len, 4), dtype=np.float32)
    return np.concatenate(batches, axis=0)


def build_results(df: pd.DataFrame, probs: np.ndarray, task: str) -> pd.DataFrame:
    motif_len = TASK_RULES[task]["length"]
    if probs.shape != (len(df), motif_len, 4):
        raise ValueError(
            f"Probability shape {probs.shape}, expected {(len(df), motif_len, 4)}"
        )
    predicted_array = np.asarray(NUCLEOTIDES)[probs.argmax(axis=2)]
    predicted = np.asarray(["".join(row) for row in predicted_array])
    references = df["motif"].astype(str).str.upper().to_numpy()
    reference_arrays = np.asarray([list(motif) for motif in references])
    token_accuracy = (predicted_array == reference_arrays).mean(axis=1)
    base_to_index = {base: i for i, base in enumerate(NUCLEOTIDES)}
    reference_probs = np.zeros((len(df), motif_len), dtype=np.float32)
    for i, motif in enumerate(references):
        for j, base in enumerate(motif):
            if base in base_to_index:
                reference_probs[i, j] = probs[i, j, base_to_index[base]]
    valid_motifs = TASK_RULES[task]["motifs"]

    # Preserve every identifying/input field except the 8192-bp sequence.
    result = df.drop(columns=["sequence"]).copy()
    result["reference_motif"] = references
    result["reference_is_canonical"] = [x in valid_motifs for x in references]
    result["predicted_motif"] = predicted
    result["token_accuracy"] = token_accuracy
    result["exact_recovery"] = predicted == references
    result["predicted_is_canonical"] = [x in valid_motifs for x in predicted]
    result["mean_reference_probability"] = reference_probs.mean(axis=1)
    result["joint_reference_probability"] = reference_probs.prod(axis=1)
    result["mean_prediction_probability"] = probs.max(axis=2).mean(axis=1)
    return result


def print_summary(result: pd.DataFrame) -> None:
    print("\nSummary")
    print(f"N = {len(result)}")
    print(f"Mean token accuracy = {result['token_accuracy'].mean():.6f}")
    print(f"Exact recovery = {result['exact_recovery'].mean():.6f}")
    print(
        "Predicted canonical fraction = "
        f"{result['predicted_is_canonical'].mean():.6f}"
    )
    print("\nRecovery by reference class")
    for is_canonical, group in result.groupby("reference_is_canonical", sort=False):
        label = "canonical" if is_canonical else "noncanonical"
        print(
            f"{label}: N={len(group)}, "
            f"token_accuracy={group['token_accuracy'].mean():.6f}, "
            f"exact_recovery={group['exact_recovery'].mean():.6f}, "
            f"predicted_canonical={group['predicted_is_canonical'].mean():.6f}"
        )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, type=Path)
    parser.add_argument("-o", "--output", required=True, type=Path)
    parser.add_argument("--task", required=True, choices=sorted(TASK_RULES))
    parser.add_argument(
        "--model", default="kuleshov-group/PlantCAD2-Small-l24-d0768"
    )
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument(
        "--dtype", choices=("auto", "float32", "float16", "bfloat16"), default="auto"
    )
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--num-workers", type=int, default=1)
    parser.add_argument("--max-samples", type=int)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    if args.batch_size < 1:
        raise ValueError("--batch-size must be at least 1")
    if args.num_workers < 0:
        raise ValueError("--num-workers cannot be negative")
    df = pd.read_csv(args.input, sep="\t")
    if args.max_samples is not None:
        if args.max_samples < 1:
            raise ValueError("--max-samples must be at least 1")
        if len(df) > args.max_samples:
            df = df.sample(n=args.max_samples, random_state=args.seed).reset_index(drop=True)
    df = prepare_dataframe(df, args.task)
    LOGGER.info("Input rows: %d; task: %s", len(df), args.task)
    model, tokenizer, device, dtype = load_model(args.model, args.device, args.dtype)
    LOGGER.info("Resolved dtype: %s", dtype)
    dataset = MotifMaskDataset(df, tokenizer)
    loader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        pin_memory=device.type == "cuda",
    )
    probs = run_inference(model, tokenizer, loader, device, TASK_RULES[args.task]["length"])
    result = build_results(df, probs, args.task)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    LOGGER.info("Output written to %s", args.output)
    print_summary(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
