#!/usr/bin/env python3
"""
Binary Classification Evaluation Script
Outputs performance metrics in TSV format for easy copying, pasting, and appending.
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate binary classification metrics and output as TSV."
    )
    parser.add_argument(
        "-t",
        "--tsv_path",
        type=str,
        required=True,
        help="Path to the ground truth TSV file.",
    )
    parser.add_argument(
        "-p",
        "--pred_path",
        type=str,
        required=True,
        help="Path to the predictions CSV file.",
    )
    parser.add_argument(
        "--label_col",
        type=str,
        default="label",
        help="Column name for ground truth label in TSV (default: 'label').",
    )
    parser.add_argument(
        "--pred_col",
        type=str,
        default="probability_positive",
        help="Column name for predicted probability in CSV (default: 'probability_positive').",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.5,
        help="Decision threshold for binary classification (default: 0.5).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Optional output TSV file path.",
    )
    return parser.parse_args()


def load_and_validate_data(
    tsv_path: str, pred_path: str, label_col: str, pred_col: str
) -> tuple[np.ndarray, np.ndarray]:
    df_tsv = pd.read_csv(tsv_path, sep="\t", usecols=[label_col])
    df_pred = pd.read_csv(pred_path, usecols=[pred_col])

    if len(df_tsv) != len(df_pred):
        raise ValueError(
            f"Row count mismatch! TSV: {len(df_tsv)}, Pred CSV: {len(df_pred)}"
        )

    s_label = df_tsv[label_col]
    s_pred = df_pred[pred_col]

    valid_mask = ~(s_label.isna() | s_pred.isna())
    n_dropped = len(s_label) - valid_mask.sum()
    if n_dropped > 0:
        print(f"[Warning] Dropped {n_dropped} row(s) containing NaN values.", file=sys.stderr)
        s_label = s_label[valid_mask]
        s_pred = s_pred[valid_mask]

    y_true = s_label.to_numpy(dtype=int)
    y_score = s_pred.to_numpy(dtype=float)

    unique_labels = set(np.unique(y_true))
    if not unique_labels.issubset({0, 1}):
        raise ValueError(f"Ground truth labels must be binary {{0, 1}}, found: {unique_labels}")

    return y_true, y_score


def calculate_metrics(y_true: np.ndarray, y_score: np.ndarray, threshold: float) -> dict[str, float | int]:
    y_pred = (y_score >= threshold).astype(int)

    n_total = len(y_true)
    n_pos = int(np.sum(y_true == 1))
    n_neg = int(np.sum(y_true == 0))
    pos_pct = (n_pos / n_total) * 100 if n_total > 0 else 0.0
    neg_pct = (n_neg / n_total) * 100 if n_total > 0 else 0.0
    pos_neg_ratio = (n_pos / n_neg) if n_neg > 0 else float("inf")

    auroc = float(roc_auc_score(y_true, y_score))
    auprc = float(average_precision_score(y_true, y_score))

    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()

    acc = float(accuracy_score(y_true, y_pred))
    prec = float(precision_score(y_true, y_pred, zero_division=0))
    rec = float(recall_score(y_true, y_pred, zero_division=0))
    spec = float(tn / (tn + fp)) if (tn + fp) > 0 else 0.0
    f1 = float(f1_score(y_true, y_pred, zero_division=0))
    mcc = float(matthews_corrcoef(y_true, y_pred))

    return {
        "total_samples": n_total,
        "n_positive": n_pos,
        "n_negative": n_neg,
        "pos_ratio_pct": round(pos_pct, 4),
        "neg_ratio_pct": round(neg_pct, 4),
        "pos_neg_ratio": round(pos_neg_ratio, 4) if pos_neg_ratio != float("inf") else -1.0,
        "auroc": round(auroc, 6),
        "auprc": round(auprc, 6),
        "threshold": threshold,
        "TP": int(tp),
        "FP": int(fp),
        "TN": int(tn),
        "FN": int(fn),
        "accuracy": round(acc, 6),
        "precision": round(prec, 6),
        "recall": round(rec, 6),
        "specificity": round(spec, 6),
        "f1_score": round(f1, 6),
        "mcc": round(mcc, 6),
    }


def main():
    args = parse_args()

    try:
        y_true, y_score = load_and_validate_data(
            tsv_path=args.tsv_path,
            pred_path=args.pred_path,
            label_col=args.label_col,
            pred_col=args.pred_col,
        )
    except Exception as e:
        print(f"[Error] {e}", file=sys.stderr)
        sys.exit(1)

    metrics = calculate_metrics(y_true, y_score, threshold=args.threshold)

    # Build a single-line TSV record.
    headers = list(metrics.keys())
    values = [str(v) for v in metrics.values()]

    tsv_header_str = "\t".join(headers)
    tsv_value_str = "\t".join(values)

    # Print the header and values to the terminal.
    print("=== TSV Output (Header & Row) ===")
    print(tsv_header_str)
    print(tsv_value_str)

    # Write output to file.
    if args.output:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Write the header for a new file; otherwise append values.
        file_exists = out_path.exists()
        with open(out_path, "a", encoding="utf-8") as f:
            if not file_exists:
                f.write(tsv_header_str + "\n")
            f.write(tsv_value_str + "\n")
        print(f"\n[Info] Appended TSV row to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
