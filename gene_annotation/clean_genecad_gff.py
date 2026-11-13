#!/usr/bin/env python3
"""Create a GeneCAD-compatible protein-coding, one-transcript-per-gene GFF3."""

from __future__ import annotations

import argparse
import gzip
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


TRANSCRIPT_TYPES = {"mRNA", "transcript"}
CHILD_TYPES = {"exon", "CDS", "five_prime_UTR", "three_prime_UTR"}


def open_text(path: Path, mode: str) -> TextIO:
    text_mode = mode if "t" in mode else mode + "t"
    kwargs = {
        "encoding": "utf-8",
        "errors": "replace" if "r" in mode else "strict",
        "newline": "",
    }
    if path.suffix.lower() == ".gz":
        return gzip.open(path, text_mode, **kwargs)
    return path.open(text_mode, **kwargs)


def parse_attributes(text: str) -> dict[str, str]:
    result = {}
    for item in text.rstrip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key.strip()] = value.strip()
    return result


def format_attributes(attributes: dict[str, str]) -> str:
    return ";".join(f"{key}={value}" for key, value in attributes.items())


@dataclass
class Record:
    fields: list[str]
    attributes: dict[str, str]
    order: int

    @property
    def feature_type(self) -> str:
        return self.fields[2]

    @property
    def feature_id(self) -> str | None:
        return self.attributes.get("ID")

    @property
    def parents(self) -> list[str]:
        return [x for x in self.attributes.get("Parent", "").split(",") if x]

    @property
    def span(self) -> int:
        return int(self.fields[4]) - int(self.fields[3]) + 1


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    records: list[Record] = []
    headers: list[str] = []
    replacement_count = 0
    with open_text(args.input, "r") as handle:
        for order, line in enumerate(handle):
            replacement_count += line.count("\ufffd")
            # Malformed bytes in Araport11 occur in free-text attributes (for
            # example an cp1252 en dash). Replace them with a safe ASCII dash;
            # feature coordinates and the Parent/ID hierarchy are unchanged.
            line = line.replace("\ufffd", "-")
            if line.startswith("#"):
                headers.append(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            records.append(Record(fields, parse_attributes(fields[8]), order))

    genes = {r.feature_id: r for r in records if r.feature_type == "gene" and r.feature_id}
    transcripts = {
        r.feature_id: r
        for r in records
        if r.feature_type in TRANSCRIPT_TYPES and r.feature_id and r.parents
    }

    tx_to_gene: dict[str, str] = {}
    gene_to_tx: dict[str, list[str]] = defaultdict(list)
    for tx_id, tx in transcripts.items():
        gene_id = tx.parents[0]
        if gene_id in genes:
            tx_to_gene[tx_id] = gene_id
            gene_to_tx[gene_id].append(tx_id)

    children_by_tx: dict[str, list[Record]] = defaultdict(list)
    for record in records:
        if record.feature_type not in CHILD_TYPES:
            continue
        for parent in record.parents:
            if parent in transcripts:
                children_by_tx[parent].append(record)

    selected_tx: dict[str, str] = {}
    for gene_id, tx_ids in gene_to_tx.items():
        candidates = []
        for tx_id in tx_ids:
            cds_length = sum(r.span for r in children_by_tx[tx_id] if r.feature_type == "CDS")
            if cds_length > 0:
                candidates.append((cds_length, transcripts[tx_id].span, tx_id))
        if candidates:
            selected_tx[gene_id] = max(candidates)[2]

    selected_tx_ids = set(selected_tx.values())
    selected_gene_ids = set(selected_tx)
    output_records: list[Record] = []
    child_segment_number: dict[tuple[str, str], int] = defaultdict(int)

    for record in records:
        keep = False
        if record.feature_type == "gene":
            keep = record.feature_id in selected_gene_ids
        elif record.feature_type in TRANSCRIPT_TYPES:
            keep = record.feature_id in selected_tx_ids
            if keep:
                record.fields[2] = "mRNA"
                record.attributes["longest"] = "1"
        elif record.feature_type in CHILD_TYPES:
            kept_parents = [p for p in record.parents if p in selected_tx_ids]
            keep = bool(kept_parents)
            if keep:
                parent = kept_parents[0]
                record.attributes["Parent"] = parent
                # Many valid GFF3 producers reuse one CDS/exon ID across all
                # segments of a transcript. GeneCAD treats feature_id as part
                # of a primary key, so every genomic segment needs a unique ID.
                counter_key = (parent, record.feature_type)
                child_segment_number[counter_key] += 1
                segment = child_segment_number[counter_key]
                record.attributes["ID"] = (
                    f"{record.feature_type}:{parent}:{record.fields[0]}:"
                    f"{record.fields[3]}:{record.fields[4]}:{record.fields[6]}:{segment}"
                )
        if keep:
            record.fields[8] = format_attributes(record.attributes)
            output_records.append(record)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.output, "w") as handle:
        if headers:
            handle.writelines(headers)
        else:
            handle.write("##gff-version 3\n")
        for record in sorted(output_records, key=lambda r: r.order):
            handle.write("\t".join(record.fields) + "\n")

    counts = defaultdict(int)
    for record in output_records:
        counts[record.feature_type] += 1
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print(f"Malformed UTF-8 bytes replaced: {replacement_count:,}")
    print(f"Protein-coding genes retained: {len(selected_gene_ids):,}")
    print(f"Representative mRNAs retained: {len(selected_tx_ids):,}")
    print("Output feature counts:")
    for feature_type in sorted(counts):
        print(f"  {feature_type}: {counts[feature_type]:,}")

    if len(selected_gene_ids) != len(selected_tx_ids):
        raise RuntimeError("Expected exactly one selected transcript per retained gene")


if __name__ == "__main__":
    main()
