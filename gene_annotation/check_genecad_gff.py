#!/usr/bin/env python3

import argparse
import gzip
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def open_text(path):
    path = Path(path)
    kwargs = {"encoding": "utf-8", "errors": "replace", "newline": ""}
    if path.suffix == ".gz":
        return gzip.open(path, "rt", **kwargs)
    return path.open("rt", **kwargs)


def parse_attributes(text):
    result = {}
    for field in text.rstrip(";").split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            result[key] = value
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", required=True)
    parser.add_argument("--fasta", required=True)
    args = parser.parse_args()

    print("Loading FASTA...")
    with open_text(args.fasta) as handle:
        genome = {
            record.id: str(record.seq).upper()
            for record in SeqIO.parse(handle, "fasta")
        }

    feature_counts = Counter()
    ids = Counter()
    genes = {}
    transcripts = {}
    children = defaultdict(list)

    errors = Counter()
    examples = defaultdict(list)
    replacement_count = 0

    with open_text(args.gff) as handle:
        for line_number, line in enumerate(handle, 1):
            replacement_count += line.count("\ufffd")
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                errors["not_9_columns"] += 1
                continue

            seqid, source, feature, start, end, score, strand, phase, attrs = fields

            try:
                start = int(start)
                end = int(end)
            except ValueError:
                errors["invalid_coordinate"] += 1
                continue

            attributes = parse_attributes(attrs)
            feature_id = attributes.get("ID")
            parents = [
                value
                for value in attributes.get("Parent", "").split(",")
                if value
            ]

            feature_counts[feature] += 1

            if seqid not in genome:
                errors["seqid_missing_from_fasta"] += 1
                continue

            if start < 1 or end < start or end > len(genome[seqid]):
                errors["coordinate_out_of_range"] += 1

            if strand not in {"+", "-"}:
                errors["invalid_strand"] += 1

            if feature_id:
                ids[feature_id] += 1

            record = {
                "seqid": seqid,
                "feature": feature,
                "start": start,
                "end": end,
                "strand": strand,
                "phase": phase,
                "id": feature_id,
                "parents": parents,
                "line": line_number,
            }

            if feature == "gene" and feature_id:
                genes[feature_id] = record
            elif feature in {"mRNA", "transcript"} and feature_id:
                transcripts[feature_id] = record

            for parent in parents:
                children[parent].append(record)

    duplicate_ids = {
        feature_id: count
        for feature_id, count in ids.items()
        if count > 1
    }

    errors["duplicate_ids"] = len(duplicate_ids)

    transcripts_per_gene = Counter()

    for transcript_id, transcript in transcripts.items():
        if len(transcript["parents"]) != 1:
            errors["transcript_parent_count_not_1"] += 1
            continue

        gene_id = transcript["parents"][0]
        transcripts_per_gene[gene_id] += 1

        gene = genes.get(gene_id)

        if gene is None:
            errors["transcript_parent_gene_missing"] += 1
            continue

        if (
            transcript["seqid"] != gene["seqid"]
            or transcript["strand"] != gene["strand"]
            or transcript["start"] < gene["start"]
            or transcript["end"] > gene["end"]
        ):
            errors["transcript_outside_gene"] += 1

    cds_statistics = Counter()
    splice_statistics = Counter()

    for transcript_id, transcript in transcripts.items():
        transcript_children = children.get(transcript_id, [])

        cds = sorted(
            [record for record in transcript_children if record["feature"] == "CDS"],
            key=lambda record: record["start"],
        )

        if not cds:
            errors["transcript_without_cds"] += 1
            continue

        for record in transcript_children:
            if (
                record["seqid"] != transcript["seqid"]
                or record["strand"] != transcript["strand"]
                or record["start"] < transcript["start"]
                or record["end"] > transcript["end"]
            ):
                errors["child_outside_transcript"] += 1

        if transcript["strand"] == "-":
            cds = list(reversed(cds))

        cds_sequence_parts = []

        for record in cds:
            sequence = genome[record["seqid"]][
                record["start"] - 1 : record["end"]
            ]

            if record["strand"] == "-":
                sequence = str(Seq(sequence).reverse_complement())

            cds_sequence_parts.append(sequence)

        cds_sequence = "".join(cds_sequence_parts)

        cds_statistics["transcripts_with_cds"] += 1

        if len(cds_sequence) % 3 == 0:
            cds_statistics["cds_length_multiple_of_3"] += 1
        else:
            cds_statistics["cds_length_not_multiple_of_3"] += 1

        if cds_sequence.startswith("ATG"):
            cds_statistics["start_ATG"] += 1
        else:
            cds_statistics["start_not_ATG"] += 1

        if cds_sequence[-3:] in {"TAA", "TAG", "TGA"}:
            cds_statistics["valid_stop"] += 1
        else:
            cds_statistics["invalid_stop"] += 1

        if len(cds_sequence) % 3 == 0:
            protein = str(Seq(cds_sequence).translate())
            internal = protein[:-1] if protein.endswith("*") else protein

            if "*" in internal:
                cds_statistics["internal_stop"] += 1
            else:
                cds_statistics["no_internal_stop"] += 1

        genomic_cds = sorted(cds, key=lambda record: record["start"])

        for left, right in zip(genomic_cds, genomic_cds[1:]):
            intron_start = left["end"] + 1
            intron_end = right["start"] - 1

            if intron_end < intron_start:
                errors["overlapping_cds_segments"] += 1
                continue

            intron = genome[transcript["seqid"]][
                intron_start - 1 : intron_end
            ]

            if transcript["strand"] == "-":
                intron = str(Seq(intron).reverse_complement())

            if len(intron) < 4:
                splice_statistics["short_intron"] += 1
                continue

            motif = f"{intron[:2]}-{intron[-2:]}"
            splice_statistics[motif] += 1

    print("\n===== Input encoding =====")
    print(f"malformed_UTF8_replacement_characters {replacement_count}")

    print("\n===== Feature counts =====")
    for key, value in feature_counts.most_common():
        print(f"{key:25s} {value}")

    print("\n===== Parent structure =====")
    print("genes:", len(genes))
    print("transcripts:", len(transcripts))
    print(
        "genes with exactly one transcript:",
        sum(value == 1 for value in transcripts_per_gene.values()),
    )
    print(
        "genes with multiple transcripts:",
        sum(value > 1 for value in transcripts_per_gene.values()),
    )

    print("\n===== Errors =====")
    for key, value in errors.most_common():
        print(f"{key:40s} {value}")

    print("\n===== CDS statistics =====")
    for key, value in cds_statistics.most_common():
        print(f"{key:40s} {value}")

    print("\n===== Most common CDS intron motifs =====")
    for motif, count in splice_statistics.most_common(15):
        print(f"{motif:15s} {count}")


if __name__ == "__main__":
    main()
