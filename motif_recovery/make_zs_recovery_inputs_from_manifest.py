#!/usr/bin/env python3
"""Generate PlantCAD2/BryoCAD zero-shot recovery inputs from FASTA + GFF3.

For each manifest row, this program keeps the longest CDS-bearing transcript
per gene and writes selected TIS, TTS, splice-donor, and splice-acceptor windows in
transcript orientation.  Model-space motif positions are intentionally fixed:

  TIS/TTS       sequence[4094:4097]
  donor/acceptor sequence[4095:4097]

GFF and reported genomic coordinates are 1-based inclusive.  motif_start and
motif_end are Python-style 0-based, end-exclusive indices in `sequence`.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


WINDOW = 8192
STOP_CODONS = {"TAA", "TAG", "TGA"}
TASKS = ("tis", "tts", "donor", "acceptor")


def open_text(path: str | Path):
    path = Path(path).expanduser()
    # Public GFF files occasionally mix UTF-8 with legacy Windows-1252 bytes
    # in free-text Note/description attributes (for example byte 0x96). IDs,
    # Parent values, coordinates, and biological sequences are ASCII, so
    # replacing only undecodable prose is preferable to aborting a species.
    kwargs = {"encoding": "utf-8", "errors": "replace", "newline": ""}
    return gzip.open(path, "rt", **kwargs) if path.suffix == ".gz" else path.open("r", **kwargs)


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn",
                                      "TGCAYRMKSWVHDBNtgcayrmkswvhdbn"))[::-1].upper()


def parse_attrs(text: str) -> dict[str, str]:
    out = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif len(item.split(None, 1)) == 2:  # limited GTF compatibility
            key, value = item.split(None, 1)
            value = value.strip('"')
        else:
            # Some otherwise usable GFF3 files contain an unescaped semicolon
            # inside a free-text Note/description.  The resulting bare fragment
            # is not an attribute and is irrelevant to ID/Parent parsing.
            continue
        out[key.strip()] = value.strip()
    return out


def attr_id(attrs: dict[str, str], *keys: str) -> str | None:
    for key in keys:
        if attrs.get(key):
            return attrs[key].split(",", 1)[0]
    return None


@dataclass
class Transcript:
    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]] = field(default_factory=list)
    cds: list[tuple[int, int]] = field(default_factory=list)
    start_codons: list[tuple[int, int]] = field(default_factory=list)
    stop_codons: list[tuple[int, int]] = field(default_factory=list)

    @property
    def cds_length(self) -> int:
        return sum(end - start + 1 for start, end in merge_intervals(self.cds))

    @property
    def exon_length(self) -> int:
        return sum(end - start + 1 for start, end in merge_intervals(self.exons))


def merge_intervals(items: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(set(items)):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(x[0], x[1]) for x in merged]


def read_gff(path: str | Path, log: logging.Logger):
    transcripts: dict[str, Transcript] = {}
    transcript_gene: dict[str, str] = {}
    transcript_meta: dict[str, tuple[str, str]] = {}
    pending: list[tuple[str, str, int, int, str, str]] = []
    gene_ids: set[str] = set()

    with open_text(path) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                log.warning("Skipping malformed GFF line %d", line_no)
                continue
            chrom, _, feature, start_s, end_s, _, strand, _, attrs_s = cols
            if strand not in {"+", "-"}:
                continue
            try:
                start, end = int(start_s), int(end_s)
            except ValueError:
                continue
            attrs = parse_attrs(attrs_s)
            feature_l = feature.lower()
            if feature_l == "gene" or feature_l.endswith("_gene"):
                gid = attr_id(attrs, "ID", "gene_id")
                if gid:
                    gene_ids.add(gid)
                continue
            if feature_l in {"mrna", "transcript", "primary_transcript"} or feature_l.endswith("rna"):
                tid = attr_id(attrs, "ID", "transcript_id")
                gid = attr_id(attrs, "Parent", "gene_id")
                if tid:
                    gid = gid or tid
                    transcript_gene[tid] = gid
                    transcript_meta[tid] = (chrom, strand)
                    transcripts.setdefault(tid, Transcript(tid, gid, chrom, strand))
                continue
            if feature_l not in {"exon", "cds", "start_codon", "stop_codon"}:
                continue
            parents = attrs.get("Parent") or attrs.get("transcript_id")
            if not parents:
                continue
            for tid in parents.split(","):
                pending.append((tid.strip(), feature_l, start, end, chrom, strand))

    for tid, feature, start, end, chrom, strand in pending:
        gid = transcript_gene.get(tid, tid)
        tx = transcripts.setdefault(tid, Transcript(tid, gid, chrom, strand))
        if (tx.chrom, tx.strand) != (chrom, strand):
            log.warning("Ignoring inconsistent child feature for transcript %s", tid)
            continue
        target = {"exon": tx.exons, "cds": tx.cds,
                  "start_codon": tx.start_codons,
                  "stop_codon": tx.stop_codons}[feature]
        target.append((start, end))

    by_gene: dict[str, list[Transcript]] = defaultdict(list)
    for tx in transcripts.values():
        if tx.cds:
            by_gene[tx.gene_id].append(tx)
    return by_gene, gene_ids


def choose_longest(by_gene: dict[str, list[Transcript]]) -> list[Transcript]:
    selected = []
    for gene_id in sorted(by_gene):
        # min with negative lengths gives deterministic lexical-ID tie breaking.
        selected.append(min(by_gene[gene_id],
                            key=lambda tx: (-tx.cds_length, -tx.exon_length,
                                            tx.transcript_id)))
    return selected


def load_fasta(path: str | Path, wanted: set[str]) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current in wanted:
                    sequences[current] = "".join(chunks).upper()
                current = line[1:].strip().split()[0]
                chunks = []
            elif current in wanted:
                chunks.append(line.strip())
        if current in wanted:
            sequences[current] = "".join(chunks).upper()
    return sequences


def extract_window(genome: str, motif_start: int, motif_end: int, strand: str,
                   input_start: int) -> str | None:
    """Extract around a 1-based inclusive motif at a fixed oriented index."""
    motif_start0 = motif_start - 1
    motif_end0 = motif_end  # exclusive in 0-based coordinates
    if strand == "+":
        left = motif_start0 - input_start
    else:
        # In reverse complement: input_start = WINDOW - (motif_end0 - left).
        left = motif_end0 - (WINDOW - input_start)
    right = left + WINDOW
    if left < 0 or right > len(genome):
        return None
    seq = genome[left:right]
    return revcomp(seq) if strand == "-" else seq.upper()


def motif_definitions(tx: Transcript, tasks=TASKS, cds_stop="included"):
    """TTS means translation termination; CDS fallback convention is explicit."""
    cds = merge_intervals(tx.cds)
    if not cds:
        return
    cds_left, cds_right = cds[0][0], cds[-1][1]
    starts = merge_intervals(tx.start_codons)
    stops = merge_intervals(tx.stop_codons)
    if tx.strand == "+":
        tis = starts[0] if starts else (cds_left, cds_left + 2)
        tts = stops[-1] if stops else (cds_right - 2, cds_right)
    else:
        tis = starts[-1] if starts else (cds_right - 2, cds_right)
        tts = stops[0] if stops else (cds_left, cds_left + 2)
    if not stops and cds_stop == "excluded":
        tts = (cds_right + 1, cds_right + 3) if tx.strand == "+" else (cds_left - 3, cds_left - 1)
    # Codons split across an intron are uncommon but valid. They cannot be
    # represented by one contiguous genomic range/window, so log them as an
    # invalid 3-bp coordinate later rather than silently joining two loci.
    if "tis" in tasks:
        yield "tis", tis[0], tis[1], 4094, "start_codon" if starts else "cds_fallback"
    if "tts" in tasks:
        yield "tts", tts[0], tts[1], 4094, "stop_codon" if stops else "cds_fallback_" + cds_stop
    if not ({"donor", "acceptor"} & set(tasks)):
        return

    blocks = merge_intervals(tx.exons) if tx.exons else cds
    source = "exon" if tx.exons else "cds_fallback"
    for left_block, right_block in zip(blocks, blocks[1:]):
        intron_left = left_block[1] + 1
        intron_right = right_block[0] - 1
        if intron_right - intron_left + 1 < 2:
            continue
        if tx.strand == "+":
            yield "donor", intron_left, intron_left + 1, 4095, source
            yield "acceptor", intron_right - 1, intron_right, 4095, source
        else:
            yield "donor", intron_right - 1, intron_right, 4095, source
            yield "acceptor", intron_left, intron_left + 1, 4095, source


def canonical(task: str, motif: str) -> bool:
    return motif == {"tis": "ATG", "donor": "GT", "acceptor": "AG"}.get(task) \
        if task != "tts" else motif in STOP_CODONS


OUTPUT_FIELDS = ["sequence", "geneid", "transcript_id", "chromosome", "start",
                 "end", "direction", "genome_range", "motif_start", "motif_end",
                 "motif", "is_canonical", "coordinate_source"]


def process_species(row: dict[str, str], out_root: Path, overwrite: bool,
                    log: logging.Logger, tasks=TASKS, cds_stop="included") -> dict[str, object]:
    short_name = row["short_name"]
    outdir = out_root / short_name
    outputs = {task: outdir / f"{task}_recovery.tsv" for task in tasks}
    if not overwrite and any(p.exists() for p in outputs.values()):
        log.info("Skipping %s: output exists (use --overwrite)", short_name)
        return {"short_name": short_name, "status": "skipped_existing"}
    outdir.mkdir(parents=True, exist_ok=True)

    by_gene, all_gene_ids = read_gff(row["gff_path"], log)
    selected = choose_longest(by_gene)
    wanted_contigs = {tx.chrom for tx in selected}
    fasta = load_fasta(row["fa_path"], wanted_contigs)
    writers = {}
    handles = []
    for task, path in outputs.items():
        handle = path.open("w", newline="")
        handles.append(handle)
        writers[task] = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writers[task].writeheader()

    counts = defaultdict(Counter)
    skipped = Counter()
    processed_rows = []
    try:
        for tx in selected:
            genome = fasta.get(tx.chrom)
            if genome is None:
                skipped["missing_fasta_contig"] += 1
                continue
            processed_rows.append({"geneid": tx.gene_id, "transcript_id": tx.transcript_id,
                                   "chromosome": tx.chrom, "strand": tx.strand,
                                   "cds_length": tx.cds_length,
                                   "exon_length": tx.exon_length,
                                   "splice_boundary_source": "exon" if tx.exons else "cds_fallback"})
            for definition in motif_definitions(tx, tasks, cds_stop):
                task, start, end, input_start, *extra = definition
                if task not in tasks:
                    continue
                expected_size = 3 if task in {"tis", "tts"} else 2
                if start < 1 or end < start or end - start + 1 != expected_size:
                    skipped[f"{task}:invalid_coordinate"] += 1
                    continue
                seq = extract_window(genome, start, end, tx.strand, input_start)
                if seq is None:
                    skipped[f"{task}:contig_boundary"] += 1
                    continue
                input_end = input_start + (end - start + 1)
                motif = seq[input_start:input_end]
                is_can = canonical(task, motif)
                counts[(task, tx.strand)][motif] += 1
                writers[task].writerow({
                    "sequence": seq, "geneid": tx.gene_id,
                    "transcript_id": tx.transcript_id, "chromosome": tx.chrom,
                    "start": start, "end": end, "direction": tx.strand,
                    "genome_range": f"{tx.chrom}:{start}-{end}",
                    "motif_start": input_start, "motif_end": input_end,
                    "motif": motif, "is_canonical": is_can,
                    "coordinate_source": extra[0] if extra else "",
                })
    finally:
        for handle in handles:
            handle.close()

    with (outdir / "processed_genes.tsv").open("w", newline="") as handle:
        fields = ["geneid", "transcript_id", "chromosome", "strand", "cds_length",
                  "exon_length", "splice_boundary_source"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(processed_rows)

    qc_fields = ["motif_type", "strand", "total", "canonical_count",
                 "noncanonical_count", "canonical_fraction", "top_observed_motifs"]
    qc_rows = []
    for task in tasks:
        for strand in ("+", "-", "all"):
            observed = (counts[(task, strand)] if strand != "all" else
                        counts[(task, "+")] + counts[(task, "-")])
            total = sum(observed.values())
            can_count = sum(n for motif, n in observed.items() if canonical(task, motif))
            fraction = can_count / total if total else 0.0
            qc_rows.append({"motif_type": task, "strand": strand, "total": total,
                            "canonical_count": can_count,
                            "noncanonical_count": total - can_count,
                            "canonical_fraction": f"{fraction:.6f}",
                            "top_observed_motifs": ",".join(f"{m}:{n}" for m, n in observed.most_common(10))})
            if strand == "all":
                log.info("%s %-8s total=%d canonical=%d fraction=%.4f",
                         short_name, task, total, can_count, fraction)
                if total and fraction == 1.0:
                    log.warning("%s %s is exactly 100%% canonical; verify annotation/QC", short_name, task)

    with (outdir / "motif_qc_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=qc_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(qc_rows)
    with (outdir / "alignment_debug_summary.tsv").open("w", newline="") as handle:
        fields = ["motif_type", "strand", "kmer", "count", "fraction",
                  "is_canonical"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for task in tasks:
            for strand in ("+", "-"):
                observed = counts[(task, strand)]
                total = sum(observed.values())
                for motif, count in observed.most_common():
                    writer.writerow({"motif_type": task, "strand": strand,
                                     "kmer": motif, "count": count,
                                     "fraction": f"{count / total:.6f}" if total else "0.000000",
                                     "is_canonical": canonical(task, motif)})
    with (outdir / "skipped_records.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["reason", "count"])
        writer.writerow(["genes_without_cds", max(0, len(all_gene_ids) - len(by_gene))])
        for reason, count in sorted(skipped.items()):
            writer.writerow([reason, count])

    result: dict[str, object] = {"short_name": short_name, "status": "completed",
                                "gff_genes": len(all_gene_ids),
                                "coding_genes_selected": len(selected)}
    for task in tasks:
        result[f"{task}_count"] = sum(counts[(task, "+")].values()) + sum(counts[(task, "-")].values())
    return result


ANGIOSPERM_DIVISIONS = {"Rosids", "Asterids", "Monocots", "Eudicots", "Magnoliids", "Basal Angio"}


def read_manifest(path: Path) -> list[dict[str, str]]:
    path = path.resolve()
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = set(reader.fieldnames or [])
        if not ({"fa_path", "genome_path"} & headers) or "gff_path" not in headers or not ({"short_name", "output"} & headers):
            raise ValueError("Need fa_path/genome_path, gff_path, and short_name/output columns")
        rows = list(reader)
    seen = set()
    for line_no, row in enumerate(rows, 2):
        if None in row:
            raise ValueError(f"Manifest line {line_no}: extra tab-separated fields")
        for target, alias in (("fa_path", "genome_path"), ("short_name", "output")):
            if row.get(target) and row.get(alias) and row[target].strip() != row[alias].strip():
                raise ValueError(f"Manifest line {line_no}: conflicting {target}/{alias}")
            row[target] = (row.get(target) or row.get(alias) or "").strip()
        for field in ("fa_path", "gff_path", "short_name"):
            if not (row.get(field) or "").strip():
                raise ValueError(f"Manifest line {line_no}: empty {field}")
        name = row["short_name"]
        if name in {".", ".."} or "/" in name or "\\" in name or name in seen:
            raise ValueError(f"Manifest line {line_no}: unsafe or duplicate short_name {name}")
        seen.add(name)
        for field in ("fa_path", "gff_path"):
            p = Path(row[field].strip()).expanduser()
            row[field] = str(p if p.is_absolute() else path.parent / p)
    return rows


def check_paths(rows, out_root):
    out_root.mkdir(parents=True, exist_ok=True)
    results = []
    for row in rows:
        for field in ("fa_path", "gff_path"):
            p = Path(row[field])
            status, detail = "ok", ""
            try:
                if not p.is_file():
                    status = "missing_or_not_file"
                elif p.stat().st_size == 0:
                    status = "empty"
                else:
                    with open_text(p) as handle:
                        handle.read(1)  # readability/header check, not a full integrity scan
            except OSError as exc:
                status, detail = "unreadable", str(exc)
            results.append(dict(short_name=row["short_name"], field=field, path=str(p), status=status, detail=detail))
    with (out_root / "path_check.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["short_name", "field", "path", "status", "detail"], delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    return all(r["status"] == "ok" for r in results)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--out", default="zs_recovery_inputs", type=Path)
    ap.add_argument("--species", nargs="*", help="Only these manifest short_name values")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--tasks", nargs="+", choices=TASKS, default=list(TASKS), help="Use --tasks tts for translation termination only")
    ap.add_argument("--angiosperms-only", action="store_true", help="Filter by manifest division")
    ap.add_argument("--check-only", action="store_true", help="Write path_check.tsv without generating sequences")
    ap.add_argument("--cds-stop", choices=("included", "excluded"), default="included", help="Fallback only: whether annotated CDS includes the stop codon; explicit stop_codon takes priority")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    log = logging.getLogger("zs_recovery")
    try:
        rows = read_manifest(args.manifest.expanduser())
    except (OSError, ValueError) as exc:
        log.error("Cannot read manifest: %s", exc)
        return 2
    if args.angiosperms_only:
        if any("division" not in row for row in rows):
            log.error("--angiosperms-only requires division column")
            return 2
        rows = [r for r in rows if r["division"] in ANGIOSPERM_DIVISIONS]
    if args.species:
        wanted = set(args.species)
        unknown = wanted - {r["short_name"] for r in rows}
        if unknown:
            log.error("Unknown/unselected species: %s", ", ".join(sorted(unknown)))
            return 2
        rows = [row for row in rows if row["short_name"] in wanted]
    if not rows:
        log.error("No manifest rows selected")
        return 2
    args.out.expanduser().mkdir(parents=True, exist_ok=True)
    paths_ok = check_paths(rows, args.out.expanduser())
    if args.check_only or not paths_ok:
        log.info("Path report: %s", args.out.expanduser() / "path_check.tsv")
        if not paths_ok:
            log.error("Some input paths are missing/unreadable on this host; no sequences generated")
        return 0 if paths_ok else 1
    summary = []
    for row in rows:
        try:
            for field in ("fa_path", "gff_path"):
                if not Path(row[field]).expanduser().is_file():
                    raise FileNotFoundError(row[field])
            summary.append(process_species(row, args.out.expanduser(), args.overwrite, log, tuple(dict.fromkeys(args.tasks)), args.cds_stop))
        except Exception as exc:
            log.exception("Failed species %s: %s", row.get("short_name", "?"), exc)
            summary.append({"short_name": row.get("short_name", "?"), "status": f"failed: {exc}"})
    fields = sorted({key for item in summary for key in item})
    with (args.out.expanduser() / "batch_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)
    return 1 if any(str(x["status"]).startswith("failed") for x in summary) else 0


if __name__ == "__main__":
    sys.exit(main())
