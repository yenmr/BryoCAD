# Gene annotation evaluation scope

Figure 5A evaluates the model-organism predictions. Figure 5D evaluates **123 bryophyte pangenome assemblies using selected chromosomes/contigs**, not full genomes. The author confirms that both panels assess CDS and transcript structure **without UTRs**.

In the bundled evaluator these are separate from `full_structure/full_tx`:

- `section=CDS, metric=cds_ex`: exact CDS segment matching.
- `section=CDS, metric=tx`: transcript matching by the complete CDS-chain structure, excluding UTR boundaries.

Use “CDS segments” and “Coding transcripts (UTRs excluded)” as unambiguous labels. Do not label the latter “Full transcripts” or substitute the full-exon/UTR-sensitive metric. Existing tables retain other metrics as original data; the coding-only plotting script selects the two fields above explicitly.

## Region selection for Figure 5D

```bash
python gene_annotation/select_genome_fraction.py \
  --manifest genome_manifest.tsv \
  --target-bp 50000000 \
  --min-contig-bp 100000 \
  --fallback-min-contig-bp 50000 25000 12500 \
  --seed 42 --outdir selected_50Mb --resume --continue-on-error
```

The TSV manifest requires `short_name`, `fa_path`, and `gff_path`. Selection is annotation-independent: eligible whole contigs are ordered by a deterministic seeded hash and accumulated until the target is reached. The final contig is not trimmed, so selected size may exceed 50 Mb. If eligible sequence totals less than 50 Mb, all eligible contigs are retained. Fallback thresholds are used only when **zero** contigs meet the previous threshold, not merely when total eligible sequence is below 50 Mb.

`--resume` reuses validated completed outputs; `--continue-on-error` processes other assemblies after failures. Archive the resulting `selected_manifest.tsv`, per-genome contig selections and selection summary with the final evaluation to establish actual membership, lengths and failures. The 123-assembly scope is distinct from the 126-assembly motif/orthogroup resource.

The provided candidate metrics table contains 123 distinct assembly IDs with available CDS-chain transcript metrics. Selection outputs are not bundled here, so exact per-assembly coverage is not independently re-established by the metrics table alone.
