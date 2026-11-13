# Bryophyte pangenome selected-contig evaluation

Manifest: 123 species (82 Mosses, 37 Liverworts, 4 Hornworts). Prediction: BryoCAD-Small / 4bryo; filename `<short_name>_GeneCAD_final.gff.gz`. Reference: `<short_name>.selected.gff3.gz`. Only manifest species are included; the three additional references in the reference directory are not included. Inputs are recorded in inputs.tsv.

Same evaluation definitions and `--structure-policy exon` as the model-organism comparison. Original CDS coordinates are preserved. No selective gene filtering, FASTA, motif, or BUSCO. Both inputs must refer to the same assembly and contig coordinate system. `--seqid-policy union` retains features unique to either input as FP/FN; it does not silently select shared seqids. Per-input differences are recorded in individual audit.json and seqid_audit_summary.tsv. Absence of a seqid in a feature table does not imply absence from the underlying FASTA; no assembly identity can be proven from GFF alone.

## Outputs

- individual/<short_name>/: metrics.tsv, evaluation.txt, audit.json, run.log, completed.json.
- pangenome.metrics.tsv: 2460 rows (123 species × 20 metrics), including taxonomy, counts, precision/recall/F1, availability and reasons.
- pangenome.f1.wide.tsv: 123 species rows, 20 metric columns.
- metric_availability.tsv: available/NA species counts and reasons per metric.
- plots/*_PR.pdf and .svg: one metric per figure; one point per species; all species on the same axes, no species/lineage panels. Legend counts are evaluable PR species.
- plots/*_F1_bar.pdf and .svg: one metric per figure, species grouped by division/order/name; NA explicitly labelled.
- plots/all_metrics_PR.pdf: 20 pages, one metric per page.
- plots/division_summary.tsv: species counts and species-equal mean/median F1. These are not count-pooled micro-F1.
- Only PDF and SVG plot formats are generated.

Colours match plot_recovery.R's default ggplot discrete hue mapping: Hornworts #F8766D, Liverworts #00BA38, Mosses #619CFF. F1 curves are analytical isocontours, not threshold PR curves. The target subset is approximately 50 Mb of whole contigs, not a guarantee of equal assembly span or equal gene counts among species.

## Availability

All 123 species have CDS segment/transcript, coding intron/chain, CDS base and coding site metrics. Full exon/UTR comparisons are available for 83 species; 40 references lack sufficient full-structure evidence under the evaluator's global availability rule. All four hornwort references fall in that unavailable group. No genes are removed to manufacture full-structure availability.

CDS/locus is available for 65 species: 58 lack explicit gene grouping for at least one coding transcript. IDs such as g12820.t1 are not treated as proof of a gene hierarchy. Full-structure locus_ic additionally requires full-structure evidence, leaving 26 species. NA is not zero; predicted descriptive counts remain in individual metric outputs when computable.

## Reproduce

```bash
python3 run_pangenome_evaluations.py \
  --manifest ../bryophytes_pangenome_manifest.tsv \
  --predictions ../bryophytes_pangenome_bryocad_4bryo_50Mb \
  --references ../selected_bryophytes_50Mb \
  --structure-policy exon --workers 4 --outdir .
Rscript plot_pangenome.R pangenome.metrics.tsv plots
```

Successful evaluations resume only if input paths, sizes, modification timestamps, evaluator SHA256 and policies match completed.json, and the expected files exist. Failures are listed in run_status.tsv; other species continue, but the command exits nonzero if any failed. Use --workers 1 for serial execution. Rerunning plotting regenerates all plots.
