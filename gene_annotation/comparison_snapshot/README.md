# Model organism comparison

18 comparisons from comparison_list.snapshot.tsv: six observed model/head combinations × Ath, Mpo, Ppa. Helixer is a baseline without a fine-tuning head. This is not a complete model × head factorial design.

All comparisons use `--structure-policy exon`. When exon and UTR conflict, exon+CDS determines terminal UTR geometry; CDS is preserved. If exon is absent but CDS+UTR exists, exons are reconstructed. No selective gene filtering. No FASTA, sequence motif scoring, or BUSCO.

The bundled evaluator adds `--seqid-policy strict|union` (default strict). The batch runner explicitly uses union: sequences carrying features in only one input are retained, producing FP/FN. This assumes reference and prediction use the same assembly and sequence names. Sequence-name differences are recorded per comparison in audit.json. A feature table alone cannot prove that the underlying assembly was identical; a completely different naming scheme must be resolved upstream.

## Outputs

- `individual/<Model>__<Head>__<Species>/`: metrics.tsv, evaluation.txt, audit.json, run.log.
- `comparison.metrics.tsv`: long matrix, 20 metrics per comparison, with Model/Head/Species/Policy metadata.
- `comparison.f1.wide.tsv`: one comparison per row, one metric per column.
- `plots/overview_PR.*`: main coding comparison, species columns, metric rows, model colour, head shape.
- `plots/<section>_PR.*`: all metrics, model columns, metric rows, species colour, head shape.
- `plots/<section>_F1_bar.*`: F1 bars, same model/metric panel layout, head on x axis, species fill.
- Figures are PDF, SVG, PNG. All panels use 0–1 axes. Grey lines are iso-F1 contours, not threshold precision–recall curves. Undefined scores are omitted and remain NA in matrices. Very low full_tx scores reflect exact complete exon-boundary matching and are not the same as CDS transcript scores.

## Metric key

| section / metric | Unit and match |
|---|---|
| CDS / cds_ex | Exact CDS segment coordinates |
| CDS / tx | Complete CDS segment set per transcript |
| CDS / locus | One-to-one gene matching via a matching CDS transcript |
| coding_structure / intron | Exact intron between CDS segments |
| coding_structure / ichain | Complete coding intron chain; multi-segment models |
| full_structure / exon | Exact effective exon coordinates, including UTR |
| full_structure / intron | Exact effective exon-derived intron |
| full_structure / ichain | Complete effective exon-derived intron chain |
| full_structure / locus_ic | Gene matching by intron chain |
| full_structure / full_tx | All effective exon boundaries must match, including outer ends |
| base_level / exonic, cds, five_prime_utr, three_prime_utr | Strand-aware base overlap; overlapping isoforms counted once |
| sites / CDS.5prime, CDS.3prime | Strand-aware CDS interval endpoints, without motif validation |
| sites / Coding.Junc.Donor, Coding.Junc.Acceptor | Coding intron boundary coordinates |
| sites / Full.Junc.Donor, Full.Junc.Acceptor | Effective exon-derived intron boundary coordinates |

## Rerun

Run from this directory; supply the actual absolute input list path:

```bash
python3 run_comparisons.py --comparison-list ../comparison_list.tsv --evaluator evaluate_compatible.py --outdir . --structure-policy exon
Rscript plot_comparisons.R comparison.metrics.tsv plots
python3 -m unittest discover -s . -q
```

The batch runner validates input paths, reruns all rows, and stops on any error. It does not silently skip failed comparisons. Input list and original GFF files are not modified. Re-running overwrites generated results in the selected output directory.
