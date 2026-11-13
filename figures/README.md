# Plotting analysis results

These workflows read existing result tables and generate plots. They do not run model inference. Run commands from the repository root; generated files are written under `results/`.

## Requirements

Python 3 and R are required. The workflows below use the R packages `ggplot2`, `readr`, `dplyr`, and `tidyr`:

```r
install.packages(c("ggplot2", "readr", "dplyr", "tidyr"))
```

## Model-organism motif recovery

Input: `source_data/model_organisms/by_task_and_species_merged.csv` in this directory.

```bash
bash figures/run_selected.sh model-organisms
```

The workflow summarizes motif recovery and plots comparisons between models. Outputs are written to `results/model_organisms/`, including `exact_recovery_long.tsv` and `exact_recovery_overview.pdf`.

## Transfer across plant lineages

Input: `source_data/lineage_transfer/combined_exact_recovery.csv`.

```bash
bash figures/run_selected.sh lineage-transfer
```

Outputs are written to `results/lineage_transfer/`. The workflow generates separate 22-species, 18-species, and 15-species views. The 15-species view excludes the seven model organisms; the 18-species view excludes the four flowering-plant model organisms. Choose the view appropriate to your comparison rather than treating these as interchangeable datasets.

## Further analysis of TTS motif recovery

These workflows examine stop-codon predictions using existing TTS results. See the [TTS analysis guide](../motif_recovery/tts/README.md) for analysis details.

For bryophytes, input tables are in `source_data/bryophyte_recovery/`:

```bash
bash figures/run_selected.sh bryophyte-tts
```

Outputs: `results/tts_bryophytes/`.

For the angiosperm panel, inputs are in `source_data/angiosperm_tts/`:

```bash
bash figures/run_selected.sh angiosperm-tts-panel
```

Outputs: `results/tts_angiosperms/`.

## Gene annotation

Inputs:

- `source_data/gene_annotation/comparison.metrics.tsv`
- `source_data/pangenome_candidates/pangenome.metrics.tsv`

```bash
bash figures/run_selected.sh gene-annotation
```

Outputs are written to `results/gene_annotation_coding/`. Comparisons use CDS segments and CDS-chain transcript metrics, excluding UTRs.

## Figure mapping and custom data

The [figure manifest](figure_manifest.tsv) maps manuscript panels to scripts and source data. Some entries are marked as candidate or unavailable; those entries do not establish exact reproduction of the final manuscript panel. Statistical methods are described in [STATISTICS.md](STATISTICS.md).

The wrapper uses the input locations listed above. To plot your own results, provide tables with the same columns and identifier conventions, and update the input paths in the corresponding script or wrapper. Keep a separate copy of the supplied source tables. Plotting these summaries does not rerun the underlying inference or annotation workflows.
