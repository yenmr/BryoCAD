# Figures and source data

`figure_manifest.tsv` maps the current six main and nine supporting figures to available data and analysis scripts. Entries marked as candidate or unavailable are not claimed to recreate the final panel exactly.

Run individual prepared workflows from the repository root:

```bash
bash figures/run_selected.sh model-organisms
bash figures/run_selected.sh lineage-transfer
bash figures/run_selected.sh bryophyte-tts
bash figures/run_selected.sh angiosperm-tts-panel
bash figures/run_selected.sh gene-annotation
```

The gene-annotation entry point selects CDS segments and CDS-chain transcript metrics (UTRs excluded) for Figures 5A and 5D. It does not use the historical full-transcript/UTR-sensitive plots. Current Figure 6 uses 15 additional species; S8 contains all 22. The 18-species subset is retained only as a historical analysis output.

Source-data directories labeled `candidates` retain analysis snapshots, not automatically final manuscript values. Figure 3/S1/S2 group comparisons use one-sided Welch tests with BH correction; see `STATISTICS.md`. Protein-generation changes and code preparation do not automatically regenerate historical orthogroups or benchmark scores.
