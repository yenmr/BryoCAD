# Stop-codon prediction bias

Separate reference codon usage, predicted codon frequency and exact recovery. A preference for predicting TAA at TAG/TGA sites is not evidence that reference TAA is abundant. GC association alone does not establish a causal mechanism.

```bash
Rscript motif_recovery/tts/analyze_tts.R   figures/source_data/bryophyte_recovery results/tts_bryophytes BryoCAD-Small
```

Inputs: `tts_codon_confusion.tsv` and `species_metadata.tsv`. The script handles canonical reference codons, predicted noncanonical outcomes and composition-preserving null comparisons. The 132-angiosperm extension has summary source tables bundled, but its processed confusion/metadata inputs still need to be located to rerun the entire analysis.

`assemble_tts_figure.R` is retained as the original panel-assembly script; final panel choices and path portability remain to be reviewed.
