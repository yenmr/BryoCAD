# Figure 5B: probabilities and annotation tracks

The author supplied `plot_gene_annotation_region.py`; the filename refers to a historical figure number. It reads GeneCAD Zarr prediction stores, predicted GFF3 and reference GFF3. Dependencies: numpy, xarray, zarr and matplotlib in the working GeneCAD environment.

Positive-strand panel (author-supplied command; coordinates are 1-based inclusive):

```bash
python gene_annotation/plot_gene_annotation_region.py \
  --predictions /path/to/prediction_result/MpTak_chr1/chr1/predictions_chr1/ \
  --prediction-gff /path/to/MpTak_chr1/MpTak_chr1_GeneCAD_final.gff \
  --annotation-gff /path/to/MpTak_genecad.gff3.gz \
  --chromosome chr1 --strand + --start 9812186 --end 9821880 \
  --output results/MpTak_chr1_positive_Mp1g09490.svg \
  --save-npz results/MpTak_chr1_positive_Mp1g09490.npz
```

The negative-strand panel shown in the PDF uses chr1:8697123–8705954, strand `-`. These coordinates were read from the figure; the exact original negative-strand invocation has not been supplied. Transcription orientation is the script default.

The plotting code is now available. To regenerate without rerunning genome prediction, retain the corresponding prediction Zarr data and GFFs. `--save-npz` exports plot data, but the current script does not accept NPZ as an input replacement for Zarr. GPU prediction and this data-dependent plot have not been rerun in this release environment.
