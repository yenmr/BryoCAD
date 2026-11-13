# Chromatin accessibility prediction

`input_template.tsv` contains 10 author-supplied validation examples: three peaks (label 1) and seven negatives (label 0), with five windows each from chr7 and chr8. All sequences are 512 bp. This small subset tests inference and output formatting; it is not a performance benchmark.

Required column: `sequence` (exactly 512 A/C/G/T/N bases). `label` is optional and must be 0/1 if supplied. Other metadata columns are retained. No Parquet conversion or manual tokenization is required. Prediction does not use labels as model input and does not require every positive example to receive a positive prediction.

```bash
bash scripts/run_atac_predict_check.sh
# Or supply your own TSV:
bash scripts/run_atac_predict_check.sh /path/to/Validation_subset.tsv
```

This downloads the pinned base model and adapter, then calls `chromatin_accessibility/lora_fine_tune.py predict`. Output `results/atac_predict_check/predictions.tsv` retains input columns and adds `probability_positive` and `predicted_label` (argmax; ties choose class 0). Logs and environment versions are saved alongside it. `run_atac_check.sh` now runs both the direct adapter-restoration check and this original prediction CLI; it accepts the same optional TSV path.

Use the [conda environment](../../environments/README.md) to install the required dependencies. No model training is needed to run this example.
