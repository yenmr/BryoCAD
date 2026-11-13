# Examples

Run these examples to learn the input formats and outputs for each task. Activate the [conda environment](../environments/README.md) and run all commands from the repository root on a GPU machine:

```bash
conda activate bryocad
export GENECAD_PYTHON="$CONDA_PREFIX/bin/python"
```

The scripts download the required models from Hugging Face and reuse cached files.

## Motif recovery

Predict donor, acceptor, translation initiation (TIS), and translation termination (TTS) motifs. Each task includes four 8,192-bp sequences.

```bash
bash examples/motif_recovery/run.sh
```

Output: `results/toy_motif_recovery/{task}.tsv`, containing predicted motifs and exact-recovery labels.

See [motif input format and output details](motif_recovery/README.md).

## SNP scoring

Score 20 SNPs using the supplied 200-kb reference sequence.

```bash
bash scripts/run_zero_shot_check.sh
```

This command also runs the motif examples. SNP output: `results/zero_shot_check/snps.vcf`, containing masked-site log-likelihood ratios.

See [SNP input format and score interpretation](snp/README.md).

## Chromatin accessibility

Predict accessibility for ten 512-bp sequences supplied in TSV format.

```bash
bash scripts/run_atac_check.sh
```

Output: `results/atac_predict_check/predictions.tsv`, containing positive-class probabilities and predicted labels alongside the input metadata.

See [ATAC input format and custom inputs](atac/README.md).

## Gene annotation

Annotate a 200-kb genomic fragment using the BryoCAD backbone and GeneCAD head.

```bash
bash scripts/run_genecad_check.sh "$HOME/software/genecad-bryocad"
```

Replace the GeneCAD source path if you installed it elsewhere. Output: `results/genecad_check/predictions/toy_mpolymorpha_GeneCAD_final.gff`.

See [GeneCAD inputs and output coordinates](gene_annotation/README.md).

## Run all examples

```bash
bash scripts/run_gpu_checks.sh "$HOME/software/genecad-bryocad"
```

This runs all four tasks and checks their output structure. The status is saved to `results/release_checks/status.tsv`; task-specific `run.log` files provide error details if a command fails. A successful run demonstrates that inference works with the supplied inputs; these small examples are not performance benchmarks.
