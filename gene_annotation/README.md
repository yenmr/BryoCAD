# Gene annotation

GeneCAD annotation heads are trained on a frozen BryoCAD-Small or PlantCAD2-Small backbone. The 2bryo set contains A. agrestis and M. polymorpha; 4bryo adds C. purpureus and S. magellanicum. Species configurations are in `configs/`.

## Software and models

Use the pinned upstream GeneCAD v0.4.0 implementation plus the supplied [modification patch](UPSTREAM_VERSION.md). Original HPC jobs are retained in `hpc_reference/`. The effective settings come from the job arguments as well as internal defaults.

The published [4bryo head](https://huggingface.co/bryocad/GeneCAD-BryoCAD-Small-4bryo) pairs with [BryoCAD-Small](https://huggingface.co/bryocad/BryoCAD-Small).

```bash
bash /path/to/patched-genecad/predict.sh \
  --input-fasta genome.fa --output-dir results/annotation --species-id sample \
  --mode plant --top-n-contigs all --cpu-workers 4 --batch-size 4 --gpus 0 \
  --model-checkpoint /path/to/model.ckpt --base-model bryocad/BryoCAD-Small
```

## Evaluation and examples

- [Figure 5A/5D metrics and selected-contig scope](EVALUATION_SCOPE.md)
- [Figure 5B probability/annotation plotting](FIGURE5B.md)
- `clean_genecad_gff.py`: representative-isoform annotation preparation.
- `comparison_snapshot/evaluate_compatible.py`: evaluator stored with the model-organism comparison outputs.
- `evaluate_compatible.py`: separately supplied pangenome evaluator; preserve evaluator/output version pairing.

M. polymorpha appears in head training; its evaluation should not be described as unseen-species testing. The 123-assembly Figure 5D analysis uses selected chromosomes/contigs, targeting roughly 50 Mb per assembly. Coding-transcript metrics exclude UTRs in both Figures 5A and 5D.

## Download predicted annotations

[GeneCAD annotations for five bryophyte genomes](https://predict.southerngenomics.org/bryocad/annotation/) are available as GFF files. These are predicted gene annotations, separate from the genome FASTA assemblies and model checkpoints.

## Use the unified conda environment

The [unified conda recipe](../environments/README.md) includes GeneCAD and BryoCAD dependencies. Activate it and set `GENECAD_PYTHON="$CONDA_PREFIX/bin/python"` before running the test. For existing installations, the runner can still use an explicit interpreter or the checkout's `.venv`; neither virtualenv nor a separate environment is required. Preflight checks Python requirements, imports and CUDA before model downloads. The smoke test makes one prediction attempt, not 20 retries.
