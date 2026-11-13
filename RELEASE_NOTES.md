# BryoCAD initial public release — v0.1.0

This release provides BryoCAD model links, analysis scripts, a unified conda environment recipe, and runnable motif recovery, SNP scoring, ATAC LoRA and GeneCAD examples.

Validated on NVIDIA RTX PRO 6000 Blackwell with PyTorch 2.7.1+cu128, Transformers 4.49.0 and Python 3.12.16. A fresh environment build and subsequent execution from its conda-pack archive passed. The 200 kb GeneCAD example completed protein refinement and generated its final GFF. Recorded example outputs are included under `examples/expected`.

A100, L40S and H100 are build targets but have not been tested with this release. Examples validate execution rather than reproduce all manuscript analyses. Known upstream package-metadata differences are documented in the environment guide.
