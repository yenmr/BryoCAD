# Installation

Create one conda environment to run BryoCAD motif recovery, SNP scoring, ATAC LoRA prediction, and GeneCAD annotation. Installation and example inference have been tested on an NVIDIA RTX PRO 6000 Blackwell GPU.

## Requirements

- Linux x86_64 with an NVIDIA GPU and a CUDA 12.8-compatible driver.
- CUDA Toolkit 12.8, including `nvcc`.
- Conda, Git, and a CUDA-compatible C++ compiler (`g++`).
- Internet access to download dependencies and model checkpoints.

Run installation on a machine with GPU access. CUDA extensions are compiled during installation, which can take time and substantial memory. On a cluster, request a GPU allocation before running the commands below.

## Create the environment

Clone the repository and enter its directory:

```bash
git clone https://github.com/yenmr/BryoCAD.git
cd BryoCAD
```

Check that the required tools are available:

```bash
conda --version
git --version
nvcc --version
g++ --version
nvidia-smi
```

`nvcc` should report CUDA Toolkit 12.8. Then create the environment and GeneCAD checkout:

```bash
bash environments/unified/build_conda_env.sh \
  bryocad \
  "$HOME/software/genecad-bryocad"
```

The environment name `bryocad` and the GeneCAD checkout directory must not already exist. Conda places the named environment in its configured default environment directory. To choose another environment name, replace `bryocad` in the installation and activation commands. The script installs the required Python packages and CUDA extensions, downloads the compatible GeneCAD source, and runs the example checks. No separate GeneCAD environment is needed.

Use this build script rather than running `conda env create` alone: the YAML file contains only the base environment, while the script completes the installation.

## Activate and use

```bash
conda activate bryocad
export GENECAD_PYTHON="$CONDA_PREFIX/bin/python"
```

Set `GENECAD_PYTHON` whenever you activate this environment to ensure that GeneCAD uses the same Python interpreter.

From the BryoCAD repository root, run all example checks:

```bash
bash scripts/run_gpu_checks.sh "$HOME/software/genecad-bryocad"
```

The checks download the required BryoCAD models and task heads from Hugging Face, then run motif recovery, SNP scoring, ATAC prediction, and GeneCAD annotation. GeneCAD also downloads the models used for protein refinement. Outputs and logs are saved under `results/`; the combined status is written to `results/release_checks/status.tsv`.

To run individual examples:

```bash
# Motif recovery and SNP scoring
bash scripts/run_zero_shot_check.sh

# ATAC adapter loading and TSV prediction
bash scripts/run_atac_check.sh

# Gene annotation
bash scripts/run_genecad_check.sh "$HOME/software/genecad-bryocad"
```

Input formats and recorded outputs are described in the [motif recovery](../examples/motif_recovery/README.md), [SNP](../examples/snp/README.md), [ATAC](../examples/atac/README.md), [GeneCAD](../examples/gene_annotation/README.md), and [example outputs](../examples/expected/README.md) guides.
