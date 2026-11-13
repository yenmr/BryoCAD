#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
ENV_TARGET=${1:?Usage: bash environments/unified/build_conda_env.sh bryocad /new/genecad/checkout}
GENECAD_SOURCE=${2:?Supply a new directory for the pinned GeneCAD checkout}
[[ "$GENECAD_SOURCE" == /* && ! -e "$GENECAD_SOURCE" ]] || { echo 'Supply a new absolute directory for the GeneCAD checkout'; exit 1; }
# Named environments use conda's configured default envs directory.
# Absolute prefixes remain supported for existing batch build wrappers.
if [[ "$ENV_TARGET" == /* ]]; then
  [[ ! -e "$ENV_TARGET" ]] || { echo 'Environment path already exists'; exit 1; }
  ENV_ARGS=(--prefix "$ENV_TARGET")
else
  [[ "$ENV_TARGET" =~ ^[A-Za-z0-9_][A-Za-z0-9_.-]*$ && "$ENV_TARGET" != base && "$ENV_TARGET" != root ]] || { echo 'Supply a valid new conda environment name'; exit 1; }
  ENV_ARGS=(--name "$ENV_TARGET")
fi
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || { echo 'Build on Linux x86_64 with NVIDIA CUDA 12.8'; exit 1; }
for cmd in conda git nvcc g++; do command -v "$cmd" >/dev/null || { echo "Missing prerequisite: $cmd"; exit 1; }; done
nvcc --version | grep -q 'release 12.8' || { echo 'Load CUDA toolkit 12.8 before building'; exit 1; }
export CUDA_HOME
CUDA_HOME=$(dirname "$(dirname "$(readlink -f "$(command -v nvcc)")")")
export TORCH_CUDA_ARCH_LIST='8.0;8.9;9.0;12.0'
export MAX_JOBS=${MAX_JOBS:-4}
export TOKENIZERS_PARALLELISM=false PYTHONUNBUFFERED=1
conda env create "${ENV_ARGS[@]}" -f "$ROOT/environments/unified/conda.yml"
# Run installation in the new env; do not alter the submitting shell's env.
conda run --no-capture-output "${ENV_ARGS[@]}" bash \
  "$ROOT/environments/unified/install_runtime.sh" "$ROOT" "$GENECAD_SOURCE"
echo "Unified build complete: $ENV_TARGET"
echo "GeneCAD checkout: $GENECAD_SOURCE"
echo "Activate with: conda activate $ENV_TARGET"
echo 'Build records: $CONDA_PREFIX/bryocad-build/ after activation'
