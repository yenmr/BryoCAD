#!/usr/bin/env bash
# Called by build_conda_env.sh through conda run.
set -euo pipefail
ROOT=${1:?}; GENECAD_SOURCE=${2:?}
: "${CONDA_PREFIX:?Must run inside the newly created conda environment}"
RECORDS="$CONDA_PREFIX/bryocad-build"
mkdir -p "$RECORDS"
exec > >(tee "$RECORDS/build.log") 2>&1
BUILD_TMP=$(mktemp -d "${TMPDIR:-/tmp}/bryocad-native.XXXXXX")
cleanup() {
  status=$?
  trap - EXIT
  set +e
  python "$ROOT/environments/unified/collect_environment.py" --output "$RECORDS/environment"
  printf '%s\n' "$status" > "$RECORDS/build_exit_code.txt"
  rm -rf "$BUILD_TMP"
  exit "$status"
}
trap cleanup EXIT
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$CUDA_HOME/lib64:${LD_LIBRARY_PATH:-}"
python -m pip install 'torch==2.7.1' 'torchvision==0.22.1' 'torchaudio==2.7.1' --index-url https://download.pytorch.org/whl/cu128
python -m pip install 'pip==26.2.1' 'setuptools==84.0.0' 'wheel==0.48.0' 'packaging==26.3' 'ninja==1.13.0' hatchling
python -m pip install -c "$ROOT/environments/unified/constraints.txt" -r "$ROOT/environments/unified/requirements.txt"
# Follow the successful dicos recipe: avoid pulling a second CUDA/NCCL stack.
python -m pip install --no-deps 'xgboost==3.4.0'
for component in causal-conv1d mamba; do
  if [[ "$component" == causal-conv1d ]]; then
    url=https://github.com/Dao-AILab/causal-conv1d.git; tag=v1.5.0.post8
  else
    url=https://github.com/state-spaces/mamba.git; tag=v2.2.4
  fi
  git clone --branch "$tag" --depth 1 --recurse-submodules "$url" "$BUILD_TMP/$component"
  git -C "$BUILD_TMP/$component" rev-parse HEAD > "$RECORDS/$component.commit.txt"
done
python "$ROOT/environments/unified/patch_cuda_architectures.py" \
  "$BUILD_TMP/causal-conv1d/setup.py" "$BUILD_TMP/mamba/setup.py" \
  --report "$RECORDS/native_architecture_patch.json"
for component in causal-conv1d mamba; do
  git -C "$BUILD_TMP/$component" diff -- setup.py > "$RECORDS/$component.patch"
done
export CAUSAL_CONV1D_FORCE_BUILD=TRUE MAMBA_FORCE_BUILD=TRUE
python -m pip install --no-build-isolation --no-cache-dir --no-deps "$BUILD_TMP/causal-conv1d"
python -m pip install --no-build-isolation --no-cache-dir --no-deps "$BUILD_TMP/mamba"
# Reuse the FlashAttention version from the working B6000 environment.
python -m pip install --no-build-isolation --no-cache-dir --no-deps 'flash-attn==2.7.4.post1'

git clone https://github.com/plantcad/genecad.git "$GENECAD_SOURCE"
git -C "$GENECAD_SOURCE" checkout 1a9baf84d61152c20a0cd910694316cd88d1d374
git -C "$GENECAD_SOURCE" apply --check "$ROOT/gene_annotation/upstream_modifications/genecad-v0.4.0-bryocad.patch"
git -C "$GENECAD_SOURCE" apply "$ROOT/gene_annotation/upstream_modifications/genecad-v0.4.0-bryocad.patch"
git -C "$GENECAD_SOURCE" rev-parse HEAD > "$RECORDS/genecad.commit.txt"
# Non-editable package installation keeps a packed environment relocatable.
# Dependency resolution is deliberately controlled by the observed CLI recipe.
python -m pip install --no-deps --no-build-isolation "$GENECAD_SOURCE"
python - <<'PY'
import importlib.metadata as md
assert md.version('torch') == '2.7.1+cu128'
assert md.version('transformers') == '4.49.0'
assert md.version('peft') == '0.14.0'
try:md.version('nvidia-nccl-cu13')
except md.PackageNotFoundError:pass
else:raise SystemExit('Unexpected nvidia-nccl-cu13 package')
PY
cp "$ROOT/environments/unified/requirements.txt" "$RECORDS/requirements.txt"
export PYTHON="$CONDA_PREFIX/bin/python" GENECAD_PYTHON="$CONDA_PREFIX/bin/python"
if [[ ${RUN_GPU_CHECKS:-1} == 1 ]]; then
  python "$ROOT/environments/unified/test_native_kernels.py" | tee "$RECORDS/native-kernel-test.log"
  export MOTIF_DTYPE=bfloat16 GENECAD_DTYPE=bfloat16
  bash "$ROOT/scripts/run_gpu_checks.sh" "$GENECAD_SOURCE"
  cp "$ROOT/results/release_checks/status.tsv" "$RECORDS/model-checks.tsv"
else
  echo 'NOT RUN: GPU validation explicitly skipped' > "$RECORDS/model-checks.tsv"
fi
