#!/bin/bash
#SBATCH --job-name=bryocad_env_build
#SBATCH --partition=b6000-al9_devel
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=192G
#SBATCH --time=06:00:00
#SBATCH --output=bryocad_env_build_%j.out
#SBATCH --error=bryocad_env_build_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch}"
SOURCE="${BRYOCAD_SOURCE:-/dicos_ui_home/yenmr/projects/BryoCAD/github/BryoCAD}"
OUTPUT_DIR="${OUTPUT_DIR:-/ceph/sharedfs/work/TYW_Lab/yenmr/envs}"
BUILD_ROOT=$(mktemp -d "/tmp/${USER}_bryocad_unified_${SLURM_JOB_ID}_XXXXXX")
ENV_PREFIX="$BUILD_ROOT/env"
mkdir -p "$OUTPUT_DIR"
exec > >(tee "$BUILD_ROOT/submission.log") 2>&1
cleanup() {
  status=$?
  trap - EXIT
  set +e
  copied=0
  records=(submission.log)
  for relative in env/bryocad-build BryoCAD/results genecad; do
    [[ ! -e "$BUILD_ROOT/$relative" ]] || records+=("$relative")
  done
  tar -czf "$OUTPUT_DIR/bryocad_unified_records_${SLURM_JOB_ID}.tar.gz" \
    --exclude='genecad/.git' --exclude='*.bin' --exclude='*.ckpt' \
    -C "$BUILD_ROOT" "${records[@]}" || copied=$?
  if [[ "$copied" == 0 ]]; then rm -rf "$BUILD_ROOT";
  else echo "[WARNING] Some records could not be copied; scratch retained: $BUILD_ROOT"; [[ "$status" != 0 ]] || status=74; fi
  exit "$status"
}
trap cleanup EXIT
module purge
module load gcc/13.1.0
module load cuda/12.8.0
module load anaconda3/2024.10-1
rsync -a --exclude='.git/' --exclude='.venv*/' --exclude='results/' \
  --exclude='checkpoints/' --exclude='hf_cache/' "$SOURCE/" "$BUILD_ROOT/BryoCAD/"
export HF_HOME="$BUILD_ROOT/huggingface" HF_HUB_CACHE="$BUILD_ROOT/huggingface/hub"
export HUGGINGFACE_HUB_CACHE="$HF_HUB_CACHE"
unset HF_HUB_OFFLINE TRANSFORMERS_OFFLINE HF_DATASETS_OFFLINE TRANSFORMERS_CACHE
export MAX_JOBS=4
bash "$BUILD_ROOT/BryoCAD/environments/unified/build_conda_env.sh" "$ENV_PREFIX" "$BUILD_ROOT/genecad"
# Separate archive name: never overwrite the previously working environment.
"$ENV_PREFIX/bin/conda-pack" -p "$ENV_PREFIX" -o "$OUTPUT_DIR/bryocad_unified_cu128_${SLURM_JOB_ID}.tar.gz"
echo "PASS: $OUTPUT_DIR/bryocad_unified_cu128_${SLURM_JOB_ID}.tar.gz"
