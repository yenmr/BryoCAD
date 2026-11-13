#!/bin/bash
#SBATCH --job-name=bryocad_checks
#SBATCH --partition=b6000-al9_devel
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=bryocad_checks_%j.out
#SBATCH --error=bryocad_checks_%j.err
set -euo pipefail
: "${SLURM_JOB_ID:?Submit this script with sbatch}"
SOURCE=${BRYOCAD_SOURCE:-/dicos_ui_home/yenmr/projects/BryoCAD/github/BryoCAD}
ENV_ARCHIVE=${ENV_ARCHIVE:-/ceph/sharedfs/work/TYW_Lab/yenmr/envs/bryocad_unified_cu128_7589721.tar.gz}
# The build records contain the exact patched GeneCAD checkout used to build/test this environment.
GENECAD_RECORDS=${GENECAD_RECORDS:-/ceph/sharedfs/work/TYW_Lab/yenmr/envs/bryocad_unified_records_7589721.tar.gz}
OUTPUT_DIR=${OUTPUT_DIR:-$SOURCE/results/dicos_checks}
[[ -d "$SOURCE/scripts" && -f "$ENV_ARCHIVE" ]] || { echo 'Missing BryoCAD source or environment archive'; exit 2; }
if [[ -n "${GENECAD_SOURCE:-}" ]]; then
  [[ -f "$GENECAD_SOURCE/predict.sh" ]] || { echo 'GENECAD_SOURCE must be a patched checkout'; exit 2; }
else
  [[ -f "$GENECAD_RECORDS" ]] || { echo 'Set GENECAD_RECORDS to the previous build records archive, or GENECAD_SOURCE to the patched checkout'; exit 2; }
fi
mkdir -p "$OUTPUT_DIR"
WORK=$(mktemp -d "/tmp/${USER}_bryocad_checks_${SLURM_JOB_ID}_XXXXXX")
exec > >(tee "$WORK/submission.log") 2>&1
cleanup() {
  status=$?
  trap - EXIT
  set +e
  printf '%s\n' "$status" > "$WORK/job_exit_code.txt"
  records=(submission.log job_exit_code.txt)
  for item in provenance.txt native-kernel-test.log BryoCAD/results; do
    [[ ! -e "$WORK/$item" ]] || records+=("$item")
  done
  archive="$OUTPUT_DIR/bryocad_checks_${SLURM_JOB_ID}.tar.gz"
  if tar -czf "$archive.partial" -C "$WORK" "${records[@]}" && mv "$archive.partial" "$archive"; then
    echo "Results: $archive (job exit $status)"
    rm -rf "$WORK"
  else
    echo "Could not save all results; scratch retained: $WORK"
    [[ "$status" != 0 ]] || status=74
  fi
  exit "$status"
}
trap cleanup EXIT
module purge
module load gcc/13.1.0
module load cuda/12.8.0
mkdir -p "$WORK/env"
tar -xzf "$ENV_ARCHIVE" -C "$WORK/env"
export PATH="$WORK/env/bin:$PATH"
unset PYTHONHOME PYTHONPATH
export PYTHONNOUSERSITE=1
export CONDA_PREFIX="$WORK/env"
"$WORK/env/bin/conda-unpack"
rsync -a --exclude='.git/' --exclude='.venv*/' --exclude='results/' \
  --exclude='checkpoints/' --exclude='hf_cache/' "$SOURCE/" "$WORK/BryoCAD/"
if [[ -n "${GENECAD_SOURCE:-}" ]]; then
  rsync -a --exclude='.git/' --exclude='.venv*/' --exclude='results/' "$GENECAD_SOURCE/" "$WORK/genecad/"
else
  tar -xzf "$GENECAD_RECORDS" -C "$WORK" genecad
fi
export PYTHON="$WORK/env/bin/python" BRYOCAD_PYTHON="$WORK/env/bin/python" GENECAD_PYTHON="$WORK/env/bin/python"
export HF_HOME="$WORK/huggingface" HF_HUB_CACHE="$WORK/huggingface/hub"
export HUGGINGFACE_HUB_CACHE="$HF_HUB_CACHE"
unset HF_HUB_OFFLINE TRANSFORMERS_OFFLINE HF_DATASETS_OFFLINE TRANSFORMERS_CACHE
export MOTIF_DTYPE=bfloat16 GENECAD_DTYPE=bfloat16
export OMP_NUM_THREADS=2
printf 'Job: %s\nSource: %s\nEnvironment: %s\nGeneCAD records: %s\nGeneCAD override: %s\n' \
  "$SLURM_JOB_ID" "$SOURCE" "$ENV_ARCHIVE" "$GENECAD_RECORDS" "${GENECAD_SOURCE:-none}" > "$WORK/provenance.txt"
cd "$WORK/BryoCAD"
# Continue to the model checks even if the standalone GPU kernel probe fails.
failed=0
"$PYTHON" environments/unified/test_native_kernels.py > "$WORK/native-kernel-test.log" 2>&1 || failed=1
cat "$WORK/native-kernel-test.log"
bash scripts/run_gpu_checks.sh "$WORK/genecad" || failed=1
exit "$failed"
