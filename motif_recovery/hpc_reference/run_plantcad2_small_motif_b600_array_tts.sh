#!/bin/bash
#SBATCH --job-name=plantcad2_small_tts
#SBATCH --partition=b6000-al9
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=5-00:00:00
#SBATCH --array=1-132%4
# Override at submission: sbatch --array="1-${n_species}%4" this_script.sh
#SBATCH -o logs/plantcad2_small_tts_%A_%a.out
#SBATCH -e logs/plantcad2_small_tts_%A_%a.err

set -euo pipefail

INPUT_ROOT="${INPUT_ROOT:-/dicos_ui_home/yenmr/projects/BryoCAD/zero_shot_analysis/zs_recovery_inputs_tts_angiosperms}"
OUTPUT_ROOT="${OUTPUT_ROOT:-/dicos_ui_home/yenmr/projects/BryoCAD/zero_shot_analysis/zs_recovery_results_tts_angiosperms_plantcad2_small}"
MODEL_NAME="${MODEL_NAME:-kuleshov-group/PlantCAD2-Small-l24-d0768}"
NAS_ENV_TAR="${NAS_ENV_TAR:-/ceph/sharedfs/work/TYW_Lab/yenmr/envs/genecad_b600_cuda128_torch271.tar.gz}"
HF_CACHE="${HF_CACHE:-/ceph/sharedfs/work/TYW_Lab/yenmr/huggingface_cache}"
DEVICE="${DEVICE:-cuda:0}"
DTYPE="${DTYPE:-auto}"
BATCH_SIZE="${BATCH_SIZE:-16}"
NUM_WORKERS="${NUM_WORKERS:-2}"

BASE_DIR="${BASE_DIR:-/dicos_ui_home/yenmr/projects/BryoCAD/zero_shot_analysis}"
SCORE_SCRIPT="${SCORE_SCRIPT:-${BASE_DIR}/score_plantcad2_motif_recovery.py}"
SPECIES_LIST="${SPECIES_LIST:-${BASE_DIR}/angiospecies_species_list.txt}"

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || {
    echo "[FATAL] Submit this script as a Slurm array job."
    exit 1
}
[[ -d "${INPUT_ROOT}" ]] || { echo "[FATAL] Missing input root: ${INPUT_ROOT}"; exit 1; }
[[ -f "${SPECIES_LIST}" ]] || { echo "[FATAL] Missing species list: ${SPECIES_LIST}"; exit 1; }
[[ -f "${SCORE_SCRIPT}" ]] || { echo "[FATAL] Missing scorer: ${SCORE_SCRIPT}"; exit 1; }
[[ -f "${NAS_ENV_TAR}" ]] || { echo "[FATAL] Missing environment: ${NAS_ENV_TAR}"; exit 1; }

species=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SPECIES_LIST}")
species="${species%$'\r'}"
[[ -n "${species}" ]] || {
    echo "[FATAL] No species at line ${SLURM_ARRAY_TASK_ID} of ${SPECIES_LIST}"
    exit 1
}
[[ "${species}" != *[[:space:]]* && "${species}" != */* && "${species}" != "." && "${species}" != ".." ]] || {
    echo "[FATAL] Invalid species name in list: ${species}"
    exit 1
}

SPECIES_INPUT_DIR="${INPUT_ROOT}/${species}"
SPECIES_OUTPUT_DIR="${OUTPUT_ROOT}/${species}"
[[ -d "${SPECIES_INPUT_DIR}" ]] || {
    echo "[FATAL] Missing species input directory: ${SPECIES_INPUT_DIR}"
    exit 1
}
mkdir -p "${SPECIES_OUTPUT_DIR}"

tasks=(
    "tts_recovery|tts"
)

LOCAL_WORKDIR=$(mktemp -d -p /tmp yenmr_plantcad2_small_tts_XXXXXX)
cleanup() {
    status=$?
    echo "[CLEANUP] Removing ${LOCAL_WORKDIR}"
    rm -rf "${LOCAL_WORKDIR}"
    exit "${status}"
}
trap cleanup EXIT

echo "============================================================"
echo "PlantCAD2-Small zero-shot TTS recovery"
echo "Array job: ${SLURM_ARRAY_JOB_ID:-unknown}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Species: ${species}"
echo "Model: ${MODEL_NAME}"
echo "Batch size: ${BATCH_SIZE}"
echo "Input: ${SPECIES_INPUT_DIR}"
echo "Output: ${SPECIES_OUTPUT_DIR}"
echo "============================================================"

module purge
module load gcc/13.1.0
module load cuda/12.8.0

tar -xzf "${NAS_ENV_TAR}" -C "${LOCAL_WORKDIR}"
set +u
source "${LOCAL_WORKDIR}/bin/activate"
if command -v conda-unpack >/dev/null 2>&1; then
    conda-unpack
fi
set -u

export PYTHONUNBUFFERED=1
export TOKENIZERS_PARALLELISM=false
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export TORCH_CUDA_ARCH_LIST="12.0"
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"
export HF_HOME="${HF_CACHE}"
export HUGGINGFACE_HUB_CACHE="${HF_CACHE}/hub"
mkdir -p "${HF_HOME}"
unset TRANSFORMERS_CACHE
unset HF_DATASETS_OFFLINE
unset TRANSFORMERS_OFFLINE
unset HF_HUB_OFFLINE

CUDA_HOME="$(dirname "$(dirname "$(readlink -f "$(command -v nvcc)")")")"
export CUDA_HOME
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${CUDA_HOME}/lib64:${LD_LIBRARY_PATH:-}"

start_time=$(date +%s)
n_completed=0
n_skipped=0
n_failed=0

for task_info in "${tasks[@]}"; do
    input_name="${task_info%%|*}"
    motif_task="${task_info##*|}"
    input_file="${SPECIES_INPUT_DIR}/${input_name}.tsv"
    output_file="${SPECIES_OUTPUT_DIR}/PlantCAD2_Small_${species}_${input_name}.tsv"

    if [[ ! -f "${input_file}" ]]; then
        echo "[FAIL] Missing input: ${input_file}"
        n_failed=$((n_failed + 1))
        continue
    fi
    if [[ -s "${output_file}" ]]; then
        echo "[SKIP] Output exists: ${output_file}"
        n_skipped=$((n_skipped + 1))
        continue
    fi

    temp_output="${output_file}.partial.${SLURM_JOB_ID:-$$}.${SLURM_ARRAY_TASK_ID}"
    echo "[RUN] ${species} ${motif_task}"
    set +e
    python "${SCORE_SCRIPT}" \
        --input "${input_file}" \
        --output "${temp_output}" \
        --task "${motif_task}" \
        --model "${MODEL_NAME}" \
        --device "${DEVICE}" \
        --dtype "${DTYPE}" \
        --batch-size "${BATCH_SIZE}" \
        --num-workers "${NUM_WORKERS}"
    task_status=$?
    set -e

    if [[ "${task_status}" -eq 0 && -s "${temp_output}" ]]; then
        mv -f "${temp_output}" "${output_file}"
        echo "[PASS] ${species} ${input_name}"
        n_completed=$((n_completed + 1))
    else
        echo "[FAIL] ${species} ${input_name}; exit=${task_status}"
        rm -f "${temp_output}"
        n_failed=$((n_failed + 1))
    fi
done

elapsed=$(($(date +%s) - start_time))
echo "============================================================"
echo "Species:   ${species}"
echo "Completed: ${n_completed}"
echo "Skipped:   ${n_skipped}"
echo "Failed:    ${n_failed}"
printf "Runtime:   %02d:%02d:%02d\n" \
    $((elapsed / 3600)) $(((elapsed % 3600) / 60)) $((elapsed % 60))
echo "============================================================"

[[ "${n_failed}" -eq 0 ]] || exit 1
