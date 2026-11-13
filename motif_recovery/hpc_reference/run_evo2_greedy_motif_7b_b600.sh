#!/bin/bash
#SBATCH --job-name=evo2_7b_greedy
#SBATCH --partition=b6000-al9_devel
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH -o evo2_7b_greedy_%j.out
#SBATCH -e evo2_7b_greedy_%j.err

set -euo pipefail

INPUT_ROOT="${INPUT_ROOT:-/dicos_ui_home/yenmr/projects/BryoCAD/zero_shot_analysis/zs_recovery_inputs_fullset_modelorganism}"
OUTPUT_ROOT="${OUTPUT_ROOT:-/dicos_ui_home/yenmr/projects/BryoCAD/zero_shot_analysis/zs_recovery_results_fullset_modelorganism_evo2_7b}"
#MODEL_NAME="${MODEL_NAME:-evo2_1b_base}"
MODEL_NAME="${MODEL_NAME:-evo2_7b}"
DEVICE="${DEVICE:-cuda:0}"
BATCH_SIZE="${BATCH_SIZE:-2}"

NAS_ENV_TAR="/ceph/sharedfs/work/TYW_Lab/yenmr/envs/evo2_b600_cuda128_torch271.tar.gz"
HF_CACHE="/ceph/sharedfs/work/TYW_Lab/yenmr/huggingface_cache"
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
SCORE_SCRIPT="${SCRIPT_DIR}/score_evo2_motif_recovery.py"

tasks=(
    "donor_recovery|donor"
    "acceptor_recovery|acceptor"
    "tis_recovery|tis"
    "tts_recovery|tts"
)

[[ -d "${INPUT_ROOT}" ]] || { echo "[FATAL] Missing input root: ${INPUT_ROOT}"; exit 1; }
[[ -f "${SCORE_SCRIPT}" ]] || { echo "[FATAL] Missing scorer: ${SCORE_SCRIPT}"; exit 1; }
[[ -f "${NAS_ENV_TAR}" ]] || { echo "[FATAL] Missing environment: ${NAS_ENV_TAR}"; exit 1; }
mkdir -p "${OUTPUT_ROOT}"

LOCAL_WORKDIR=$(mktemp -d -p /tmp yenmr_evo2_greedy_XXXXXX)
cleanup() {
    status=$?
    echo "[CLEANUP] Removing ${LOCAL_WORKDIR}"
    rm -rf "${LOCAL_WORKDIR}"
    exit "${status}"
}
trap cleanup EXIT

echo "============================================================"
echo "Evo2 forward + reverse-complement greedy motif recovery"
echo "============================================================"
echo "[INFO] Input root:  ${INPUT_ROOT}"
echo "[INFO] Output root: ${OUTPUT_ROOT}"
echo "[INFO] Model:       ${MODEL_NAME}"
echo "[INFO] Batch size:  ${BATCH_SIZE}"

module purge
module load gcc/13.1.0
module load cuda/12.8.0
nvidia-smi --query-gpu=name,driver_version,memory.total,compute_cap --format=csv,noheader

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

CUDA_HOME="$(dirname "$(dirname "$(readlink -f "$(which nvcc)")")")"
CUDNN_ROOT="${CONDA_PREFIX}/lib/python3.11/site-packages/nvidia/cudnn"
export CUDA_HOME
export NVTE_CUDA_INCLUDE_PATH="${CUDA_HOME}/include"
export LD_LIBRARY_PATH="${CUDNN_ROOT}/lib:${CONDA_PREFIX}/lib:${CUDA_HOME}/lib64:${LD_LIBRARY_PATH:-}"

python - <<'PY'
import importlib.metadata
import torch
print("Torch:", torch.__version__)
print("Torch CUDA:", torch.version.cuda)
print("GPU:", torch.cuda.get_device_name(0))
print("Capability:", torch.cuda.get_device_capability(0))
print("Evo2:", importlib.metadata.version("evo2"))
assert torch.cuda.is_available()
assert torch.cuda.get_device_capability(0) == (12, 0)
PY

start_time=$(date +%s)
n_groups=0
n_completed=0
n_skipped=0
n_failed=0

for group_dir in "${INPUT_ROOT}"/*/; do
    [[ -d "${group_dir}" ]] || continue
    group_dir="${group_dir%/}"
    group_name=$(basename "${group_dir}")
    n_groups=$((n_groups + 1))

    echo "============================================================"
    echo "Processing group: ${group_name}"
    echo "============================================================"

    for task_info in "${tasks[@]}"; do
        input_name="${task_info%%|*}"
        evo2_task="${task_info##*|}"
        input_file="${group_dir}/${input_name}.tsv"
        output_file="${OUTPUT_ROOT}/Evo2_${group_name}_${input_name}_greedy.tsv"

        if [[ ! -f "${input_file}" ]]; then
            echo "[SKIP] Missing input: ${input_file}"
            n_skipped=$((n_skipped + 1))
            continue
        fi

        if [[ -s "${output_file}" ]]; then
            echo "[SKIP] Output exists: ${output_file}"
            n_skipped=$((n_skipped + 1))
            continue
        fi

        echo "[RUN] ${group_name} ${evo2_task}"
        set +e
        python "${SCORE_SCRIPT}" \
            --input "${input_file}" \
            --output "${output_file}" \
            --task "${evo2_task}" \
            --model "${MODEL_NAME}" \
            --device "${DEVICE}" \
            --batch-size "${BATCH_SIZE}" \
	    --max-samples 200
        task_status=$?
        set -e

        if [[ "${task_status}" -eq 0 && -s "${output_file}" ]]; then
            echo "[PASS] ${group_name} ${input_name}"
            n_completed=$((n_completed + 1))
        else
            echo "[FAIL] ${group_name} ${input_name}; exit=${task_status}"
            rm -f "${output_file}"
            n_failed=$((n_failed + 1))
        fi
    done
done

elapsed=$(($(date +%s) - start_time))
echo "============================================================"
echo "FINAL SUMMARY"
echo "============================================================"
echo "Groups found:    ${n_groups}"
echo "Tasks completed: ${n_completed}"
echo "Tasks skipped:   ${n_skipped}"
echo "Tasks failed:    ${n_failed}"
echo "Output:          ${OUTPUT_ROOT}"
printf "Runtime:         %02d:%02d:%02d\n" \
    $((elapsed / 3600)) $(((elapsed % 3600) / 60)) $((elapsed % 60))

[[ "${n_failed}" -eq 0 ]] || exit 1
echo "[FINAL] All available benchmarks completed successfully"
