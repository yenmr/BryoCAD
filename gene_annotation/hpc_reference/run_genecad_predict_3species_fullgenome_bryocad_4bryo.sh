#!/bin/bash
#SBATCH --job-name=genecad_3species_full
#SBATCH --partition=b6000-al9
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=5-00:00:00
#SBATCH --output=genecad_3species_full_bryocad_4bryo_%j.out
#SBATCH --error=genecad_3species_full_bryocad_4bryo_%j.err

set -euo pipefail

# =============================================================================
# Software and model
# =============================================================================

ENV_TAR="/ceph/sharedfs/work/TYW_Lab/yenmr/envs/genecad_b600_cuda128_torch271.tar.gz"
GENECAD_SOURCE="/dicos_ui_home/yenmr/projects/BryoCAD/annotation/genecad"
HF_CACHE="/ceph/sharedfs/work/TYW_Lab/yenmr/huggingface_cache"

#BASE_MODEL="kuleshov-group/PlantCAD2-Small-l24-d0768"
BASE_MODEL="/ceph/sharedfs/work/TYW_Lab/yenmr/projects/BryoCAD/warehouse/models/bryocad_small_final"

#HEAD_CHECKPOINT="/dicos_ui_home/yenmr/projects/BryoCAD/annotation/genecad_output/bryophyte_training/plantcad2_small_d768_4bryo_v1/checkpoints/best.ckpt"
HEAD_CHECKPOINT="/dicos_ui_home/yenmr/projects/BryoCAD/annotation/genecad_output/bryophyte_training/bryocad_small_d768_4bryo_v1/checkpoints/best.ckpt"

INPUT_ROOT="/dicos_ui_home/yenmr/projects/BryoCAD/warehouse/library/genomes3"

OUTPUT_ROOT="/dicos_ui_home/yenmr/projects/BryoCAD/annotation/genecad_output/prediction/full_genome/bryocad_small_d768_4bryo_3species"

# Batch size 67 was previously confirmed to work on this B600.
PREDICT_BATCH_SIZE="${PREDICT_BATCH_SIZE:-67}"

# Set to yes only when intentionally continuing an interrupted prediction.
ALLOW_RESUME="${ALLOW_RESUME:-no}"

# =============================================================================
# Input definitions
#
# Format:
# species_id|absolute FASTA path|absolute reference GFF path
#
# The array order is the execution order.
# =============================================================================

TESTS=(
    "Ath|${INPUT_ROOT}/Araport11.fa.gz|${INPUT_ROOT}/Araport11_genecad.gff3.gz"
    "Mpo|${INPUT_ROOT}/MpTak_v6.fa.gz|${INPUT_ROOT}/MpTak_genecad.gff3.gz"
    "Ppa|${INPUT_ROOT}/Ppatens_318_v3.softmasked.fa.gz|${INPUT_ROOT}/Ppatens_genecad.gff3.gz"
)

# =============================================================================
# Validate common inputs
# =============================================================================

for path in \
    "${ENV_TAR}" \
    "${GENECAD_SOURCE}/pyproject.toml" \
    "${GENECAD_SOURCE}/predict.sh" \
    "${HEAD_CHECKPOINT}"
do
    [[ -e "${path}" ]] || {
        echo "[FATAL] Missing: ${path}"
        exit 1
    }
done

for test_info in "${TESTS[@]}"; do
    IFS='|' read -r species_id fasta_path reference_path <<< "${test_info}"

    [[ -s "${fasta_path}" ]] || {
        echo "[FATAL] Missing FASTA for ${species_id}: ${fasta_path}"
        exit 1
    }

    [[ -s "${reference_path}" ]] || {
        echo "[FATAL] Missing reference GFF for ${species_id}: ${reference_path}"
        exit 1
    }
done

mkdir -p "${OUTPUT_ROOT}" "${HF_CACHE}"

# =============================================================================
# Local runtime
# =============================================================================

LOCAL_WORKDIR=$(mktemp -d -p /tmp "${USER}_genecad_3species_XXXXXX")
LOCAL_ENV="${LOCAL_WORKDIR}/env"
LOCAL_REPO="${LOCAL_WORKDIR}/genecad"

cleanup() {
    status=$?
    echo "[CLEANUP] Removing ${LOCAL_WORKDIR}"
    rm -rf "${LOCAL_WORKDIR}"
    exit "${status}"
}
trap cleanup EXIT

module purge
module load gcc/13.1.0
module load cuda/12.8.0
module load anaconda3/2024.10-1

mkdir -p "${LOCAL_ENV}"

tar -xzf "${ENV_TAR}" -C "${LOCAL_ENV}"
cp -a "${GENECAD_SOURCE}" "${LOCAL_REPO}"

set +u
source "${LOCAL_ENV}/bin/activate"
if command -v conda-unpack >/dev/null 2>&1; then
    conda-unpack
fi
set -u

export GENECAD_PYTHON="${LOCAL_ENV}/bin/python"
export PYTHONUNBUFFERED=1
export TOKENIZERS_PARALLELISM=false
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export TORCH_CUDA_ARCH_LIST="12.0"
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"
export WANDB_MODE=disabled

export HF_HOME="${HF_CACHE}"
export HUGGINGFACE_HUB_CACHE="${HF_CACHE}/hub"

unset TRANSFORMERS_CACHE
unset HF_DATASETS_OFFLINE
unset TRANSFORMERS_OFFLINE
unset HF_HUB_OFFLINE

CUDA_HOME=$(dirname "$(dirname "$(readlink -f "$(which nvcc)")")")
export CUDA_HOME
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${CUDA_HOME}/lib64:${LD_LIBRARY_PATH:-}"

cd "${LOCAL_REPO}"

python -m pip install -e . --no-deps

# =============================================================================
# Runtime preflight
# =============================================================================

echo "============================================================"
echo "GeneCAD three-species sequential full-genome evaluation"
echo "============================================================"
echo "Base model:       ${BASE_MODEL}"
echo "Head checkpoint:  ${HEAD_CHECKPOINT}"
echo "Input root:       ${INPUT_ROOT}"
echo "Output root:      ${OUTPUT_ROOT}"
echo "Prediction batch: ${PREDICT_BATCH_SIZE}"
echo "Execution order:  Ath full genome -> Mpo full genome -> Ppa full genome"
echo "============================================================"

python - <<'PY'
import torch

print("Torch:", torch.__version__)
print("Torch CUDA:", torch.version.cuda)
print("CUDA available:", torch.cuda.is_available())

assert torch.cuda.is_available(), "CUDA is unavailable"

print("GPU:", torch.cuda.get_device_name(0))
print("Capability:", torch.cuda.get_device_capability(0))
print("BF16:", torch.cuda.is_bf16_supported())
PY

# =============================================================================
# Run one species completely before starting the next
# =============================================================================

run_species() {
    local species_id="$1"
    local fasta="$2"
    local reference_gff="$3"

    local output_dir="${OUTPUT_ROOT}/${species_id}"
    local local_reference="${LOCAL_WORKDIR}/${species_id}.reference.gff3"

    local raw_gff="${output_dir}/${species_id}_GeneCAD_raw.gff"
    local orf_gff="${output_dir}/${species_id}_GeneCAD_orf.gff"
    local final_gff="${output_dir}/${species_id}_GeneCAD_final.gff"
    local report="${output_dir}/${species_id}_final_evaluation_no_busco.txt"

    echo
    echo "################################################################"
    echo "STARTING: ${species_id}"
    echo "################################################################"
    echo "FASTA:       ${fasta}"
    echo "Reference:   ${reference_gff}"
    echo "Base model:  ${BASE_MODEL}"
    echo "Checkpoint:  ${HEAD_CHECKPOINT}"
    echo "Output:      ${output_dir}"
    echo "Start time:  $(date)"
    echo "################################################################"

    if [[ "${ALLOW_RESUME}" != "yes" && -d "${output_dir}" ]]; then
        if find "${output_dir}" -mindepth 1 -print -quit | grep -q .; then
            echo "[FATAL] Output directory is not empty:"
            echo "        ${output_dir}"
            echo
            echo "Use a new OUTPUT_ROOT or submit with:"
            echo "  sbatch --export=ALL,ALLOW_RESUME=yes $0"
            return 1
        fi
    fi

    mkdir -p "${output_dir}"

    # -------------------------------------------------------------------------
    # Prediction
    # -------------------------------------------------------------------------

    bash predict.sh \
        --input-fasta "${fasta}" \
        --output-dir "${output_dir}" \
        --species-id "${species_id}" \
        --mode plant \
        --top-n-contigs all \
        --cpu-workers "${SLURM_CPUS_PER_TASK}" \
        --batch-size "${PREDICT_BATCH_SIZE}" \
        --gpus 0 \
        --model-checkpoint "${HEAD_CHECKPOINT}" \
        --base-model "${BASE_MODEL}"

    [[ -s "${raw_gff}" ]] || {
        echo "[FATAL] Missing raw GFF: ${raw_gff}"
        return 1
    }

    [[ -s "${orf_gff}" ]] || {
        echo "[FATAL] Missing ORF GFF: ${orf_gff}"
        return 1
    }

    [[ -s "${final_gff}" ]] || {
        echo "[FATAL] Missing final GFF: ${final_gff}"
        return 1
    }

    # -------------------------------------------------------------------------
    # Reference decompression
    # -------------------------------------------------------------------------

    if [[ "${reference_gff}" == *.gz ]]; then
        gzip -cd "${reference_gff}" > "${local_reference}"
    else
        cp "${reference_gff}" "${local_reference}"
    fi

    [[ -s "${local_reference}" ]] || {
        echo "[FATAL] Failed to prepare reference: ${reference_gff}"
        return 1
    }

    # -------------------------------------------------------------------------
    # Evaluation
    # -------------------------------------------------------------------------

    genecad evaluate \
        --ref "${local_reference}" \
        --pred "${final_gff}" \
        --fasta "${fasta}" \
        --skip-busco \
        --cpu "${SLURM_CPUS_PER_TASK}" \
        --output "${report}"

    [[ -s "${report}" ]] || {
        echo "[FATAL] Evaluation report was not created: ${report}"
        return 1
    }

    # -------------------------------------------------------------------------
    # Print this species' complete evaluation before moving to the next species
    # -------------------------------------------------------------------------

    echo
    echo "================================================================"
    echo "EVALUATION RESULT: ${species_id}"
    echo "================================================================"
    cat "${report}"
    echo "================================================================"
    echo "END OF RESULT: ${species_id}"
    echo "================================================================"
    echo "Raw GFF:    ${raw_gff}"
    echo "ORF GFF:    ${orf_gff}"
    echo "Final GFF:  ${final_gff}"
    echo "Report:     ${report}"
    echo "Finish time: $(date)"
    echo "[PASS] ${species_id} prediction and evaluation completed"
    echo "================================================================"
}

# =============================================================================
# Sequential execution
# =============================================================================

completed=0

for test_info in "${TESTS[@]}"; do
    IFS='|' read -r species_id fasta_path reference_path <<< "${test_info}"

    run_species \
        "${species_id}" \
        "${fasta_path}" \
        "${reference_path}"

    completed=$((completed + 1))

    echo
    echo "[PROGRESS] ${completed}/${#TESTS[@]} species completed"
    echo
done

# =============================================================================
# Final concise comparison
# =============================================================================

echo
echo "################################################################"
echo "FINAL COMPARISON: MAIN METRICS"
echo "################################################################"

for test_info in "${TESTS[@]}"; do
    IFS='|' read -r species_id fasta_path reference_path <<< "${test_info}"

    report="${OUTPUT_ROOT}/${species_id}/${species_id}_final_evaluation_no_busco.txt"

    echo
    echo "===== ${species_id} ====="

    grep -E -A4 \
        "CDS-exon-level|Locus-level|Transcript-level|Site +Ref" \
        "${report}" || true
done

echo
echo "################################################################"
echo "[PASS] All ${#TESTS[@]} species completed successfully"
echo "Output root: ${OUTPUT_ROOT}"
echo "Finish time: $(date)"
echo "################################################################"
