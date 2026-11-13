#!/bin/bash
#SBATCH --job-name=genecad_2bryo
#SBATCH --partition=b6000-al9
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=5-00:00:00
#SBATCH --output=genecad_plantcad2_2bryo_%j.out
#SBATCH --error=genecad_plantcad2_2bryo_%j.err

set -euo pipefail

# =============================================================================
# Paths
# =============================================================================

ENV_TAR="/ceph/sharedfs/work/TYW_Lab/yenmr/envs/genecad_b600_cuda128_torch271.tar.gz"

GENECAD_SOURCE="/dicos_ui_home/yenmr/projects/BryoCAD/annotation/genecad"

HF_CACHE="/ceph/sharedfs/work/TYW_Lab/yenmr/huggingface_cache"

GENOME_DIR="/dicos_ui_home/yenmr/projects/BryoCAD/warehouse/library/genomes3"

SPECIES_CONFIG="${GENOME_DIR}/bryophytes_2train_1eval.yaml"

OUTPUT_DIR="/dicos_ui_home/yenmr/projects/BryoCAD/annotation/genecad_output/bryophyte_training/plantcad2_small_d768_2bryo_v1"

BASE_MODEL="kuleshov-group/PlantCAD2-Small-l24-d0768"

TRAINING_SPECIES="AagrBONN MpTak"

# =============================================================================
# Training parameters
# =============================================================================

EPOCHS=5
BATCH_SIZE=8
EFFECTIVE_BATCH=40
LEARNING_RATE="1e-4"

# =============================================================================
# Input validation
# =============================================================================

for path in \
    "${ENV_TAR}" \
    "${GENECAD_SOURCE}/pyproject.toml" \
    "${GENECAD_SOURCE}/train.sh" \
    "${SPECIES_CONFIG}"
do
    [[ -e "${path}" ]] || {
        echo "[FATAL] Missing: ${path}"
        exit 1
    }
done

# This run must start from a new GeneCAD head.
if [[ -e "${OUTPUT_DIR}/pipeline" ||
      -e "${OUTPUT_DIR}/checkpoints" ||
      -e "${OUTPUT_DIR}/training.log" ]]; then
    echo "[FATAL] Output directory already contains training artifacts:"
    echo "        ${OUTPUT_DIR}"
    echo
    echo "Use a new OUTPUT_DIR."
    echo "This job will not resume an earlier checkpoint."
    exit 1
fi

# =============================================================================
# Local runtime
# =============================================================================

LOCAL_WORKDIR=$(mktemp -d -p /tmp "${USER}_genecad_2bryo_XXXXXX")
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

mkdir -p \
    "${LOCAL_ENV}" \
    "${OUTPUT_DIR}" \
    "${HF_CACHE}"

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

# Register GeneCAD without reinstalling dependencies.
python -m pip install -e . --no-deps

# =============================================================================
# Verify the source modifications
# =============================================================================

echo "============================================================"
echo "Checking GeneCAD legacy-head settings"
echo "============================================================"

grep -n \
    -E 'TOKEN_EMBED_DIM=|BASE_FROZEN=|learning-rate-decay|checkpoint-frequency' \
    train.sh

grep -n 'weight_decay=' src/modeling.py

grep -q 'TOKEN_EMBED_DIM=128' train.sh || {
    echo "[FATAL] train.sh does not use TOKEN_EMBED_DIM=128"
    exit 1
}

grep -q 'weight_decay=1e-5' src/modeling.py || {
    echo "[FATAL] src/modeling.py does not use weight_decay=1e-5"
    exit 1
}

grep -Eq -- '--learning-rate-decay[[:space:]]+none' train.sh || {
    echo "[FATAL] train.sh does not use constant learning rate"
    exit 1
}

# =============================================================================
# Environment check
# =============================================================================

python - <<'PY'
import torch

print("Torch:", torch.__version__)
print("Torch CUDA:", torch.version.cuda)
print("CUDA available:", torch.cuda.is_available())

assert torch.cuda.is_available(), "CUDA is unavailable"
assert torch.cuda.get_device_capability(0) == (12, 0)

print("GPU:", torch.cuda.get_device_name(0))
print("Capability:", torch.cuda.get_device_capability(0))
print("BF16:", torch.cuda.is_bf16_supported())
PY

# =============================================================================
# Training summary
# =============================================================================

echo "============================================================"
echo "GeneCAD 2-bryophyte training"
echo "============================================================"
echo "Source commit:         $(git rev-parse HEAD 2>/dev/null || echo unknown)"
echo "Species config:        ${SPECIES_CONFIG}"
echo "Training species:      ${TRAINING_SPECIES}"
echo "Base model:            ${BASE_MODEL}"
echo "Base frozen:           yes"
echo "Token embedding:       128"
echo "Calculated head dim:   768"
echo "Calculated projection: 640"
echo "Head layers:           8"
echo "Attention heads:       8"
echo "Epochs:                ${EPOCHS}"
echo "Batch size:            ${BATCH_SIZE}"
echo "Effective batch:       ${EFFECTIVE_BATCH}"
echo "Gradient accumulation: $((EFFECTIVE_BATCH / BATCH_SIZE))"
echo "Learning rate:         ${LEARNING_RATE}"
echo "Learning-rate decay:   none"
echo "Weight decay:          1e-5"
echo "Auto class weights:    yes"
echo "Output:                 ${OUTPUT_DIR}"
echo "============================================================"

# =============================================================================
# Train from a newly initialized GeneCAD head
#
# No checkpoint is supplied here.
# BASE_MODEL loads only the frozen pretrained PlantCAD2 backbone.
# =============================================================================

bash train.sh \
    --domain plant \
    --base-model "${BASE_MODEL}" \
    --base-frozen yes \
    --species-config "${SPECIES_CONFIG}" \
    --species-ids "${TRAINING_SPECIES}" \
    --plant-input-dir "${GENOME_DIR}" \
    --plant-gff-dir "${GENOME_DIR}" \
    --output "${OUTPUT_DIR}" \
    --run-name plantcad2-small-d768-2bryo-v1 \
    --project genecad-bryophytes \
    --gpus 1 \
    --batch-size "${BATCH_SIZE}" \
    --effective "${EFFECTIVE_BATCH}" \
    --lr "${LEARNING_RATE}" \
    --epochs "${EPOCHS}" \
    --auto-class-weights yes \
    --window-size 8192 \
    --intergenic-proportion 0.15 \
    --gff-parallel 1

# =============================================================================
# Validate output
# =============================================================================

BEST_CHECKPOINT="${OUTPUT_DIR}/checkpoints/best.ckpt"
SELECTION_JSON="${OUTPUT_DIR}/qc/best_checkpoint_summary.json"
METRICS_CSV=$(find "${OUTPUT_DIR}/logs/csv" -name metrics.csv -type f |
    sort |
    tail -n 1)

[[ -e "${BEST_CHECKPOINT}" ]] || {
    echo "[FATAL] best.ckpt was not created"
    exit 1
}

[[ -s "${SELECTION_JSON}" ]] || {
    echo "[FATAL] Missing checkpoint selection summary"
    exit 1
}

[[ -n "${METRICS_CSV}" && -s "${METRICS_CSV}" ]] || {
    echo "[FATAL] metrics.csv was not found"
    exit 1
}

# =============================================================================
# Confirm head architecture
# =============================================================================

echo "============================================================"
echo "Checking trained checkpoint architecture"
echo "============================================================"

CHECKPOINT="${BEST_CHECKPOINT}" python - <<'PY'
import os
import torch

checkpoint = torch.load(
    os.environ["CHECKPOINT"],
    map_location="cpu",
    weights_only=False,
)
state = checkpoint["state_dict"]

expected = {
    "token_embedding.weight": (16, 128),
    "embedding_projection.weight": (640, 1536),
    "head_encoder.layers.0.attn.Wqkv.weight": (2304, 768),
    "classifier.Wi.weight": (6144, 768),
    "classifier.Wo.weight": (17, 3072),
}

for name, expected_shape in expected.items():
    actual_shape = tuple(state[name].shape)
    print(f"{name}: {actual_shape}")
    assert actual_shape == expected_shape, (
        f"{name}: expected {expected_shape}, got {actual_shape}"
    )

print("Checkpoint epoch:", checkpoint.get("epoch"))
print("Global step:", checkpoint.get("global_step"))
print("[PASS] Checkpoint uses the legacy d768 GeneCAD head")
PY

# =============================================================================
# Print validation history
# =============================================================================

echo "============================================================"
echo "Validation history"
echo "============================================================"

METRICS_CSV="${METRICS_CSV}" python - <<'PY'
import os
import pandas as pd

path = os.environ["METRICS_CSV"]
df = pd.read_csv(path)

columns = [
    "epoch",
    "step",
    "train/loss_epoch",
    "valid/loss_epoch",
    "valid__token__overall/f1",
    "valid__entity__overall/f1",
]

available = [column for column in columns if column in df.columns]
summary = df.loc[
    df["valid__entity__overall/f1"].notna(),
    available,
].copy()

summary = summary.drop_duplicates(
    subset=["epoch"],
    keep="last",
)

print(summary.to_string(index=False))
PY

echo "============================================================"
echo "Best checkpoint selection"
echo "============================================================"
cat "${SELECTION_JSON}"

echo
echo "============================================================"
echo "Saved checkpoints"
echo "============================================================"
find "${OUTPUT_DIR}/checkpoints" \
    -maxdepth 1 \
    -type f \
    -printf '%f %s bytes\n' |
    sort

echo "============================================================"
echo "[PASS] GeneCAD 2-bryophyte training completed"
echo "Best checkpoint: ${BEST_CHECKPOINT}"
echo "Metrics:         ${METRICS_CSV}"
echo "Output:          ${OUTPUT_DIR}"
echo "============================================================"
