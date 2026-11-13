#!/bin/bash
#SBATCH --job-name=LoRA_TIS
#SBATCH --partition=b6000-al9
##SBATCH --partition=a100-al9
#SBATCH --gres=gpu:1                 
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
##SBATCH --time=24:00:00
#SBATCH --time=5-00:00:00
#SBATCH -o BryoCAD_lofa_tis_512.out


set -eo pipefail

# ==========================================
# 1. Define shared-storage and local paths.
# ==========================================
# Source paths on shared storage.
NAS_ENV_TAR="/ceph/sharedfs/work/TYW_Lab/yenmr/envs/plantcad2_evo2_b600_sm120_cuda130.tar.gz"
#NAS_DATA_DIR="/ceph/sharedfs/work/TYW_Lab/yenmr/projects/BryoCAD/warehouse/training_data"
#NAS_OUTPUT_DIR="/ceph/sharedfs/work/TYW_Lab/yenmr/projects/BryoCAD/warehouse/models/bryocad_large_v1"

# Create a temporary workspace on local storage.
LOCAL_WORKDIR=$(mktemp -d -p /tmp yenmr_bryocad_XXXXXX)

# ==========================================
# 2. Register cleanup on exit.
# ==========================================
# Run cleanup when the shell exits.
trap 'echo "[CLEANUP] Removing local temporary directory: $LOCAL_WORKDIR..."; rm -rf "$LOCAL_WORKDIR"; echo "[CLEANUP] Done."' EXIT

echo "[INFO] Job Started on Node: $SLURMD_NODENAME"
echo "[INFO] Local Workspace: $LOCAL_WORKDIR"

# ==========================================
# 3. Copy and extract inputs.
# ==========================================
echo "[INFO] Copying and extracting the conda environment..."
tar -xzf "$NAS_ENV_TAR" -C "$LOCAL_WORKDIR"
source "$LOCAL_WORKDIR/bin/activate"
conda-unpack

#echo "[INFO] Copying the Arrow dataset to local storage (/tmp)..."
# Use rsync to preserve the directory structure.
#rsync -a "$NAS_DATA_DIR/" "$LOCAL_WORKDIR/dataset/"

# ==========================================
# 4. Configure hardware settings.
# ==========================================
# Enable TF32 through environment variables.
export TORCH_CUDNN_V8_API_ENABLED=1
export TORCH_CUDA_ARCH_LIST="8.0;12.0" # Target A100 (8.0) and Blackwell (12.0).

# Disable online Hugging Face lookups.
export HF_DATASETS_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export PYTHONUNBUFFERED=1

# ==========================================
# 5. Launch the Python task.
# ==========================================
echo "[INFO] Launching the Python script..."


python lora_fine_tune.py train \
    --train_dir="BryoCAD_TIS_LoRA_input/Train_512.parquet" \
    --valid_dir="BryoCAD_TIS_LoRA_input/Validation_512.parquet" \
    --model_name="../warehouse/models/bryocad_small_final/" \
    --task_type="classification" \
    --output_dir="./tis_fine_tuned_lora_checkpoint_512" \
    --train_batch_size=8 \
    --eval_batch_size=32 \
    --gradient_accumulation_steps=16 \
    --max_steps=2000 \
    --eval_steps=200 \
    --save_steps=200 \
    --bf16=True


echo "[SUCCESS] Training session finished."
