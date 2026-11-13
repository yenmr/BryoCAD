#!/usr/bin/env bash
set -euo pipefail
# Tiny inference checks use one GPU. Keep the first scheduler-assigned device.
# An explicitly empty mask remains empty (do not escape a CPU-only allocation).
if [[ ${CUDA_VISIBLE_DEVICES+x} ]]; then
  export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES%%,*}"
else
  export CUDA_VISIBLE_DEVICES=0
fi
ROOT=$(cd "$(dirname "$0")/.." && pwd)
INPUT=${1:-"$ROOT/examples/atac/input_template.tsv"}
# Resolve user-supplied relative paths before changing directory.
PYTHON=${PYTHON:-python}
INPUT=$("$PYTHON" -c 'import sys; from pathlib import Path; print(Path(sys.argv[1]).resolve())' "$INPUT")
cd "$ROOT"
export PYTHONFAULTHANDLER=1 PYTHONUNBUFFERED=1
OUT=results/atac_predict_check
mkdir -p "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'code=$?; printf "%s\n" "$code" > "$OUT/exit_code.txt"' EXIT
"$PYTHON" scripts/check_environment.py --require-cuda --profile atac --output "$OUT/environment.json"
"$PYTHON" -c 'import fire, datasets, sklearn, scipy'
BASE=$("$PYTHON" scripts/download_test_model.py BryoCAD-Small)
ADAPTER=$("$PYTHON" scripts/download_test_model.py BryoCAD-Small-ATAC-Mpolymorpha-Ctrl-LoRA)
"$PYTHON" chromatin_accessibility/lora_fine_tune.py predict \
  --checkpoint_dir="$ADAPTER" --model_name="$BASE" --data_dir="$INPUT" \
  --output_file="$OUT/predictions.tsv" --task_type=classification \
  --batch_size=1 --bf16=False
"$PYTHON" scripts/check_atac_predictions.py "$INPUT" "$OUT/predictions.tsv"
