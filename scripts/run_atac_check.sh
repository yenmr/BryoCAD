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
cd "$ROOT"
PYTHON=${PYTHON:-python}
export PYTHONFAULTHANDLER=1 PYTHONUNBUFFERED=1
OUT=${CHECK_OUTPUT_DIR:-results/atac_check}
mkdir -p "$OUT"
exec > >(tee "$OUT/run.log") 2>&1
trap 'code=$?; printf "%s\n" "$code" > "$OUT/exit_code.txt"' EXIT
"$PYTHON" scripts/check_environment.py --require-cuda --profile atac --output "$OUT/environment.json"
cp environments/model_revisions.json "$OUT/model_revisions.json"
BASE=$("$PYTHON" scripts/download_test_model.py BryoCAD-Small)
ADAPTER=$("$PYTHON" scripts/download_test_model.py BryoCAD-Small-ATAC-Mpolymorpha-Ctrl-LoRA)
"$PYTHON" scripts/check_atac_loading.py --base "$BASE" --adapter "$ADAPTER" --output "$OUT/atac.json"
# Exercise the original prediction CLI as well, using readable raw TSV input.
bash scripts/run_atac_predict_check.sh "${1:-$ROOT/examples/atac/input_template.tsv}"
echo "PASS: ATAC restoration and original predict CLI; output: $OUT and results/atac_predict_check"
