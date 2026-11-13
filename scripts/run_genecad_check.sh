#!/usr/bin/env bash
set -euo pipefail
export PYTHONFAULTHANDLER=1 PYTHONUNBUFFERED=1
ROOT=$(cd "$(dirname "$0")/.." && pwd)
GENECAD_DIR=${1:?Usage: bash scripts/run_genecad_check.sh /absolute/path/to/patched-genecad}
GENECAD_DIR=$(cd "$GENECAD_DIR" && pwd)
[[ -f "$GENECAD_DIR/predict.sh" ]] || { echo 'predict.sh not found'; exit 1; }
source "$ROOT/scripts/select_genecad_python.sh"
PYTHON=$(select_genecad_python "$GENECAD_DIR")
export GENECAD_PYTHON="$PYTHON"
cd "$ROOT"
OUT="$ROOT/results/genecad_check"
mkdir -p "$OUT"
if [[ -d "$OUT/predictions" ]]; then
  mv "$OUT/predictions" "$OUT/predictions.previous.$(date +%Y%m%d_%H%M%S).$$"
fi
exec > >(tee "$OUT/run.log") 2>&1
TEMP_PREDICT=""
trap 'code=$?; if [[ -n "$TEMP_PREDICT" ]]; then rm -f "$TEMP_PREDICT"; fi; printf "%s\n" "$code" > "$OUT/exit_code.txt"' EXIT
echo "GeneCAD interpreter: $PYTHON"
"$PYTHON" scripts/check_genecad_environment.py --checkout "$GENECAD_DIR" --output "$OUT/genecad_environment.json"
"$PYTHON" scripts/check_environment.py --output "$OUT/environment.json"
HEAD=$("$PYTHON" scripts/download_test_model.py GeneCAD-BryoCAD-Small-4bryo)
BASE=$("$PYTHON" scripts/download_test_model.py BryoCAD-Small)
export GENECAD_PYTHON
GENECAD_PYTHON=$(command -v "$PYTHON")
cd "$GENECAD_DIR"
# Upstream predict.sh hard-codes BF16. Use a temporary wrapper in the same
# directory, retaining SCRIPT_DIR and the original checkout unchanged.
TEMP_PREDICT=$(mktemp "$GENECAD_DIR/.bryocad_predict.XXXXXX.sh")
trap 'code=$?; rm -f "$TEMP_PREDICT"; printf "%s\n" "$code" > "$OUT/exit_code.txt"' EXIT
"$PYTHON" "$ROOT/scripts/prepare_genecad_predict.py" "$GENECAD_DIR/predict.sh" "$TEMP_PREDICT" --dtype "${GENECAD_DTYPE:-float32}"
bash "$TEMP_PREDICT" --input-fasta "$ROOT/examples/gene_annotation/input.fa.gz" \
  --output-dir "$OUT/predictions" --species-id toy_mpolymorpha \
  --mode plant --top-n-contigs all --cpu-workers 2 --batch-size 1 --gpus 0 \
  --model-checkpoint "$HEAD/best.ckpt" --base-model "$BASE"
"$PYTHON" "$ROOT/scripts/check_genecad_output.py" "$OUT/predictions" --fasta "$ROOT/examples/gene_annotation/input.fa.gz"
echo 'PASS: GeneCAD execution check; inspect results/genecad_check/'
