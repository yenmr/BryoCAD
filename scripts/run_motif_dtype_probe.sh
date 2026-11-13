#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
PYTHON=${PYTHON:-python}
export PYTHONFAULTHANDLER=1 PYTHONUNBUFFERED=1
OUT=results/motif_dtype_probe
mkdir -p "$OUT"
"$PYTHON" scripts/check_environment.py --require-cuda --profile zero-shot --output "$OUT/environment.json"
BASE=$("$PYTHON" scripts/download_test_model.py BryoCAD-Small)
printf 'dtype\texit_code\n' > "$OUT/status.tsv"
failed=0
for dtype in float32 float16; do
  if "$PYTHON" motif_recovery/score_plantcad2_motif_recovery.py \
    --input examples/motif_recovery/donor.tsv --output "$OUT/$dtype.tsv" \
    --task donor --model "$BASE" --batch-size 1 --num-workers 0 \
    --max-samples 1 --device cuda:0 --dtype "$dtype" > "$OUT/$dtype.log" 2>&1; then
    printf '%s\t0\n' "$dtype" >> "$OUT/status.tsv"
  else
    code=$?
    printf '%s\t%s\n' "$dtype" "$code" >> "$OUT/status.tsv"
    failed=1
  fi
done
cat "$OUT/status.tsv"
exit "$failed"
