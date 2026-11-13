#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
PYTHON=${PYTHON:-python}
export PYTHONFAULTHANDLER=1 PYTHONUNBUFFERED=1
# Optional first argument: patched GeneCAD checkout. Individual tests remain usable.
mkdir -p results/release_checks
SUMMARY="$ROOT/results/release_checks/status.tsv"
printf 'test\texit_code\tstatus\n' > "$SUMMARY"
failed=0
run_check() {
  local name=$1
  shift
  if "$@"; then
    printf '%s\t0\tPASS\n' "$name" >> "$SUMMARY"
  else
    local code=$?
    printf '%s\t%s\tFAIL\n' "$name" "$code" >> "$SUMMARY"
    failed=1
  fi
}
run_check zero-shot env PYTHON="${BRYOCAD_PYTHON:-$PYTHON}" CHECK_OUTPUT_DIR="$ROOT/results/zero_shot_check" bash scripts/run_zero_shot_check.sh
run_check atac env PYTHON="${BRYOCAD_PYTHON:-$PYTHON}" CHECK_OUTPUT_DIR="$ROOT/results/atac_check" bash scripts/run_atac_check.sh
if [[ $# -ge 1 ]]; then
  run_check genecad bash scripts/run_genecad_check.sh "$1"
else
  printf 'genecad\tNA\tNOT_REQUESTED (supply GeneCAD directory to include)\n' >> "$SUMMARY"
fi
cat "$SUMMARY"
exit "$failed"
