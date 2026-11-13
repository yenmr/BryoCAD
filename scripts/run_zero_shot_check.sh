#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
PYTHON=${PYTHON:-python}
export PYTHONFAULTHANDLER=1 PYTHONUNBUFFERED=1
OUT=${CHECK_OUTPUT_DIR:-results/zero_shot_check}
mkdir -p "$OUT/motif"
exec > >(tee "$OUT/run.log") 2>&1
trap 'code=$?; printf "%s\n" "$code" > "$OUT/exit_code.txt"' EXIT
"$PYTHON" scripts/check_environment.py --require-cuda --profile zero-shot --output "$OUT/environment.json"
cp environments/model_revisions.json "$OUT/model_revisions.json"
BASE=$("$PYTHON" scripts/download_test_model.py BryoCAD-Small)
# FP32 is intentional: do not infer Triton FP16 compatibility from GPU capability.
for task in donor acceptor tis tts; do
  "$PYTHON" motif_recovery/score_plantcad2_motif_recovery.py \
    --input "examples/motif_recovery/$task.tsv" --output "$OUT/motif/$task.tsv" \
    --task "$task" --model "$BASE" --batch-size 1 --num-workers 0 \
    --device cuda:0 --dtype "${MOTIF_DTYPE:-float32}"
done
"$PYTHON" examples/motif_recovery/check_examples.py --results "$OUT/motif"
"$PYTHON" variant_analysis/zero_shot_score.py \
  -input-vcf examples/snp/input.vcf -input-fasta examples/gene_annotation/input.fa.gz \
  -model "$BASE" -output "$OUT/snps.vcf" -device cuda:0 -batchSize 1 -contextSize 8192
"$PYTHON" scripts/check_snp_output.py "$OUT/snps.vcf" --input examples/snp/input.vcf
echo "PASS: zero-shot checks; output: $OUT"
