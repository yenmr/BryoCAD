#!/usr/bin/env bash
set -euo pipefail
for task in donor acceptor tis tts; do
  python motif_recovery/score_plantcad2_motif_recovery.py \
    --input "examples/motif_recovery/${task}.tsv" \
    --output "results/toy_motif_recovery/${task}.tsv" \
    --task "$task" --model "${MODEL:-bryocad/BryoCAD-Small}" \
    --batch-size 1 --num-workers 0
done
