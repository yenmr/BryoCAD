#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
PYTHON="${PYTHON:-python3}"
case "${1:-}" in
  model-organisms)
    "$PYTHON" figures/summarize_model_organisms.py
    Rscript figures/plot_model_organisms.R ;;
  lineage-transfer)
    Rscript figures/plot_lineage_transfer.R ;;
  bryophyte-tts)
    Rscript motif_recovery/tts/analyze_tts.R figures/source_data/bryophyte_recovery results/tts_bryophytes BryoCAD-Small ;;
  angiosperm-tts-panel)
    Rscript motif_recovery/tts/assemble_tts_figure.R figures/source_data/angiosperm_tts results/tts_angiosperms PlantCAD2-Small ;;
  gene-annotation)
    Rscript figures/plot_gene_annotation_coding.R ;;
  *) echo 'Usage: bash figures/run_selected.sh {model-organisms|lineage-transfer|bryophyte-tts|angiosperm-tts-panel|gene-annotation}' >&2; exit 2 ;;
esac
