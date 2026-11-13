#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
PYTHON=${PYTHON:-python3}
"$PYTHON" chromatin_accessibility/test_forward_replication.py
"$PYTHON" chromatin_accessibility/test_atac_tsv.py
"$PYTHON" scripts/test_cuda_preflight.py
"$PYTHON" scripts/test_peft_forward_compat.py
"$PYTHON" variant_analysis/test_snp_scoring.py
"$PYTHON" orthogroups/test_extract_proteins.py
"$PYTHON" examples/motif_recovery/check_examples.py
"$PYTHON" scripts/check_toy_inputs.py
"$PYTHON" scripts/check_environment.py
if "$PYTHON" -c 'import vcf' 2>/dev/null; then
  "$PYTHON" variant_analysis/test_snp_vcf_io.py
else
  echo 'Optional VCF I/O test skipped: install PyVCF3 to enable it.'
fi
