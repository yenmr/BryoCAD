#!/usr/bin/env bash
# Sourced by the GeneCAD runner; independent of BryoCAD's PYTHON override.
select_genecad_python() {
  local checkout=$1 candidate
  if [[ -n ${GENECAD_PYTHON:-} ]]; then
    candidate=$GENECAD_PYTHON
  elif [[ -x "$checkout/.venv/bin/python" ]]; then
    candidate="$checkout/.venv/bin/python"
  else
    candidate=python
  fi
  command -v "$candidate" || {
    echo "GeneCAD Python not found: $candidate. Set GENECAD_PYTHON to its environment's executable." >&2
    return 1
  }
}
