#!/usr/bin/env bash
# Read-only to the original checkout. Creates a NEW review directory.
set -euo pipefail
if [[ $# != 2 ]]; then
  echo 'Usage: bash capture_upstream_state.sh /path/to/genecad /path/to/new-review-directory' >&2
  exit 2
fi
repo=$(cd "$1" && pwd)
out=$2
if [[ -e "$out" ]]; then echo 'Output already exists; use a new directory.' >&2; exit 2; fi
mkdir -p "$out"
out=$(cd "$out" && pwd)
pin=1a9baf84d61152c20a0cd910694316cd88d1d374
if git -C "$repo" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  git -C "$repo" rev-parse HEAD > "$out/HEAD.txt"
  git -C "$repo" status --short > "$out/status.txt"
  if git -C "$repo" cat-file -e "$pin^{commit}" 2>/dev/null; then
    git -C "$repo" diff "$pin" -- train.sh predict.sh src scripts pyproject.toml > "$out/changes_vs_v0.4.0.patch"
  else
    echo 'Pinned upstream commit absent locally; no fetch performed.' > "$out/diff_unavailable.txt"
  fi
else
  echo 'No Git metadata; working source files captured for comparison.' > "$out/no_git_metadata.txt"
fi
# Capture only implementation files; no environments, data or model weights.
python3 - "$repo" "$out" <<'PY'
import hashlib,json,shutil,sys
from pathlib import Path
repo,out=map(Path,sys.argv[1:]);records=[]
paths=[repo/p for p in ['train.sh','predict.sh','pyproject.toml','LICENSE']]
for folder in ['src','scripts']:
    paths.extend((repo/folder).rglob('*.py'))
    paths.extend((repo/folder).rglob('*.sh'))
for p in sorted(set(paths)):
    if not p.is_file():continue
    rel=p.relative_to(repo);dst=out/'source'/rel;dst.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(p,dst)
    records.append({'file':str(rel),'sha256':hashlib.sha256(dst.read_bytes()).hexdigest()})
(out/'source_hashes.json').write_text(json.dumps(records,indent=2)+'\n')
print(f'Captured {len(records)} source files for local review: {out}')
PY
