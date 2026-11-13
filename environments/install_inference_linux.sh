#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
exec bash "$ROOT/environments/unified/build_conda_env.sh" "$@"
