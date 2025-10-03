#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT="${RSCRIPT:-$(command -v Rscript || true)}"
if [[ -z "${RSCRIPT}" ]]; then
  echo "ERROR: Rscript not found in PATH" >&2
  exit 127
fi

exec "${RSCRIPT}" "${SCRIPT_DIR}/bayesCPQuery.r" "$@"