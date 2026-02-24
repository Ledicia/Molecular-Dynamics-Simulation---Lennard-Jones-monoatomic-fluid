#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BIN_DIR="${ROOT_DIR}/build/one_run/bin"
OUT_DIR="${ROOT_DIR}/outputs"

EXE="${BIN_DIR}/md_initial_config_program"

mkdir -p "${OUT_DIR}"

if [[ ! -x "${EXE}" ]]; then
  echo "Executable not found: ${EXE}"
  echo "Run: build/one_run/compile.sh"
  exit 1
fi

cd "${ROOT_DIR}"
"${EXE}"

echo "=> initial_config finished."
echo "=> Expected: outputs/rv_init.dat"
