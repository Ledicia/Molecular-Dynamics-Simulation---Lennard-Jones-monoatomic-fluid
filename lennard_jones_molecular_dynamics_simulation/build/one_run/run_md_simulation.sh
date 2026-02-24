#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BIN_DIR="${ROOT_DIR}/build/one_run/bin"
OUT_DIR="${ROOT_DIR}/outputs"

EXE="${BIN_DIR}/md_simulation_program"

if [[ ! -x "${EXE}" ]]; then
  echo "Executable not found: ${EXE}"
  echo "Run: build/one_run/compile.sh"
  exit 1
fi

if [[ ! -f "${OUT_DIR}/rv_init.dat" ]]; then
  echo "Missing outputs/rv_init.dat"
  echo "Run: build/one_run/run_initial_config.sh"
  exit 1
fi

cd "${ROOT_DIR}"
"${EXE}"

echo "=> molecular_dynamics_simulation finished."
echo "=> Check outputs/ for generated files."
