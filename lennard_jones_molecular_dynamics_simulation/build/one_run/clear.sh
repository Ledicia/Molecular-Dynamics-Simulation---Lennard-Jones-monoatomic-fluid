#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

"${ROOT_DIR}/build/one_run/compile.sh"
"${ROOT_DIR}/build/one_run/run_initial_config.sh"
"${ROOT_DIR}/build/one_run/run_md_simulation.sh"
"${ROOT_DIR}/build/one_run/run_results.sh"

echo "=> FULL PIPELINE finished."
