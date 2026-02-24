#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

ONE_RUN_DIR="${ROOT_DIR}/build/one_run"
OBJ_DIR="${ONE_RUN_DIR}/obj"
MOD_DIR="${ONE_RUN_DIR}/mod"
BIN_DIR="${ONE_RUN_DIR}/bin"

SRC_BASE="${ROOT_DIR}/scripts/base"
SRC_PHYS="${ROOT_DIR}/scripts/physics"
SRC_STATS="${ROOT_DIR}/scripts/stats"

# Programas principales (ajusta rutas si en tu repo están en otro sitio)
MAIN_INIT="${ROOT_DIR}/scripts/md_initial_config_program.f90"
MAIN_MD="${ROOT_DIR}/scripts/md_simulation_program.f90"

FC="${FC:-gfortran}"

FFLAGS=(
  -O2
  -ffree-form
  -ffree-line-length-none
  -Wall -Wextra
  -Wno-line-truncation
  -Wno-error=line-truncation
  -fimplicit-none
  -J "${MOD_DIR}"
  -I "${MOD_DIR}"
)

mkdir -p "${OBJ_DIR}" "${MOD_DIR}" "${BIN_DIR}"

echo "[build/one_run] Cleaning old .o/.mod..."
rm -f "${OBJ_DIR}"/*.o "${MOD_DIR}"/*.mod 2>/dev/null || true

echo "[build/one_run] Compiling shared modules..."

# --- base (ORDER MATTERS!)
"${FC}" "${FFLAGS[@]}" -c "${SRC_BASE}/define_precision.f90"  -o "${OBJ_DIR}/define_precision.o"
"${FC}" "${FFLAGS[@]}" -c "${SRC_BASE}/random_numbers.f90"   -o "${OBJ_DIR}/random_numbers.o"
"${FC}" "${FFLAGS[@]}" -c "${SRC_BASE}/md_types.f90"         -o "${OBJ_DIR}/md_types.o"

if [[ ! -f "${MOD_DIR}/md_types.mod" ]]; then
  echo "[build/one_run][ERROR] md_types.mod not generated."
  echo "  -> Check ${SRC_BASE}/md_types.f90: module name must be: module md_types"
  exit 1
fi

"${FC}" "${FFLAGS[@]}" -c "${SRC_BASE}/read_input_files.f90" -o "${OBJ_DIR}/read_input_files.o"

# --- physics
"${FC}" "${FFLAGS[@]}" -c "${SRC_PHYS}/geometry_pbc.f90"        -o "${OBJ_DIR}/geometry_pbc.o"
"${FC}" "${FFLAGS[@]}" -c "${SRC_PHYS}/lj_potential_energy.f90" -o "${OBJ_DIR}/lj_potential_energy.o"
"${FC}" "${FFLAGS[@]}" -c "${SRC_PHYS}/verlet.f90"              -o "${OBJ_DIR}/verlet.o"
"${FC}" "${FFLAGS[@]}" -c "${SRC_PHYS}/thermodynamic_coefs.f90" -o "${OBJ_DIR}/thermodynamic_coefs.o"

# --- stats
"${FC}" "${FFLAGS[@]}" -c "${SRC_STATS}/stats_math.f90"      -o "${OBJ_DIR}/stats_math.o"
"${FC}" "${FFLAGS[@]}" -c "${SRC_STATS}/md_means.f90"        -o "${OBJ_DIR}/md_means.o"
"${FC}" "${FFLAGS[@]}" -c "${SRC_STATS}/md_correlations.f90" -o "${OBJ_DIR}/md_correlations.o"

OBJS=(
  "${OBJ_DIR}/define_precision.o"
  "${OBJ_DIR}/random_numbers.o"
  "${OBJ_DIR}/md_types.o"
  "${OBJ_DIR}/read_input_files.o"
  "${OBJ_DIR}/geometry_pbc.o"
  "${OBJ_DIR}/lj_potential_energy.o"
  "${OBJ_DIR}/verlet.o"
  "${OBJ_DIR}/thermodynamic_coefs.o"
  "${OBJ_DIR}/stats_math.o"
  "${OBJ_DIR}/md_means.o"
  "${OBJ_DIR}/md_correlations.o"
)

echo "[build/one_run] Linking md_initial_config_program -> ${BIN_DIR}/md_initial_config_program"
"${FC}" "${FFLAGS[@]}" "${OBJS[@]}" "${MAIN_INIT}" -o "${BIN_DIR}/md_initial_config_program"

echo "[build/one_run] Linking md_simulation_program -> ${BIN_DIR}/md_simulation_program"
"${FC}" "${FFLAGS[@]}" "${OBJS[@]}" "${MAIN_MD}" -o "${BIN_DIR}/md_simulation_program"


echo "[build/one_run] OK -> ${BIN_DIR}/ (md_initial_config_program, md_simulation_program)"
