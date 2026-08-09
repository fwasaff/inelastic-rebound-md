#!/usr/bin/env bash
# run_Nc_dumps.sh — corridas cortas con dump por-átomo para medir N_c(N).
# Todas a la misma altura de plateau (h0=41.47, v_imp≈1.47) para que N_c(N)
# se mida a velocidad de impacto constante. Salida en Nc_runs/N{N}/.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
OUTDIR="${SCRIPT_DIR}/data/outputs/terminado/Nc_runs"
LAMMPS_SCRIPT="${SCRIPT_DIR}/lammps/in.bounce_cluster"
CLUSTER_DIR="${SCRIPT_DIR}/lammps/clusters"
TABLE="${SCRIPT_DIR}/data/outputs/terminado/reprod_2012/table_repulsiva.txt"
H0=41.47
N_VALUES="${1:-7 13 55 79 87}"
JOBS="${2:-5}"

run_one() {
  local n="$1"
  local rundir="${OUTDIR}/N${n}"
  mkdir -p "${rundir}"
  cp "${TABLE}" "${rundir}/table_repulsiva.txt"
  cd "${rundir}" && \
  lmp -in "${LAMMPS_SCRIPT}" \
      -var cluster_file "${CLUSTER_DIR}/cluster_N${n}_min.dat" \
      -var h0 "${H0}" -var t_total 120.0 -var dump_freq 20 \
      -log "${rundir}/lammps.log" -screen none 2>"${rundir}/stderr.log" \
    && echo "  OK N=${n}" || echo "  [ERROR] N=${n}"
}
export -f run_one
export OUTDIR LAMMPS_SCRIPT CLUSTER_DIR TABLE H0
printf '%s\n' $N_VALUES | xargs -P "${JOBS}" -I{} bash -c 'run_one {}'
echo "N_c dumps listos en ${OUTDIR}"
