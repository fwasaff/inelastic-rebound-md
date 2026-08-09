#!/usr/bin/env bash
# run_orient_sweep.sh — sweep N × altura × orientación para barras de error.
#
# Responde al punto M3 del arbitraje: cada (N, h0) se corre para 5
# orientaciones iniciales del cúmulo (una cristalina + 4 rotaciones
# aleatorias fijas), de modo que cada fracción α se reporta como media ± σ.
#
# Salida: data/outputs/terminado/cluster_orient/N{N}/h0_{h0}/o{k}/bounce.dat

set -euo pipefail

N_VALUES="1 2 7 13 55 79 87"
H0_VALUES="11.03 18.64 28.15 41.47 58.6 79.53"
JOBS=8
T_TOTAL=200.0
SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
OUTDIR="${SCRIPT_DIR}/data/outputs/terminado/cluster_orient"
LAMMPS_SCRIPT="${SCRIPT_DIR}/lammps/in.bounce_cluster"
CLUSTER_DIR="${SCRIPT_DIR}/lammps/clusters"
TABLE="${SCRIPT_DIR}/data/outputs/terminado/reprod_2012/table_repulsiva.txt"

# Ángulos de rotación (grados) por orientación: ax ay az. Semilla fija.
ANGLES=(
  "0.0 0.0 0.0"
  "327.646 128.500 20.874"
  "216.787 127.396 282.301"
  "278.714 268.623 312.305"
  "152.471 211.222 295.981"
)

while [[ $# -gt 0 ]]; do
  case "$1" in
    --n)    N_VALUES="$2";  shift 2 ;;
    --h0)   H0_VALUES="$2"; shift 2 ;;
    --jobs) JOBS="$2";      shift 2 ;;
    *) echo "Argumento desconocido: $1"; exit 1 ;;
  esac
done

run_sim() {
  local n="$1" h0="$2" k="$3" ax="$4" ay="$5" az="$6"
  local cluster_file="${CLUSTER_DIR}/cluster_N${n}_min.dat"
  local rundir="${OUTDIR}/N${n}/h0_${h0}/o${k}"
  [[ -f "${cluster_file}" ]] || { echo "  [SKIP] falta ${cluster_file}"; return; }
  mkdir -p "${rundir}"
  cp "${TABLE}" "${rundir}/table_repulsiva.txt"
  cd "${rundir}" && \
  lmp -in "${LAMMPS_SCRIPT}" \
      -var cluster_file "${cluster_file}" \
      -var h0 "${h0}" -var t_total "${T_TOTAL}" \
      -var ax "${ax}" -var ay "${ay}" -var az "${az}" \
      -log "${rundir}/lammps.log" -screen none \
      2>"${rundir}/stderr.log" \
    && echo "  OK N=${n} h0=${h0} o${k}" \
    || echo "  [ERROR] N=${n} h0=${h0} o${k} (ver stderr.log)"
}
export -f run_sim
export TABLE CLUSTER_DIR OUTDIR LAMMPS_SCRIPT T_TOTAL

echo "Sweep orientaciones: N={${N_VALUES}} h0={${H0_VALUES}} × ${#ANGLES[@]} orient."

PAIRS=()
for n in $N_VALUES; do
  for h0 in $H0_VALUES; do
    for k in "${!ANGLES[@]}"; do
      PAIRS+=("${n} ${h0} ${k} ${ANGLES[$k]}")
    done
  done
done

printf '%s\n' "${PAIRS[@]}" | xargs -P "${JOBS}" -I '{}' bash -c 'run_sim $@' _ {}

echo "Sweep de orientaciones completado: ${OUTDIR}"
