#!/usr/bin/env bash
# run_slab_test.sh — test de convergencia de tamaño finito de losa (punto M4).
# Corre N=79 a h0 intermedia con la losa estándar (8x8x6) y una losa el doble
# de grande en las tres direcciones (16x16x12), con idénticos parámetros, para
# ver si alpha_surf/alpha_ball dependen del tamaño de la losa (reflexión de
# fonones en el fondo fijo / bordes periódicos).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
OUTDIR="${SCRIPT_DIR}/data/outputs/terminado/slab_test"
LAMMPS_SCRIPT="${SCRIPT_DIR}/lammps/in.bounce_cluster"
CLUSTER="${SCRIPT_DIR}/lammps/clusters/cluster_N79_min.dat"
TABLE="${SCRIPT_DIR}/data/outputs/terminado/reprod_2012/table_repulsiva.txt"
H0=41.47

run_slab() {
  local tag="$1" nx="$2" ny="$3" nz="$4"
  local rundir="${OUTDIR}/${tag}"
  mkdir -p "${rundir}"
  cp "${TABLE}" "${rundir}/table_repulsiva.txt"
  cd "${rundir}" && \
  lmp -in "${LAMMPS_SCRIPT}" \
      -var cluster_file "${CLUSTER}" -var h0 "${H0}" -var t_total 200.0 \
      -var Nx "${nx}" -var Ny "${ny}" -var Nz "${nz}" \
      -log "${rundir}/lammps.log" -screen none 2>"${rundir}/stderr.log" \
    && echo "  OK ${tag} (${nx}x${ny}x${nz})" || { echo "  [ERROR] ${tag}"; tail -5 "${rundir}/stderr.log"; }
}

run_slab "slab_8x8x6"    8  8  6  &
run_slab "slab_16x16x12" 16 16 12 &
wait
echo "Test de losa listo en ${OUTDIR}"
