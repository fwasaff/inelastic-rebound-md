#!/usr/bin/env bash
# run_cluster_sweep.sh — corre in.bounce_cluster para múltiples N y alturas.
#
# Uso:
#   bash scripts/run_cluster_sweep.sh [--jobs N_PARALLEL]
#   bash scripts/run_cluster_sweep.sh --n 1           # solo N=1
#   bash scripts/run_cluster_sweep.sh --n "1 2 7"     # varios N
#
# Los resultados se guardan en:
#   data/outputs/cluster_sweep/N{N}/h0_{h0}/bounce.dat

set -euo pipefail

# Valores por defecto
N_VALUES="1"
H0_VALUES="11.03 18.64 28.15 41.47 58.6 79.53"
JOBS=4
SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
OUTDIR="${SCRIPT_DIR}/data/outputs/cluster_sweep"
LAMMPS_SCRIPT="${SCRIPT_DIR}/lammps/in.bounce_cluster"
CLUSTER_DIR="${SCRIPT_DIR}/lammps/clusters"
TABLE="${SCRIPT_DIR}/data/outputs/reprod_2012/table_repulsiva.txt"

# --- Parse argumentos ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --n)     N_VALUES="$2";  shift 2 ;;
        --h0)    H0_VALUES="$2"; shift 2 ;;
        --jobs)  JOBS="$2";      shift 2 ;;
        *)       echo "Argumento desconocido: $1"; exit 1 ;;
    esac
done

echo "=================================================="
echo "Sweep de nanocúmulos LJ"
echo "N: ${N_VALUES}"
echo "h0: ${H0_VALUES}"
echo "Salida: ${OUTDIR}"
echo "=================================================="

# Verificar tabla repulsiva
if [[ ! -f "${TABLE}" ]]; then
    echo "[1/2] Generando tabla de potencial repulsivo..."
    mkdir -p "$(dirname "${TABLE}")"
    python3 "${SCRIPT_DIR}/scripts/gen_tabla_repulsiva.py" --output "${TABLE}"
else
    echo "[1/2] Tabla ya existe: ${TABLE}"
fi

# --- Función para correr una simulación individual ---------------------------
run_sim() {
    local n="$1"
    local h0="$2"
    local cluster_file="${CLUSTER_DIR}/cluster_N${n}_min.dat"
    local rundir="${OUTDIR}/N${n}/h0_${h0}"

    if [[ ! -f "${cluster_file}" ]]; then
        echo "  [SKIP] Cluster no encontrado: ${cluster_file}"
        return
    fi

    mkdir -p "${rundir}"
    cp "${TABLE}" "${rundir}/table_repulsiva.txt"

    # Correr desde el directorio de salida para que bounce.dat quede ahí
    # Usar paths absolutos para que LAMMPS encuentre los archivos
    cd "${rundir}" && \
    lmp -in "${LAMMPS_SCRIPT}" \
        -var cluster_file "${cluster_file}" \
        -var h0 "${h0}" \
        -log "${rundir}/lammps.log" \
        -screen none \
        2>"${rundir}/stderr.log" || {
        echo "  [ERROR] N=${n}, h0=${h0} falló. Ver: ${rundir}/stderr.log"
        return
    }

    echo "  Completado: N=${n}, h0=${h0}"
}

export -f run_sim
export TABLE CLUSTER_DIR OUTDIR LAMMPS_SCRIPT

# --- Lanzar simulaciones en paralelo -----------------------------------------
echo "[2/2] Corriendo simulaciones (${JOBS} jobs en paralelo)..."

# Construir lista de pares (N, h0)
PAIRS=()
for n in $N_VALUES; do
    for h0 in $H0_VALUES; do
        PAIRS+=("${n} ${h0}")
    done
done

# Correr en paralelo
printf '%s\n' "${PAIRS[@]}" | \
    xargs -P "${JOBS}" -I '{}' bash -c 'run_sim $@' _ {}

echo ""
echo "=================================================="
echo "Sweep completado."
echo "Analizar con:"
echo "  python3 scripts/analizar_cluster_sweep.py"
echo "=================================================="
