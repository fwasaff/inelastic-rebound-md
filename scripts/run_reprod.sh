#!/bin/bash
# ==============================================================================
# run_reprod.sh
# Corre la reproducción de Gutiérrez et al. 2012 para las 6 alturas del paper.
#
# Uso: bash scripts/run_reprod.sh [--jobs N]
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
OUTDIR="$ROOT_DIR/data/outputs/reprod_2012"
LAMMPS_IN="$ROOT_DIR/lammps/in.bounce_reprod"
TABLE="$OUTDIR/table_repulsiva.txt"
LAMMPS="${LAMMPS_EXE:-lmp}"
JOBS=1

# Parsear argumentos
while [[ $# -gt 0 ]]; do
    case "$1" in
        --jobs) JOBS="$2"; shift 2 ;;
        *) echo "Uso: $0 [--jobs N]"; exit 1 ;;
    esac
done

# Alturas del paper original (Gutiérrez et al. 2012, Fig. 7)
H0_VALUES=(11.03 18.64 28.15 41.47 58.6 79.53)

mkdir -p "$OUTDIR"

echo "=================================================="
echo "Reproducción Gutiérrez et al. 2012"
echo "Alturas: ${H0_VALUES[*]}"
echo "Salida: $OUTDIR"
echo "=================================================="

# --- Paso 1: Generar tabla de potencial repulsivo ----------------------------
echo ""
echo "[1/2] Generando tabla de potencial LJ repulsivo (c=0)..."
python3 "$SCRIPT_DIR/gen_tabla_repulsiva.py" --output "$TABLE"

echo "Tabla guardada en: $TABLE"

# --- Paso 2: Correr simulaciones ---------------------------------------------
echo ""
echo "[2/2] Corriendo simulaciones (${#H0_VALUES[@]} alturas, $JOBS jobs en paralelo)..."

run_one() {
    local h0="$1"
    local rundir="$OUTDIR/h0_${h0}"
    mkdir -p "$rundir"

    # Copiar tabla al directorio de trabajo (LAMMPS la busca localmente)
    cp "$TABLE" "$rundir/"

    echo "  Iniciando h0=${h0} sigma → $rundir"
    cd "$rundir"
    $LAMMPS -in "$LAMMPS_IN" -var h0 "$h0" \
            -log "lammps_h0_${h0}.log" \
            -screen none
    echo "  Completado: h0=${h0}"
}

export -f run_one
export OUTDIR TABLE LAMMPS_IN LAMMPS

if command -v parallel &>/dev/null && [[ $JOBS -gt 1 ]]; then
    printf '%s\n' "${H0_VALUES[@]}" | parallel -j "$JOBS" run_one {}
else
    for h0 in "${H0_VALUES[@]}"; do
        run_one "$h0"
    done
fi

echo ""
echo "=================================================="
echo "Todas las simulaciones completadas."
echo "Analizar con:"
echo "  python3 scripts/analizar_reprod.py --outdir $OUTDIR"
echo "=================================================="
