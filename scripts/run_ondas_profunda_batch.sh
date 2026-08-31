#!/bin/bash
# ==============================================================================
# run_ondas_profunda_batch.sh
#
# Corre las 3 alturas ya establecidas de la línea onda/calor (h28.15, h56.93,
# h79.03) EN SECUENCIA (no en paralelo -- son 4 procesos MPI cada una, en un
# CPU de 4 cores físicos; correrlas juntas solo las haría competir por el
# mismo cómputo) con el dump SIN recortar (para poder medir el eco de fondo
# real, ver la nota 2026-08-12 en run_ondas_profunda.sh) y --margen 20.
#
# Motivo (2026-08-12): la bitácora (docs/bitacora/00_pregunta_activa.md) pide
# el discriminador tau_th(v) -- tres velocidades, no una. Las tres corridas
# "profunda" que ya existen usan --recorte (z_dump_bot=35), que no deja ver
# si el eco de fondo real (~14.2 tau predicho) está o no dentro de la ventana.
# Este batch las relanza sin ese recorte, con nombre "_zfull" para no pisar
# las corridas viejas (que igual sirven para lo lateral, ya medido).
#
# Estimado (medido en la corrida real de h79.03, 2.06 pasos/s): ~7.6 h para
# las 3 juntas -- pensado para dejarlo corriendo toda la noche.
#
# IMPORTANTE -- que la laptop NO se suspenda durante la corrida:
#   systemd-inhibit --what=idle:sleep:handle-lid-switch \
#       --why="corrida LAMMPS overnight" \
#       bash scripts/run_ondas_profunda_batch.sh
# Así, aunque cierres la tapa, el sistema no se suspende mientras el batch
# siga vivo (systemd-inhibit libera el bloqueo solo al terminar el comando).
#
# Uso:
#   bash scripts/run_ondas_profunda_batch.sh                # margen 20, sufijo _zfull
#   bash scripts/run_ondas_profunda_batch.sh --margen 25
#   bash scripts/run_ondas_profunda_batch.sh --sufijo _test
#   bash scripts/run_ondas_profunda_batch.sh --alturas "39.42 67.55"   # otras alturas
#
# 2026-08-13: el sistema es DETERMINISTA (in.bounce_cluster hace
# `velocity all set 0.0 0.0 0.0`, sin ruido térmico inicial) -- repetir la
# misma altura da bit a bit la misma trayectoria. Para más puntos, siempre
# alturas NUEVAS, nunca repetir una ya corrida esperando otra semilla.
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

MARGEN=20
SUFIJO="_zfull"
ALTURAS=(28.15 56.93 79.03)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --margen)  MARGEN="$2"; shift 2 ;;
        --sufijo)  SUFIJO="$2"; shift 2 ;;
        --alturas) read -r -a ALTURAS <<< "$2"; shift 2 ;;
        *) echo "Opción desconocida: $1"; exit 1 ;;
    esac
done

INICIO_BATCH=$(date +%s)
echo "=============================================================="
echo " BATCH ondas_N79_profunda -- ${#ALTURAS[@]} alturas, margen=$MARGEN, sufijo=$SUFIJO"
echo " Empezó: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=============================================================="

for H0 in "${ALTURAS[@]}"; do
    NOMBRE="ondas_N79_profunda_h${H0}${SUFIJO}"
    echo
    echo ">>> [$( date '+%H:%M:%S' )] Lanzando h0=$H0 -> $NOMBRE"
    t0=$(date +%s)
    bash "$SCRIPT_DIR/run_ondas_profunda.sh" --h0 "$H0" --margen "$MARGEN" \
         --nombre "$NOMBRE" --fg
    t1=$(date +%s)
    echo ">>> [$( date '+%H:%M:%S' )] Terminó h0=$H0 en $(( (t1-t0)/60 )) min"
done

FIN_BATCH=$(date +%s)
echo
echo "=============================================================="
echo " BATCH TERMINADO: $(date '+%Y-%m-%d %H:%M:%S')"
echo " Duración total: $(( (FIN_BATCH-INICIO_BATCH)/3600 )) h $(( ((FIN_BATCH-INICIO_BATCH)%3600)/60 )) min"
echo "=============================================================="
echo
echo "Siguiente paso para cada corrida (verificar la ventana antes de confiar"
echo "en cualquier medición de frente/estela):"
for H0 in "${ALTURAS[@]}"; do
    NOMBRE="ondas_N79_profunda_h${H0}${SUFIJO}"
    OUTDIR="data/outputs/activo/$NOMBRE"
    echo "  python3 scripts/medir_eco_fondo.py --traj $OUTDIR/traj.lammpstrj \\"
    echo "      --bounce $OUTDIR/bounce.dat --out $OUTDIR/fig_eco_fondo.png --tmax $MARGEN"
done
