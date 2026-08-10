#!/bin/bash
# ==============================================================================
# run_ondas_profunda.sh
# Corrida de losa PROFUNDA para la línea onda/calor.
#
# Motivo: con la losa de 4.65 sigma (Nz=6) el eco del fondo rígido vuelve a
# 1.8 tau, y el forzamiento (contacto bola-superficie) dura 3.2 tau. La ventana
# limpia es más corta que el forzamiento, así que no se puede medir ninguna
# relajación. Con Nz=40 (losa de 62 sigma) el eco se va a ~14 tau.
#
# Ventanas limpias resultantes, ya balanceadas:
#   fondo   2*58.9/8.28 = 14.2 tau   (antes: 1.8)
#   lateral 62.0/4.6    = 13.5 tau
# No conviene ir más profundo: a partir de ahí el lateral pasa a ser el límite.
#
# Uso:
#   bash scripts/run_ondas_profunda.sh                  # h0=28.15, en 2º plano
#   bash scripts/run_ondas_profunda.sh --h0 79.03       # otra velocidad
#   bash scripts/run_ondas_profunda.sh --fg             # en 1er plano (bloquea)
#   bash scripts/run_ondas_profunda.sh --np 4           # menos procesos
#   bash scripts/run_ondas_profunda.sh --recorte        # dump solo capas de
#                                                       # arriba (~2.7 GB vs 5)
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
LAMMPS_IN="$ROOT_DIR/lammps/in.bounce_cluster"
CLUSTER="$ROOT_DIR/lammps/clusters/cluster_N79_min.dat"
TABLE_SRC="$ROOT_DIR/data/outputs/activo/ondas_N79_libre/table_repulsiva.txt"
LAMMPS="${LAMMPS_EXE:-lmp}"

H0=28.15
NZ=40

# Cores FÍSICOS, no hilos lógicos. OpenMPI asigna un slot por core físico y
# aborta con "not enough slots available" si se le piden más -- que es lo que
# pasaba usando nproc (en un i7-8550U nproc dice 8, pero son 4 cores con
# hyperthreading). Y para MD, que es cómputo en coma flotante saturado, los
# hilos extra del HT no aportan casi nada: compiten por la misma unidad
# vectorial.
NP="$(lscpu -p=CORE 2>/dev/null | grep -v '^#' | sort -u | wc -l)"
[[ -z "$NP" || "$NP" -lt 1 ]] && NP="$(grep '^core id' /proc/cpuinfo 2>/dev/null | sort -u | wc -l)"
[[ -z "$NP" || "$NP" -lt 1 ]] && NP=1
FG=0
RECORTE=0
NOMBRE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --h0)      H0="$2"; shift 2 ;;
        --np)      NP="$2"; shift 2 ;;
        --nz)      NZ="$2"; shift 2 ;;
        --nombre)  NOMBRE="$2"; shift 2 ;;
        --fg)      FG=1; shift ;;
        --recorte) RECORTE=1; shift ;;
        *) echo "Opción desconocida: $1"; exit 1 ;;
    esac
done

[[ -z "$NOMBRE" ]] && NOMBRE="ondas_N79_profunda_h${H0}"
OUTDIR="$ROOT_DIR/data/outputs/activo/$NOMBRE"

# --- Parámetros derivados -----------------------------------------------------
# La bola arranca en z = Nz*a0 + h0 y contacta en ~Nz*a0 + 1.73, así que cae
# (h0 - 1.73) bajo g. De ahí t_imp; el dump arranca 2 tau antes, y la corrida
# termina 16 tau después del impacto (la ventana limpia es ~13.5).
#
# La fórmula SUBESTIMA t_imp en ~1.5% (medido: predice 45.1, real 45.8), porque
# el mínimo de Zcm llega algo después del primer contacto: la bola sigue bajando
# mientras se deforma. El error va en la dirección segura -- el dump arranca
# antes de lo pedido (2.7 tau en vez de 2), nunca después.
read -r T_IMP T_TOTAL DUMP_DELAY <<< "$(python3 -c "
h0=$H0; g=0.026; dt=0.005
t_imp=(2*(h0-1.73)/g)**0.5
print(f'{t_imp:.1f} {t_imp+16:.0f} {round((t_imp-2)/dt)}')
")"

EXTRA=""
if [[ "$RECORTE" == "1" ]]; then
    # 27 sigma de arriba en una losa de 62: contiene de sobra el modo de
    # Rayleigh, que decae en ~1 lambda (~17 sigma).
    EXTRA="-var z_dump_bot 35.0"
fi

echo "=============================================================="
echo " Corrida profunda: $NOMBRE"
echo "   Nz=$NZ  ->  losa de $(python3 -c "print(f'{$NZ*1.5496:.1f}')") sigma"
echo "   h0=$H0  ->  t_imp = $T_IMP tau,  t_total = $T_TOTAL tau"
echo "   dump desde el paso $DUMP_DELAY (2 tau antes del impacto)"
echo "   $NP procesos MPI"
[[ "$RECORTE" == "1" ]] && echo "   dump RECORTADO a z > 35 sigma"
echo "   salida: $OUTDIR"
echo "=============================================================="

mkdir -p "$OUTDIR"
cp "$TABLE_SRC" "$OUTDIR/table_repulsiva.txt"
cd "$OUTDIR"

CMD=(mpirun -np "$NP" "$LAMMPS" -in "$LAMMPS_IN"
     -var cluster_file "$CLUSTER"
     -var bc_xy s
     -var Nx 80 -var Ny 80 -var Nz "$NZ"
     -var h0 "$H0"
     -var t_total "$T_TOTAL"
     -var dump_freq 50
     -var dump_delay "$DUMP_DELAY"
     $EXTRA
     -log lammps.log -screen none)

if [[ "$FG" == "1" ]]; then
    "${CMD[@]}"
    echo "Terminado: $OUTDIR"
else
    # nohup + & : la corrida sobrevive a que cierres la terminal.
    nohup "${CMD[@]}" > corrida.out 2>&1 &
    PID=$!
    echo "$PID" > corrida.pid
    echo
    echo "Corriendo en SEGUNDO PLANO (PID $PID). La terminal queda libre."
    echo
    echo "  ver avance:     tail -f $OUTDIR/lammps.log"
    echo "  ¿sigue viva?:   ps -p $PID"
    echo "  tamaño del dump: du -h $OUTDIR/traj.lammpstrj"
    echo "  abortar:        kill $PID"
    echo
    echo "Al terminar, verifica que la ventana se abrió:"
    echo "  cd $ROOT_DIR && python3 scripts/medir_eco_fondo.py \\"
    echo "      --traj $OUTDIR/traj.lammpstrj \\"
    echo "      --bounce $OUTDIR/bounce.dat \\"
    echo "      --out data/outputs/activo/fig_eco_${NOMBRE}.png --tmax 16"
    echo "  (el eco debe salir cerca de 14 tau, no de 1.8)"
fi
