"""
gen_cluster_init.py  —  Genera geometría inicial de un nanocúmulo FCC LJ

Crea un archivo LAMMPS data file con N átomos dispuestos en una red FCC
perfecta, centrados en el origen, listos para minimización posterior con
in.minimize_cluster.

Uso:
    python3 scripts/gen_cluster_init.py --N 79
    python3 scripts/gen_cluster_init.py --N 55
    python3 scripts/gen_cluster_init.py --N 87

Shells FCC (centrado en sitio de red, a0=1.5496σ):
    N=13 (r<1.2σ), N=19, N=43, N=55, N=79, N=87, N=135, ...
"""

import argparse
import numpy as np
from pathlib import Path

A0 = 1.5496   # parámetro de red FCC argón LJ [σ]
PAD = 6.0     # padding de caja alrededor del cluster [σ]

# Shells FCC completos (átomos por shell, distancia al centro):
# Solo shells exactos son reproducibles y físicamente limpios.
COMPLETE_SHELLS = {1, 13, 19, 43, 55, 79, 87, 135}


def fcc_lattice_sites(n_cells, a0=A0):
    """Devuelve todas las posiciones FCC en una grilla de ±n_cells celdas."""
    basis = a0 * np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ])
    pts = []
    for nx in range(-n_cells, n_cells + 1):
        for ny in range(-n_cells, n_cells + 1):
            for nz in range(-n_cells, n_cells + 1):
                for b in basis:
                    pts.append(np.array([nx, ny, nz]) * a0 + b)
    return np.array(pts)


def select_cluster(N, a0=A0):
    """
    Selecciona los N átomos más cercanos al origen en la red FCC.
    Para N = shell completo, todos los átomos del último shell tienen
    la misma distancia → resultado único y reproducible.
    """
    n_cells = int(3.5 * a0) + 3
    all_pts  = fcc_lattice_sites(n_cells, a0)
    dists    = np.linalg.norm(all_pts, axis=1)
    idx_sort = np.argsort(dists)

    # Verificar si N es un número de shell completo
    # (todos los átomos del último shell están incluidos)
    d_Nth    = dists[idx_sort[N - 1]]
    d_next   = dists[idx_sort[N]]     # primero excluido
    if abs(d_Nth - d_next) < 1e-6:
        print(f"ADVERTENCIA: N={N} cae en medio de un shell degenerado.")
        print(f"  d[{N}-1]={d_Nth:.4f}σ  d[{N}]={d_next:.4f}σ")
        print(f"  Shells completos recomendados: {sorted(COMPLETE_SHELLS)}")
        print(f"  Procediendo igualmente con los {N} más cercanos al origen.")
    else:
        R_cluster = d_Nth
        print(f"N={N}: todos en shells completos. Radio exterior = {R_cluster:.4f}σ")

    pos = all_pts[idx_sort[:N]]
    # Centrar exactamente en el origen (debería serlo, pero por simetría numérica)
    pos -= pos.mean(axis=0)
    return pos


def write_lammps_data(pos, outfile, comment=""):
    """Escribe un LAMMPS data file en formato atomic con átomos tipo 2."""
    N   = len(pos)
    ext = pos.max(axis=0) - pos.min(axis=0)
    L   = ext.max() / 2 + PAD   # semilado de caja cúbica con padding

    lines = []
    lines.append(f"# FCC LJ cluster  N={N}  a0={A0}  {comment}\n")
    lines.append(f"\n{N} atoms\n2 atom types\n")
    lines.append(f"\n{-L:.8f} {L:.8f} xlo xhi")
    lines.append(f"\n{-L:.8f} {L:.8f} ylo yhi")
    lines.append(f"\n{-L:.8f} {L:.8f} zlo zhi\n")
    lines.append("\nMasses\n\n1 1.0\n2 1.0\n")
    lines.append("\nAtoms # atomic\n\n")
    for i, (x, y, z) in enumerate(pos, start=1):
        lines.append(f"{i} 2 {x:.16e} {y:.16e} {z:.16e} 0 0 0\n")

    Path(outfile).parent.mkdir(parents=True, exist_ok=True)
    with open(outfile, 'w') as f:
        f.writelines(lines)
    print(f"Escrito: {outfile}  ({N} átomos)")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--N', type=int, default=79,
                        help='Número de átomos del cluster (default: 79)')
    parser.add_argument('--a0', type=float, default=A0,
                        help=f'Parámetro de red FCC [σ] (default: {A0})')
    parser.add_argument('--outdir', default='lammps/clusters',
                        help='Directorio de salida')
    args = parser.parse_args()

    print(f"\nGenerando cluster FCC  N={args.N}  a0={args.a0}σ")
    pos = select_cluster(args.N, args.a0)

    R_max = np.linalg.norm(pos, axis=1).max()
    PE_approx = 0.0  # solo geometría, sin potencial calculado aquí
    print(f"  Radio máximo:  {R_max:.4f}σ")
    print(f"  CM desplazado: {pos.mean(axis=0)}")

    outfile = f"{args.outdir}/cluster_N{args.N}_init.dat"
    write_lammps_data(pos, outfile,
                      comment=f"generado por gen_cluster_init.py  a0={args.a0}")


if __name__ == '__main__':
    main()
