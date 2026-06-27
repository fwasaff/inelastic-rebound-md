"""
Genera la tabla de potencial puramente repulsivo para la interacción
bola-superficie en la reproducción de Gutiérrez et al. 2012.

V(r) = 4*eps*(sigma/r)^12   (término c=0, sin atracción)
F(r) = -dV/dr = 48*eps*sigma^12 / r^13

En unidades LJ: eps=1, sigma=1.
"""

import numpy as np
from pathlib import Path

def generar_tabla(r_min=0.75, r_max=2.5, N=5000,
                  output="table_repulsiva.txt"):
    r = np.linspace(r_min, r_max, N)
    V = 4.0 / r**12
    F = 48.0 / r**13   # F = -dV/dr, positivo = repulsivo (convenio LAMMPS)

    out = Path(output)
    with out.open("w") as f:
        f.write("# Potencial LJ puramente repulsivo (c=0)\n")
        f.write("# V(r) = 4*(1/r)^12  [unidades LJ]\n")
        f.write("# Ref: Gutiérrez et al. 2012, interacción bola-superficie\n")
        f.write("#\n")
        f.write("REPULSIVE_LJ\n")
        f.write(f"N {N} R {r_min:.4f} {r_max:.4f}\n\n")
        for i, (ri, Vi, Fi) in enumerate(zip(r, V, F), start=1):
            f.write(f"{i:6d}  {ri:.8f}  {Vi:.10f}  {Fi:.10f}\n")

    print(f"Tabla generada: {out.resolve()} ({N} puntos, r ∈ [{r_min}, {r_max}] σ)")

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--output", default="table_repulsiva.txt")
    args = p.parse_args()
    generar_tabla(output=args.output)
