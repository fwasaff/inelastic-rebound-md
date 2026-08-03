"""
medir_Nc.py — Mide directamente el número efectivo de átomos de superficie
activados en el contacto, N_c, a partir de un dump por-átomo de LAMMPS.

Responde al punto M2 del arbitraje: en vez de dejar N_c como parámetro de
ajuste (con tres calibraciones distintas), se mide en la simulación.

Dos estimadores independientes, evaluados en la ventana de contacto:

  1. N_c^part  — número de participación de la energía cinética de superficie
                 (inverse participation ratio, sin umbrales libres):
                     N_c^part = (Σ_i K_i)^2 / Σ_i K_i^2
                 sobre los átomos de superficie. Cuenta cuántos átomos
                 comparten efectivamente la energía transferida.

  2. N_c^geom  — huella geométrica: nº de átomos de superficie a distancia
                 < r_c = 2.5 σ de algún átomo de la bola en la compresión
                 máxima (átomos que "sienten" a la bola).

Uso:
    python3 scripts/medir_Nc.py <dump.lammpstrj> [--rc 2.5]
    python3 scripts/medir_Nc.py --sweep     # todos los N de data/outputs/Nc_runs
"""
import sys
import numpy as np
from pathlib import Path

RC = 2.5
DT = 0.005


def iter_frames(path):
    """Generador: (timestep, ids, types, xyz(N,3), ke(N)) por frame."""
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                return
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            ts = int(f.readline())
            f.readline()                       # ITEM: NUMBER OF ATOMS
            nat = int(f.readline())
            f.readline()                       # ITEM: BOX BOUNDS
            for _ in range(3):
                f.readline()
            cols = f.readline().split()[2:]    # ITEM: ATOMS id type x y z ...
            ix = {c: i for i, c in enumerate(cols)}
            data = np.empty((nat, len(cols)))
            for i in range(nat):
                data[i] = np.fromstring(f.readline(), sep=" ")
            yield ts, data, ix


def analizar_dump(path, rc=RC, verbose=True):
    """Recorre el dump, localiza la compresión máxima y devuelve N_c."""
    zcm_ball = []
    ksurf_tot = []
    npart = []
    ngeom = []
    tss = []
    for ts, d, ix in iter_frames(path):
        typ = d[:, ix["type"]].astype(int)
        z = d[:, ix["z"]]
        xyz = d[:, [ix["x"], ix["y"], ix["z"]]]
        ke = d[:, ix["c_keratom"]]
        is_ball = typ == 2
        is_surf = typ == 1
        tss.append(ts)
        zcm_ball.append(z[is_ball].mean())

        ks = ke[is_surf]
        ktot = ks.sum()
        ksurf_tot.append(ktot)
        # número de participación
        npart.append((ktot ** 2 / np.sum(ks ** 2)) if ktot > 0 else 0.0)

        # huella geométrica: átomos de superficie con algún átomo de bola a < rc
        surf_xyz = xyz[is_surf]
        ball_xyz = xyz[is_ball]
        # solo superficie en la zona superior (z > mitad de la caja de losa)
        # para no contar por accidente; el criterio de distancia ya lo asegura
        d2min = np.full(surf_xyz.shape[0], np.inf)
        for b in ball_xyz:
            dd = np.sum((surf_xyz - b) ** 2, axis=1)
            d2min = np.minimum(d2min, dd)
        ngeom.append(int(np.sum(d2min < rc ** 2)))

    tss = np.array(tss); zcm_ball = np.array(zcm_ball)
    ksurf_tot = np.array(ksurf_tot); npart = np.array(npart); ngeom = np.array(ngeom)

    i_comp = int(np.argmin(zcm_ball))          # compresión máxima
    i_kmax = int(np.argmax(ksurf_tot))         # máxima KE de superficie

    res = {
        "t_comp": tss[i_comp] * DT,
        "zcm_min": float(zcm_ball[i_comp]),
        "Nc_part_comp": float(npart[i_comp]),
        "Nc_geom_comp": int(ngeom[i_comp]),
        "Nc_part_kmax": float(npart[i_kmax]),
        "Nc_geom_kmax": int(ngeom[i_kmax]),
        # promedio de la participación en la ventana de contacto (KE > 5% del pico)
        "Nc_part_win": float(np.mean(npart[ksurf_tot > 0.05 * ksurf_tot.max()])),
    }
    if verbose:
        print(f"  compresión máx en t={res['t_comp']:.1f} τ (Zcm_min={res['zcm_min']:.2f})")
        print(f"  N_c^part (compresión) = {res['Nc_part_comp']:.1f}")
        print(f"  N_c^part (KE_surf máx) = {res['Nc_part_kmax']:.1f}")
        print(f"  N_c^part (ventana)     = {res['Nc_part_win']:.1f}")
        print(f"  N_c^geom (compresión)  = {res['Nc_geom_comp']}")
    return res


def sweep(base):
    """Analiza data/outputs/Nc_runs/N{N}/traj.lammpstrj para varios N."""
    base = Path(base)
    filas = []
    for d in sorted(base.glob("N*"), key=lambda p: int(p.name[1:])):
        N = int(d.name[1:])
        dump = d / "traj.lammpstrj"
        if not dump.exists():
            continue
        print(f"N={N}:")
        r = analizar_dump(dump)
        r["N"] = N
        filas.append(r)
        print()

    print("=" * 74)
    print(f"{'N':>4} {'Nc_part(comp)':>14} {'Nc_part(win)':>13} {'Nc_geom':>9} {'Nc/N':>7}")
    print("-" * 74)
    for r in filas:
        print(f"{r['N']:>4} {r['Nc_part_comp']:>14.1f} {r['Nc_part_win']:>13.1f} "
              f"{r['Nc_geom_comp']:>9} {r['Nc_part_win']/r['N']:>7.2f}")

    # Ajuste de ley de potencias N_c = A N^beta (estimador de ventana)
    Ns = np.array([r["N"] for r in filas], float)
    Nc = np.array([r["Nc_part_win"] for r in filas], float)
    m = Ns >= 7                                # excluir N=1,2 (sin huella real)
    if m.sum() >= 2:
        beta, logA = np.polyfit(np.log(Ns[m]), np.log(Nc[m]), 1)
        A = np.exp(logA)
        print("-" * 74)
        print(f"Ajuste (N>=7):  N_c = {A:.1f} · N^{beta:.3f}")
    np.save(base / "Nc_medido.npy", filas, allow_pickle=True)
    return filas


if __name__ == "__main__":
    if "--sweep" in sys.argv:
        base = "data/outputs/Nc_runs"
        sweep(base)
    else:
        path = sys.argv[1]
        rc = RC
        if "--rc" in sys.argv:
            rc = float(sys.argv[sys.argv.index("--rc") + 1])
        analizar_dump(path, rc=rc)
