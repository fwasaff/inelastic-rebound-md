"""
leer_dump.py — lector en STREAMING de dumps de texto de LAMMPS.

Motivo: los dumps de las losas grandes (Nx=Ny=80, ~1.5e5 átomos, 360 frames)
pesan ~4 GB. Cargarlos con `f.readlines()` y guardar los frames en una lista
pedía >20 GB de RAM (cada línea es un objeto str de Python, con ~50 bytes de
overhead cada uno, más los arrays de todos los frames a la vez) — con 15 GB
físicos el sistema entra en swap thrashing y se congela.

Aquí nunca hay más de UN frame en memoria: se lee, se procesa, se descarta.
Medido sobre el dump de 4 GB de ondas_N79_esponja_v2: 154 MB de RSS, 74 s.

Convención de columnas: el dump se escribe en in.bounce_cluster como
    dump viz all custom N traj.lammpstrj id type x y z vx vy vz c_keratom c_pe_atom

Uso típico (dos pasadas: la 1 para el centro de impacto, la 2 para el análisis):

    import leer_dump as ld

    cx, cy, info = ld.centro_bola(traj)
    for fr in ld.iter_frames(traj, ("type", "x", "y", "ke")):
        surf = fr.data[fr.data[:, fr.col["type"]] == 1]
        ...

Leer el archivo dos veces del disco es mucho más barato que hundirse en swap.
"""

from collections import namedtuple
from itertools import islice

import numpy as np

# Índices de columna en el dump tal como lo escribe in.bounce_cluster
COL_DUMP = dict(id=0, type=1, x=2, y=3, z=4, vx=5, vy=6, vz=7, ke=8, pe=9)

# fr.col mapea nombre -> índice DENTRO de fr.data (ya reducido a las columnas
# pedidas), que en general no coincide con el índice en el archivo.
Frame = namedtuple("Frame", "step xlo xhi ylo yhi data col")


def iter_frames(path, cols=("type", "x", "y", "ke"), dtype=np.float32):
    """Genera un Frame por cada instantánea del dump, sin acumular nada.

    cols  : nombres de columna a extraer (ver COL_DUMP). Pedir solo las que se
            usan ahorra la mayor parte del tiempo de parseo: convertir texto a
            float es el costo dominante, y las columnas no pedidas ni se miran.
    dtype : float32 alcanza y sobra para posiciones/velocidades/energías que
            van a un gráfico, y ocupa la mitad que float64.
    """
    usecols = tuple(COL_DUMP[c] for c in cols)
    # np.loadtxt respeta el orden de usecols, así que data[:, k] es cols[k].
    col = {nombre: k for k, nombre in enumerate(cols)}

    with open(path) as f:
        while True:
            head = f.readline()
            if not head:
                return
            if not head.startswith("ITEM: TIMESTEP"):
                raise ValueError(f"{path}: se esperaba 'ITEM: TIMESTEP', hay {head[:40]!r}")
            step = int(f.readline())
            f.readline()                                   # ITEM: NUMBER OF ATOMS
            n = int(f.readline())
            f.readline()                                   # ITEM: BOX BOUNDS
            xlo, xhi = (float(v) for v in f.readline().split()[:2])
            ylo, yhi = (float(v) for v in f.readline().split()[:2])
            f.readline()                                   # cotas z
            f.readline()                                   # ITEM: ATOMS ...
            # islice(f, n) le entrega a loadtxt un iterador sobre las n líneas
            # del frame, sin materializarlas en una lista.
            data = np.loadtxt(islice(f, n), usecols=usecols, dtype=dtype)
            yield Frame(step, xlo, xhi, ylo, yhi, data, col)


InfoTraj = namedtuple("InfoTraj", "steps n_atoms xlo xhi ylo yhi Lx Ly z_surf0")


def centro_bola(path, verbose=True):
    """Pasada 1: centroide (x,y) de la bola promediado en el tiempo.

    Devuelve (cx, cy, info) con info.steps, info.Lx/Ly, info.z_surf0 (el z más
    alto de la superficie sin perturbar, del primer frame — el 'nivel del
    suelo' que usan los observables 3D).
    """
    if verbose:
        print(f"Leyendo {path} (streaming, pasada 1/2: centro de impacto) ...")
    steps, ball_xy = [], []
    info0 = None
    for fr in iter_frames(path, ("type", "x", "y", "z")):
        c = fr.col
        if info0 is None:
            es_surf = fr.data[:, c["type"]] == 1
            info0 = (fr.data.shape[0], fr.xlo, fr.xhi, fr.ylo, fr.yhi,
                     float(fr.data[es_surf, c["z"]].max()))
        steps.append(fr.step)
        es_bola = fr.data[:, c["type"]] == 2
        ball_xy.append(fr.data[es_bola][:, [c["x"], c["y"]]].mean(axis=0))

    if info0 is None:
        raise ValueError(f"{path}: no se leyó ningún frame")
    n_atoms, xlo, xhi, ylo, yhi, z_surf0 = info0
    info = InfoTraj(np.array(steps), n_atoms, xlo, xhi, ylo, yhi,
                    xhi - xlo, yhi - ylo, z_surf0)
    cx, cy = (float(v) for v in np.mean(ball_xy, axis=0))
    if verbose:
        print(f"  {len(info.steps)} frames, {n_atoms} átomos, "
              f"caja x=[{xlo:.2f},{xhi:.2f}] y=[{ylo:.2f},{yhi:.2f}]")
        print(f"  centro de impacto (x,y) = ({cx:.2f}, {cy:.2f}) sigma")
    return cx, cy, info


def media_por_bin(r, valores, r_edges):
    """Promedio de `valores` por capa radial, en una sola pasada (bincount).

    Reemplaza el patrón `[valores[(r>=e[k])&(r<e[k+1])].mean() for k in ...]`,
    que arma nbins máscaras booleanas de N elementos cada una: con N=1.5e5 y
    nbins=120 son 18 millones de comparaciones por frame. Los bins vacíos
    quedan en NaN (no en 0, que se confundiría con "medido y da cero").
    """
    nbins = len(r_edges) - 1
    idx = np.digitize(r, r_edges) - 1
    ok = (idx >= 0) & (idx < nbins)
    cnt = np.bincount(idx[ok], minlength=nbins)
    tot = np.bincount(idx[ok], weights=valores[ok], minlength=nbins)
    with np.errstate(invalid="ignore"):
        return np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan)


def t_impacto(bounce_dat):
    """Instante del impacto = mínimo de Zcm en bounce.dat."""
    d = np.loadtxt(bounce_dat, comments="#")
    return float(d[np.argmin(d[:, 1]), 0])
