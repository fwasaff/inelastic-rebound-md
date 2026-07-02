"""
estilo_figuras.py — Estilo matplotlib común para las figuras del paper.

Convenciones de revista:
  - ancho de figura = \textwidth del artículo (390 pt = 5.42 in), escala 1:1
  - fuentes 8-9 pt al tamaño final, serif + mathtext CM (coherente con LaTeX)
  - sin títulos sobre los paneles (esa información va en el caption)
  - etiquetas de panel (a), (b), ... en negrita dentro del eje
  - ticks hacia adentro en los cuatro lados
"""

import numpy as np
from pathlib import Path
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

G = 0.026
COLS = "time Zcm Vcmz KB TB VB Uc UF KB3d KS PE_surf TS".split()

ANCHO = 5.42          # in = 390 pt = \textwidth (article 12pt a4)

# Colores consistentes por N en todas las figuras
COLOR_N = {1: "0.55", 2: "#8c564b", 7: "#2ca02c", 13: "#1f77b4",
           55: "#ff7f0e", 79: "#d62728", 87: "#9467bd"}
MARKER_N = {1: "x", 2: "+", 7: "v", 13: "s", 55: "^", 79: "o", 87: "D"}


def aplicar_estilo():
    plt.rcParams.update({
        "font.size": 9,
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.5,
        "lines.linewidth": 1.2,
        "lines.markersize": 4.5,
        "axes.linewidth": 0.7,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "savefig.bbox": "tight",
        "legend.framealpha": 0.9,
        "legend.edgecolor": "0.8",
        "legend.handlelength": 1.6,
    })


def panel_label(ax, s, x=0.03, y=0.96):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", zorder=10)


def cargar_run(path, M):
    """Partición del primer rebote de un bounce.dat (M = masa de la bola).

    Misma definición que analizar_reprod_fisica.cargar_rebotes:
      a_CM = h1/h0_eff, a_surf = Δ(KS+PE_surf)/E0 (directo),
      a_ball = 1 − a_CM − a_surf, a_TB = TB(ap1)/E0, a_VB = a_ball − a_TB.
    """
    d = np.loadtxt(Path(path), comments="#")
    arr = {c: d[:, i] for i, c in enumerate(COLS[:d.shape[1]])}
    Zcm, TB = arr["Zcm"], arr["TB"]
    KS = arr.get("KS", np.zeros_like(Zcm))
    PEs = arr.get("PE_surf", np.zeros_like(Zcm))

    mins, _ = find_peaks(-Zcm, prominence=0.5, distance=50)
    maxs, _ = find_peaks(Zcm, prominence=0.5, distance=50)
    if len(mins) == 0 or len(maxs) == 0:
        return None

    z0 = float(Zcm[mins].mean())
    h0_eff = float(Zcm[0] - z0)
    E0 = M * G * h0_eff
    m1 = maxs[0]

    a_CM = float((Zcm[m1] - z0) / h0_eff)
    a_surf = float((KS[m1] + PEs[m1] - KS[0] - PEs[0]) / E0)
    a_ball = 1.0 - a_CM - a_surf
    a_TB = float(TB[m1] / E0)

    return {
        "arr": arr, "mins": mins, "maxs": maxs, "z0": z0,
        "h0_eff": h0_eff, "E0": E0,
        "v_imp": float(np.sqrt(2 * G * h0_eff)),
        "e1": float(np.sqrt(max(a_CM, 0.0))),
        "a_CM": a_CM, "a_TB": a_TB, "a_surf": a_surf,
        "a_ball": a_ball, "a_VB": a_ball - a_TB,
    }
