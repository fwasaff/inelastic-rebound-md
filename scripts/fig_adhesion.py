"""
fig_adhesion.py — Figura del barrido con adhesión (c_ps=1) para el paper.

Genera data/outputs/fig_paso_C_adhesion.pdf:
  (a) Z_cm(t) para N=79, 6 alturas: la bola pega en todos los casos
  (b) energía de adhesión E_adh(N) y umbral extrapolado h0*(N)

Uso:  python3 scripts/fig_adhesion.py
"""

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, str(Path(__file__).parent))
from estilo_figuras import aplicar_estilo, panel_label, ANCHO, COLS, G

LANG = "en" if "--en" in sys.argv else "es"
T = {"es": {"tope": "tope de la superficie",
            "umbral": r"$h_0^*$ (ec. umbral)",
            "ylabel_b": r"$h_0^*$ estimado $(\sigma)$"},
     "en": {"tope": "surface top",
            "umbral": r"$h_0^*$ (threshold eq.)",
            "ylabel_b": r"estimated $h_0^*$ $(\sigma)$"}}[LANG]

BASE = Path("data/outputs")
NS = [13, 55, 79, 87]
A_CM = 0.35          # fracción CM del caso repulsivo (ver ec. umbral)
OUT = BASE / ("fig_paso_C_adhesion.pdf" if LANG == "es"
              else "fig_paso_C_adhesion_en.pdf")

aplicar_estilo()

# --- E_adh por N: |<Uc>| sobre las últimas 50 salidas, promediado en h0 --------
E_adh = {}
for N in NS:
    vals = []
    for d in sorted((BASE / f"adhesion_N{N}").iterdir()):
        dat = np.loadtxt(d / "bounce.dat", comments="#")
        arr = {c: dat[:, i] for i, c in enumerate(COLS[:dat.shape[1]])}
        vals.append(abs(arr["Uc"][-50:].mean()))
    E_adh[N] = (np.mean(vals), np.std(vals))

Nv = np.array(NS, float)
Ev = np.array([E_adh[N][0] for N in NS])
Es = np.array([E_adh[N][1] for N in NS])

# Ajustes: E_adh = C N^p  →  h0* = E_adh/(0.35 N g)
p, logC = np.polyfit(np.log(Nv), np.log(Ev), 1)
C = np.exp(logC)
h0s = Ev / (A_CM * Nv * G)

# --- Figura ---------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(ANCHO, 2.55))
fig.subplots_adjust(wspace=0.32)

# (a) trayectorias N=79
runs = sorted((BASE / "adhesion_N79").iterdir(),
              key=lambda x: float(x.name.split("_")[1]))
cmap = plt.cm.viridis(np.linspace(0.0, 0.85, len(runs)))
for d, col in zip(runs, cmap):
    dat = np.loadtxt(d / "bounce.dat", comments="#")
    h0 = d.name.split("_")[1]
    ax1.plot(dat[:, 0], dat[:, 1], color=col, lw=0.9,
             label=rf"$h_0={h0}\,\sigma$")
ax1.axhline(9.3, color="0.5", ls="--", lw=0.7)
ax1.text(118, 8.4, T["tope"], fontsize=7, color="0.35",
         ha="right", va="top")
ax1.set_xlim(0, 120)
ax1.set_ylim(5, 42)
ax1.set_xlabel(r"$t$ $(\tau)$")
ax1.set_ylabel(r"$Z_{\rm cm}$ $(\sigma)$")
ax1.legend(fontsize=6.5, loc="upper right", ncol=2,
           columnspacing=0.8, handletextpad=0.4, handlelength=1.2)
panel_label(ax1, "(a)")

# (b) E_adh(N) y h0*(N)
Ngrid = np.linspace(10, 100, 200)
ax2.errorbar(Nv, h0s, yerr=Es / (A_CM * Nv * G), fmt="o",
             color="#d62728", ms=5, capsize=2.5,
             label=T["umbral"])
ax2.plot(Ngrid, C * Ngrid**p / (A_CM * Ngrid * G), "k--", lw=1.1,
         label=rf"$E_{{\rm adh}} = {C:.1f}\,N^{{{p:.2f}}}$")
ax2.plot(Ngrid, h0s[-1] * (Ngrid / Nv[-1])**(-1 / 3), ":", color="#5B8DB8",
         lw=1.1, label=r"$\propto N^{-1/3}$")
ax2.set_xlabel(r"$N$")
ax2.set_ylabel(T["ylabel_b"])
ax2.set_xlim(0, 100)
ax2.legend(loc="upper right", fontsize=7)
panel_label(ax2, "(b)")

fig.savefig(OUT)
print(f"Figura guardada: {OUT}")
for N in NS:
    print(f"  N={N:3d}: E_adh = {E_adh[N][0]:6.1f} ± {E_adh[N][1]:.1f} eps"
          f"   E/N = {E_adh[N][0]/N:.2f}   h0* = {E_adh[N][0]/(A_CM*N*G):.0f}")
print(f"  Ajuste: E_adh = {C:.1f} N^{p:.3f}  →  h0* ∝ N^{p-1:.2f}")
