"""
fig_largo.py — Figura de la simulación larga (5000 τ) para el paper.

Genera data/outputs/terminado/fig_paso_B_convergencia.pdf:
  (a) Z_cm(t) − z0, primeros 1000 τ (inset: 5000 τ completos);
      línea roja discontinua = colapso a través de la red
  (b) alturas de apogeo h_n vs predicción h*_n = T_B/(Mg)
  (c) razón h_n / h*_n por rebote

Uso:  python3 scripts/fig_largo.py
"""

import numpy as np
from pathlib import Path
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, str(Path(__file__).parent))
from estilo_figuras import aplicar_estilo, panel_label, ANCHO, COLS, G

LANG = "en" if "--en" in sys.argv else "es"
T = {"es": {"colapso": "colapso", "medido": r"$h_n$ medido",
            "n_rebote": r"número de rebote $n$", "altura": r"altura $(\sigma)$",
            "primer": "primer rebote"},
     "en": {"colapso": "collapse", "medido": r"$h_n$ measured",
            "n_rebote": r"bounce number $n$", "altura": r"height $(\sigma)$",
            "primer": "first bounce"}}[LANG]

DATA = Path("data/outputs/terminado/cluster_N79_long/h0_28.15/bounce.dat")
OUT = Path("data/outputs/terminado/fig_paso_B_convergencia.pdf" if LANG == "es"
           else "data/outputs/terminado/fig_paso_B_convergencia_en.pdf")
M = 79.0

aplicar_estilo()

d = np.loadtxt(DATA, comments="#")
arr = {c: d[:, i] for i, c in enumerate(COLS[:d.shape[1]])}
t, Zcm, TB = arr["time"], arr["Zcm"], arr["TB"]

# nivel de contacto: media de los mínimos pre-colapso (misma convención
# que el análisis de la partición: z0 ≈ 11.5 σ)
mins, _ = find_peaks(-Zcm, prominence=0.5, distance=50)
mins_pre = mins[Zcm[mins] > Zcm[mins[0]] - 1.0]
z0 = float(Zcm[mins_pre].mean())
# colapso: primer instante con Zcm < z0 - 1 sigma
i_col = np.argmax(Zcm < z0 - 3.0)
t_col = t[i_col] if i_col > 0 else None

# apogeos pre-colapso
maxs, _ = find_peaks(Zcm, prominence=0.3, distance=50)
maxs = maxs[t[maxs] < (t_col if t_col else t[-1])]
h_n = Zcm[maxs] - z0
hstar_n = TB[maxs] / (M * G)
n = np.arange(1, len(maxs) + 1)

# --- Figura ---------------------------------------------------------------------
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(ANCHO, 5.6))
fig.subplots_adjust(hspace=0.42)

# (a) trayectoria
ax1.plot(t, Zcm - z0, color="#1f77b4", lw=0.7)
ax1.axhline(0, color="0.6", lw=0.6, ls="--")
if t_col:
    ax1.axvline(t_col, color="#d62728", ls="--", lw=1.0)
    ax1.text(t_col - 15, -7.5,
             rf"{T['colapso']} $t\approx{t_col:.0f}\,\tau$",
             fontsize=7.5, color="#d62728", ha="right")
ax1.set_xlim(0, 1000)
ax1.set_ylim(-12, 28)
ax1.set_xlabel(r"$t$ $(\tau)$")
ax1.set_ylabel(r"$Z_{\rm cm} - z_0$ $(\sigma)$")
panel_label(ax1, "(a)")
# inset: los 5000 tau completos
axi = ax1.inset_axes([0.55, 0.55, 0.42, 0.40])
axi.plot(t, Zcm - z0, color="#1f77b4", lw=0.4)
axi.axhline(0, color="0.6", lw=0.4, ls="--")
axi.set_xlim(0, 5000)
axi.tick_params(labelsize=6.5, pad=1.5)
axi.set_xticks([0, 2500, 5000])
axi.text(0.95, 0.85, r"$5000\,\tau$", transform=axi.transAxes,
         fontsize=6.5, ha="right")

# (b) alturas de apogeo
ax2.plot(n, h_n, "o-", color="#1f77b4", ms=4, lw=1.0,
         label=T["medido"])
ax2.plot(n, hstar_n, "s--", color="#d62728", ms=4, lw=1.0, mfc="none",
         label=r"$h_n^* = T_B^{(n)}/(Mg)$")
ax2.set_xlabel(T["n_rebote"])
ax2.set_ylabel(T["altura"])
ax2.set_xlim(0.5, len(n) + 0.5)
ax2.legend(loc="upper right")
panel_label(ax2, "(b)")

# (c) razón
ax3.plot(n, h_n / hstar_n, "D-", color="#7b2d8b", ms=3.5, lw=1.0)
ax3.axhline(1.0, color="0.5", ls="--", lw=0.8)
ax3.annotate(rf"{T['primer']}: {h_n[0]/hstar_n[0]:.2f}",
             (1, h_n[0] / hstar_n[0]), textcoords="offset points",
             xytext=(14, -9), fontsize=7.5)
ax3.set_xlabel(T["n_rebote"])
ax3.set_ylabel(r"$h_n / h_n^*$")
ax3.set_xlim(0.5, len(n) + 0.5)
ax3.set_ylim(0, 1.15)
panel_label(ax3, "(c)")

fig.savefig(OUT)
print(f"Figura guardada: {OUT}")
print(f"  z0={z0:.2f}, t_col={t_col}, rebotes={len(n)}")
print(f"  h1={h_n[0]:.2f}, h*1={hstar_n[0]:.2f}, ratio={h_n[0]/hstar_n[0]:.3f}")
print(f"  ratio final (n={len(n)}): {h_n[-1]/hstar_n[-1]:.2f}")
