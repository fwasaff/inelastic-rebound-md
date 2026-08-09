"""
fig_sweep_N.py — Figura del barrido en N para el paper (estilo revista).

Genera data/outputs/terminado/fig_alpha_TB_full_sweep.pdf:
  (a) α_TB vs v_imp para N ∈ {1,2,7,13,55,79,87}
  (b) ⟨α_TB⟩_plateau vs N, comparado con la fórmula de equipartición
      para Nc=290 (azul sólida) y Nc(N) = A·N^β ajustado (negra discontinua);
      inset: Nc_eff(N) y su ajuste de ley de potencias.

Uso:  python3 scripts/fig_sweep_N.py
"""

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, str(Path(__file__).parent))
from estilo_figuras import (aplicar_estilo, panel_label, cargar_run,
                            ANCHO, COLOR_N, MARKER_N)

LANG = "en" if "--en" in sys.argv else "es"
T = {"es": {"xlabel_b": r"$N$ (átomos en la bola)", "sin_plateau": "sin plateau", "nc_ef": r"$N_c^{\rm ef}$"},
     "en": {"xlabel_b": r"$N$ (atoms in the ball)", "sin_plateau": "no plateau", "nc_ef": r"$N_c^{\rm eff}$"}}[LANG]

BASE = Path("data/outputs")
NS = [1, 2, 7, 13, 55, 79, 87]
V_PLATEAU = 1.1
OUT = BASE / ("fig_alpha_TB_full_sweep.pdf" if LANG == "es"
              else "fig_alpha_TB_full_sweep_en.pdf")


def alpha_eq(N, Nc):
    return (3 * N - 3) / (3 * N - 3 + 3 * Nc)


aplicar_estilo()

# --- Cargar datos --------------------------------------------------------------
datos = {}
for N in NS:
    runs = []
    for d in sorted((BASE / f"cluster_N{N}").iterdir(),
                    key=lambda p: float(p.name.split("_")[1])):
        r = cargar_run(d / "bounce.dat", M=N)
        if r is not None:
            runs.append(r)
    datos[N] = runs

# --- Figura ---------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(ANCHO, 2.75))
fig.subplots_adjust(wspace=0.28, top=0.88)

# (a) α_TB vs v_imp
for N in NS:
    v = [r["v_imp"] for r in datos[N]]
    a = [r["a_TB"] for r in datos[N]]
    ax1.plot(v, a, marker=MARKER_N[N], color=COLOR_N[N],
             label=f"$N={N}$", lw=1.0, ms=3.8,
             mfc="none" if N == 13 else None)
ax1.axvline(V_PLATEAU, color="0.5", ls=":", lw=0.8)
ax1.text(V_PLATEAU + 0.03, 0.005, r"$v^*$", fontsize=8, color="0.35")
ax1.set_xlabel(r"$v_{\rm imp}$ $(\sigma/\tau)$")
ax1.set_ylabel(r"$\alpha_{T_B} = T_B^{(1)}/E_0$")
ax1.set_xlim(0.6, 2.1)
ax1.set_ylim(-0.01, 0.35)
# leyenda única para toda la figura, en una fila superior
handles, labels = ax1.get_legend_handles_labels()
fig.legend(handles, labels, ncol=7, loc="upper center",
           bbox_to_anchor=(0.5, 1.00), columnspacing=0.9,
           handletextpad=0.35, fontsize=7.5, frameon=False)
panel_label(ax1, "(a)")

# (b) plateau vs N + fórmula
Nn, mean_a, std_a = [], [], []
for N in NS:
    plateau = [r["a_TB"] for r in datos[N] if r["v_imp"] > V_PLATEAU]
    if plateau:
        Nn.append(N)
        mean_a.append(np.mean(plateau))
        std_a.append(np.std(plateau))
Nn, mean_a, std_a = map(np.array, (Nn, mean_a, std_a))

# Nc_eff por inversión de la fórmula (solo N con modos internos y plateau)
N_fit = np.array([13, 55, 79, 87])
a_fit = np.array([mean_a[list(Nn).index(n)] for n in N_fit])
Nc_eff = (3 * N_fit - 3) * (1 - a_fit) / (3 * a_fit)
# Ajuste de ley de potencias Nc = A N^beta
beta, logA = np.polyfit(np.log(N_fit), np.log(Nc_eff), 1)
A = np.exp(logA)

Ngrid = np.linspace(2, 95, 300)
ax2.plot(Ngrid, alpha_eq(Ngrid, 290.0), "-", color="#5B8DB8", lw=1.3,
         label=r"$N_c = 290$")
ax2.plot(Ngrid, alpha_eq(Ngrid, A * Ngrid**beta), "k--", lw=1.3,
         label=rf"$N_c = {A:.0f}\,N^{{{beta:.2f}}}$")
for N, m, s in zip(Nn, mean_a, std_a):
    ax2.errorbar(N, m, yerr=s, marker=MARKER_N[N], color=COLOR_N[N],
                 ms=5, capsize=2, lw=1.0,
                 mfc="none" if N == 13 else None)
# marcar N=13 como fuera de plateau
i13 = list(Nn).index(13)
ax2.annotate(T["sin_plateau"], (13, mean_a[i13] + std_a[i13]),
             textcoords="offset points", xytext=(2, 6), fontsize=7,
             color=COLOR_N[13], ha="left")
ax2.set_xlabel(T["xlabel_b"])
ax2.set_ylabel(r"$\langle\alpha_{T_B}\rangle_{\rm plateau}$")
ax2.set_xlim(0, 95)
ax2.set_ylim(0, 0.32)
ax2.legend(loc="upper left", bbox_to_anchor=(0.02, 0.99))
panel_label(ax2, "(b)", x=0.86)

# inset: Nc_eff(N)
axi = ax2.inset_axes([0.52, 0.12, 0.44, 0.38])
axi.plot(N_fit, Nc_eff, "ko", ms=3)
axi.plot(Ngrid, A * Ngrid**beta, "k--", lw=0.8)
axi.set_xlabel(r"$N$", fontsize=7, labelpad=1)
axi.set_ylabel(T["nc_ef"], fontsize=7, labelpad=1)
axi.tick_params(labelsize=6.5, pad=1.5)
axi.set_xlim(0, 95)

fig.savefig(OUT)
print(f"Figura guardada: {OUT}")
print(f"  Nc_eff = {dict(zip(N_fit, Nc_eff.round(0)))}")
print(f"  Ajuste: Nc = {A:.1f} N^{beta:.3f}")
for N, m, s in zip(Nn, mean_a, std_a):
    print(f"  N={N:3d}: <a_TB> = {m:.3f} ± {s:.3f}")
