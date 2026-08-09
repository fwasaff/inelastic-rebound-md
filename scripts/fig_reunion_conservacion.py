"""
fig_reunion_conservacion.py — Figuras de control de calidad para la reunión:
  1. fig_conservacion.pdf : (a) drift de energía total conservada,
                            (b) canales de energía vs tiempo (primer rebote)
  2. fig_termalizacion.pdf: (a) Zcm(t) y T_B(t) — escalera de termalización,
                            (b) T_ef vs v_imp con banda de fusión

Datos: data/outputs/activo/viz_N79_h28 (N=79, h0=28.15σ, v_imp≈1.16) y
       data/outputs/terminado/cluster_orient/tabla_orient.txt (N=79, 5 orientaciones).

Uso:  python3 scripts/fig_reunion_conservacion.py
Salida: docs/figuras_reunion/
"""

import re
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).parent))
from estilo_figuras import aplicar_estilo, ANCHO, plt

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "data/outputs/activo/viz_N79_h28"
OUT = ROOT / "docs/figuras_reunion"
OUT.mkdir(exist_ok=True)

M, G, N = 79.0, 0.026, 79
Z0 = 11.7          # nivel de contacto efectivo [sigma]

aplicar_estilo()

# --- leer termo del log (Step Time Zcm Vcmz KB TB VB Uc UF TotEng) ----------
rows = []
inrun = False
for line in (RUN / "lammps.log").read_text().splitlines():
    if line.startswith("   Step          Time"):
        inrun = True
        continue
    if inrun:
        p = line.split()
        if len(p) == 10:
            try:
                rows.append([float(x) for x in p])
            except ValueError:
                inrun = False
        else:
            inrun = False
log = np.array(rows)
t_log = log[:, 1]
E_cons = log[:, 9] + log[:, 8]          # TotEng + U_grav
E0 = M * G * (log[0, 2] - Z0)           # energía de impacto

# --- leer bounce.dat --------------------------------------------------------
# time Zcm Vcmz KB TB VB Uc UF KB3d KS PE_surf
d = np.loadtxt(RUN / "bounce.dat")
t, Zcm, KB, TB, VB, Uc, KS, PEs = (d[:, 0], d[:, 1], d[:, 3], d[:, 4],
                                   d[:, 5], d[:, 6], d[:, 9], d[:, 10])
E_surf = KS + (PEs - PEs[0])
dVB = VB - VB[0]

i_imp = np.argmin(d[:, 2])              # velocidad más negativa ~ impacto
t_imp = t[i_imp]

# ============================================================================
# Figura 1: conservación
# ============================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(ANCHO, 2.5))

ax1.plot(t_log, (E_cons - E_cons[0]) / E0 * 1e4, color="#1f77b4")
ax1.axhline(0, color="0.6", lw=0.6, ls=":")
ax1.set_xlabel(r"$t\ [\tau]$")
ax1.set_ylabel(r"$(E(t)-E(0))/E_0\ \times 10^{4}$")
ax1.text(0.03, 0.93, "(a)", transform=ax1.transAxes, weight="bold")
ax1.text(0.97, 0.05,
         (r"$|\Delta E|_{\max}=%.3f\,\varepsilon=%.1f\times10^{-4}E_0$"
          % (np.abs(E_cons - E_cons[0]).max(),
             np.abs(E_cons - E_cons[0]).max() / E0 * 1e4)),
         transform=ax1.transAxes, ha="right", fontsize=7.5)

sel = (t > t_imp - 12) & (t < t_imp + 55)
ax2.plot(t[sel], KB[sel] / E0, label=r"$K_{\rm cm}$", color="#1f77b4")
ax2.plot(t[sel], TB[sel] / E0, label=r"$T_B$", color="#d62728")
ax2.plot(t[sel], dVB[sel] / E0, label=r"$\Delta V_B$", color="#ff7f0e")
ax2.plot(t[sel], E_surf[sel] / E0, label=r"$\Delta E_{\rm surf}$",
         color="#2ca02c")
ax2.plot(t[sel], Uc[sel] / E0, label=r"$U_c$", color="0.5", ls="--")
ax2.axvline(t_imp, color="0.8", lw=0.6)
ax2.set_xlabel(r"$t\ [\tau]$")
ax2.set_ylabel(r"$E/E_0$")
ax2.legend(ncol=2, loc="upper right")
ax2.text(0.03, 0.93, "(b)", transform=ax2.transAxes, weight="bold")

fig.tight_layout()
fig.savefig(OUT / "fig_conservacion.pdf")
plt.close(fig)
print("OK fig_conservacion.pdf  |dE|max = %.4f eps (%.1e E0)"
      % (np.abs(E_cons - E_cons[0]).max(),
         np.abs(E_cons - E_cons[0]).max() / E0))

# ============================================================================
# Figura 2: termalización
# ============================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(ANCHO, 2.5))

ax1.plot(t, Zcm - Z0, color="0.4", lw=0.8, label=r"$z_{\rm cm}-z_0$")
ax1.set_xlabel(r"$t\ [\tau]$")
ax1.set_ylabel(r"$z_{\rm cm}-z_0\ [\sigma]$")
ax1b = ax1.twinx()
Tef = 2 * TB / (3 * N - 3)
ax1b.plot(t, Tef, color="#d62728", lw=1.0, label=r"$T_{\rm ef}$")
ax1b.set_ylabel(r"$T_{\rm ef}=2T_B/(3N-3)\ [\varepsilon]$",
                color="#d62728")
ax1b.tick_params(axis="y", colors="#d62728")
ax1.text(0.03, 0.93, "(a)", transform=ax1.transAxes, weight="bold")

# --- T_ef vs v_imp desde tabla_orient (N=79) --------------------------------
vs, atb, satb = [], [], []
for line in (ROOT / "data/outputs/terminado/cluster_orient/tabla_orient.txt"
             ).read_text().splitlines():
    p = line.split()
    if len(p) == 10 and p[0] == "79":
        vs.append(float(p[3]))
        m = re.match(r"([-\d.]+)±([\d.]+)", p[6])
        atb.append(float(m.group(1)))
        satb.append(float(m.group(2)))
vs, atb, satb = map(np.array, (vs, atb, satb))
# T_B = a_TB * E0, E0 = M v^2/2  ->  T_ef = a_TB M v^2 / (3N-3)
Tef_v = atb * M * vs**2 / (3 * N - 3)
sTef_v = satb * M * vs**2 / (3 * N - 3)

ax2.errorbar(vs, Tef_v, yerr=sTef_v, fmt="o", color="#d62728",
             capsize=2, label=r"$N=79$ (5 orientaciones)")
ax2.axhspan(0.30, 0.35, color="#ff7f0e", alpha=0.25, lw=0,
            label=r"fusión cúmulo LJ $\sim$80 át.")
ax2.axvline(1.7, color="0.5", ls=":", lw=0.8)
ax2.text(1.68, 0.02, r"límite plateau $v\lesssim1.7$", rotation=90,
         fontsize=7, ha="right", va="bottom", color="0.35")
ax2.set_xlabel(r"$v_{\rm imp}\ [\sigma/\tau]$")
ax2.set_ylabel(r"$T_{\rm ef}$ en el apogeo $[\varepsilon]$")
ax2.set_ylim(0, 0.40)
ax2.legend(loc="upper left")
ax2.text(0.9, 0.93, "(b)", transform=ax2.transAxes, weight="bold")

fig.tight_layout()
fig.savefig(OUT / "fig_termalizacion.pdf")
plt.close(fig)
print("OK fig_termalizacion.pdf")
