"""
fig_Nc_escalas.py — Las tres escalas de N_c medidas directamente (M1+M2).

Panel (a): N_c^geom (huella de contacto) y N_c^part (nº de participación de la
           KE de superficie) medidos en el dump vs N, más N_c^eff (capacidad
           neta de la superficie, de la partición) — muestra la separación
           N_c^eff << N_c^geom << N_c^part.
Panel (b): f_ball(N) con barras de error y el ajuste de equipartición sobre la
           energía disipada con N_c^eff constante (reemplaza el límite 1/2).
"""
import numpy as np
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from estilo_figuras import aplicar_estilo, panel_label, COLOR_N, MARKER_N, ANCHO
from analizar_orient import main as run_orient
import matplotlib.pyplot as plt

# --- N_c medido directamente (medir_Nc.py --sweep, v_imp≈1.47) ----------------
NC = {  # N: (Nc_geom, Nc_part)
    7:  (47, 284.4),
    13: (43, 317.6),
    55: (67, 340.9),
    79: (65, 413.7),
    87: (93, 397.4),
}


def main():
    aplicar_estilo()
    plateau = run_orient()          # reusa el análisis con barras de error

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(ANCHO, 2.8))
    fig.subplots_adjust(wspace=0.32)

    # (a) escalas de N_c
    Ns = np.array(sorted(NC))
    geom = np.array([NC[n][0] for n in Ns])
    part = np.array([NC[n][1] for n in Ns])
    ax1.plot(Ns, part, "o-", color="#C62828", ms=4, lw=1.0,
             label=r"$N_c^{\rm part}$ (participación KE)")
    ax1.plot(Ns, geom, "s-", color="#1565C0", ms=4, lw=1.0,
             label=r"$N_c^{\rm geom}$ (huella $<r_c$)")
    # N_c^eff de la partición (constante ajustado)
    ax1.axhspan(17, 33, color="0.85", zorder=0)
    ax1.axhline(25, color="k", ls="--", lw=1.0,
                label=r"$N_c^{\rm ef}\approx25$ (partición)")
    ax1.plot(Ns, Ns, ":", color="0.5", lw=0.9, label=r"$N_c=N$ (hipótesis previa)")
    ax1.set_xlabel(r"$N$ (átomos en la bola)")
    ax1.set_ylabel(r"$N_c$ (átomos de superficie)")
    ax1.set_yscale("log")
    ax1.set_xlim(0, 95)
    ax1.legend(loc="center right", fontsize=6.3)
    panel_label(ax1, "(a)")

    # (b) f_ball con ajuste
    from scipy.optimize import curve_fit
    NS = [n for n in [2, 7, 13, 55, 79, 87] if n in plateau]
    fb = np.array([plateau[n]["fball"][0] for n in NS])
    fbs = np.array([plateau[n]["fball"][1] for n in NS])
    NSa = np.array(NS)
    model = lambda N, Nc: (3*N-3)/(3*N-3+3*Nc)
    m = NSa >= 7
    popt, _ = curve_fit(model, NSa[m], fb[m], p0=[20],
                        sigma=np.maximum(fbs[m], 1e-3), absolute_sigma=True)
    Nc_fit = popt[0]
    Ng = np.logspace(np.log10(1.6), np.log10(400), 300)
    ax2.plot(Ng, model(Ng, Nc_fit), "k--", lw=1.2,
             label=rf"$\frac{{3N-3}}{{3N-3+3N_c^{{\rm ef}}}}$, $N_c^{{\rm ef}}={Nc_fit:.0f}$")
    ax2.plot(Ng, 0.5*np.ones_like(Ng), ":", color="0.5", lw=0.9,
             label=r"límite previo $1/2$")
    for n in NS:
        ax2.errorbar(n, plateau[n]["fball"][0], yerr=plateau[n]["fball"][1],
                     marker=MARKER_N[n], color=COLOR_N[n], ms=5, capsize=2, lw=1.0,
                     mfc="none" if n == 13 else None)
    ax2.set_xscale("log")
    ax2.set_xlabel(r"$N$ (átomos en la bola)")
    ax2.set_ylabel(r"$f_{\rm ball}=\dfrac{\alpha_{\rm ball}}{\alpha_{\rm ball}+\alpha_{\rm surf}}$")
    ax2.set_xlim(1.6, 400); ax2.set_ylim(0, 1.05)
    ax2.set_xticks([2, 7, 13, 55, 87]); ax2.set_xticklabels(["2","7","13","55","87"])
    ax2.axhline(1.0, color="0.6", ls=":", lw=0.7)
    ax2.legend(loc="lower right", fontsize=6.8)
    panel_label(ax2, "(b)")

    out = Path("data/outputs/fig_Nc_escalas.pdf")
    fig.savefig(out); print(f"Figura guardada: {out}")
    print(f"  N_c^eff (ajuste f_ball, N>=7) = {Nc_fit:.1f}")


if __name__ == "__main__":
    main()
