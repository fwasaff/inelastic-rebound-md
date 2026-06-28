"""
Verificación de la anharmonicidad LJ: ΔV_B vs T_{B1}

Para cada N ∈ {2,7,13,55,83} y las 6 alturas:
  · TB1  = KE interna del CM de la bola en el primer apogeo
  · ΔVB  = VB(apogeo) − VB(t=0)  (cambio en PE intra-bola)

Ajuste global: ΔVB = a·TB + b·TB²  (forzado a pasar por el origen)

Predicción armónica:          a = 1, b = 0  (equipartición)
Predicción LJ anarmónica:     a ≈ 1, b > 0  (cuártico domina al cúbico)

Genera fig_C1_dVB_vs_TB.pdf en data/outputs/cluster_sweep/
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from pathlib import Path

G    = 0.026
H0S  = [11.03, 18.64, 28.15, 41.47, 58.6, 79.53]
PROJ = Path(__file__).resolve().parent.parent
SWEEP = PROJ / "data/outputs/cluster_sweep"
REPRO = PROJ / "data/outputs/reprod_2012"

COLORES = {2:"#ff7f0e", 7:"#2ca02c", 13:"#d62728", 55:"#9467bd", 83:"#1f77b4"}
MARKERS = {2:"s", 7:"^", 13:"D", 55:"v", 83:"o"}

plt.rcParams.update({
    "font.size": 11, "axes.labelsize": 12, "legend.fontsize": 9,
    "figure.dpi": 150, "axes.grid": True, "grid.alpha": 0.3,
})


def cargar_TB_dVB(N):
    """Devuelve arrays (TB1, ΔVB) para un N dado, todas las alturas."""
    TBs, dVBs = [], []
    for h0 in H0S:
        path = REPRO/f"h0_{h0}/bounce.dat" if N == 83 \
               else SWEEP/f"N{N}/h0_{h0}/bounce.dat"
        if not path.exists():
            continue
        d = np.loadtxt(path, comments="#")
        Zcm = d[:, 1]; TB = d[:, 4]; VB = d[:, 5]
        mins, _ = find_peaks(-Zcm, prominence=0.5, distance=50)
        maxs, _ = find_peaks( Zcm, prominence=0.5, distance=50)
        if not len(maxs):
            continue
        TB1 = float(TB[maxs[0]])
        dVB = float(VB[maxs[0]]) - float(VB[0])
        if TB1 > 0.1:
            TBs.append(TB1)
            dVBs.append(dVB)
    return np.array(TBs), np.array(dVBs)


def teorica_LJ():
    """
    Estima los coeficientes de la corrección anarmónica de 1er orden para LJ.

    V(r) = V(r₀) + ½k·x² + (γ/6)·x³ + (λ/24)·x⁴ + ...
    con x = r − r₀, r₀ = 2^(1/6)σ

    A orden (kBT)², la corrección a ΔVB−TB por par es:
       c_pair = λ/(8k²) − 5γ²/(12k³)    [en unidades de 1/ε]

    El coeficiente cuártico domina al cúbico por ~1%.
    """
    r0    = 2 ** (1 / 6)                  # mínimo LJ en unidades de σ
    k     = 72 / r0 ** 2                  # V''(r₀) [ε/σ²]
    gamma = -756 / r0 ** 3                # V'''(r₀) [ε/σ³]  < 0
    lam   = 26712 / r0 ** 4              # V''''(r₀) [ε/σ⁴] > 0
    c_cubic   = -5 * gamma**2 / (12 * k**3)   # [1/ε]  < 0
    c_quartic =  lam / (8 * k**2)              # [1/ε]  > 0
    return c_cubic, c_quartic, c_cubic + c_quartic


def main():
    # --- Cargar datos -----------------------------------------------------------
    data = {N: cargar_TB_dVB(N) for N in [2, 7, 13, 55, 83]}

    def quad(x, a, b):
        return a * x + b * x**2

    # --- Ajuste por N -----------------------------------------------------------
    print("=== ΔVB = a·TB + b·TB²  (ajuste por N) ===")
    print(f"{'N':>4}  {'a':>7} {'b':>9}  {'R²':>7}")
    print("-" * 40)
    for N in [13, 55, 83]:
        TB, dVB = data[N]
        if len(TB) < 3:
            continue
        popt, _ = curve_fit(quad, TB, dVB, p0=[1.0, 0.01])
        resid = dVB - quad(TB, *popt)
        SS_tot = np.sum((dVB - dVB.mean()) ** 2)
        R2 = 1 - np.sum(resid**2) / SS_tot
        print(f"{N:>4}  {popt[0]:>7.4f} {popt[1]:>9.5f}  {R2:>7.5f}")

    # --- Ajuste global (N=55 y 83, v > 0.9 σ/τ) --------------------------------
    TB_all = np.concatenate([data[N][0] for N in [55, 83]])
    dVB_all = np.concatenate([data[N][1] for N in [55, 83]])
    popt_g, pcov_g = curve_fit(quad, TB_all, dVB_all, p0=[1.0, 0.01])
    perr_g = np.sqrt(np.diag(pcov_g))
    a_fit, b_fit = popt_g
    SS_tot = np.sum((dVB_all - dVB_all.mean()) ** 2)
    R2_q = 1 - np.sum((dVB_all - quad(TB_all, *popt_g))**2) / SS_tot
    a_lin = np.dot(dVB_all, TB_all) / np.dot(TB_all, TB_all)
    R2_l = 1 - np.sum((dVB_all - a_lin * TB_all)**2) / SS_tot

    print("\n=== Ajuste global (N=55 + N=83) ===")
    print(f"  a = {a_fit:.4f} ± {perr_g[0]:.4f}  (armónico: a=1)")
    print(f"  b = {b_fit:.5f} ± {perr_g[1]:.5f}  (anarmónico: b>0)")
    print(f"  R² cuadrático = {R2_q:.5f}")
    print(f"  R² lineal     = {R2_l:.5f}")
    print(f"  Significancia de b: {abs(b_fit)/perr_g[1]:.1f}σ")

    # --- Predicción teórica (perturbación 1er orden) ----------------------------
    c_cub, c_qua, c_net = teorica_LJ()
    print("\n=== Predicción teórica (perturbación 1er orden) ===")
    print(f"  c_cubic   = {c_cub:+.4f}/ε   (V''' < 0 → reduce ΔVB)")
    print(f"  c_quartic = {c_qua:+.4f}/ε   (V'''' > 0 → aumenta ΔVB)")
    print(f"  c_neto    = {c_net:+.5f}/ε   (cuártico domina ~1%)")
    print(f"  → La perturbación de 1er orden es suficiente para el SIGNO")
    print(f"    pero no para el valor de b (régimen kBT/ε ~ 0.2, no << 1)")

    # --- Figura -----------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    TB_range = np.linspace(0, max(data[83][0]) * 1.05, 300)

    # Panel izq: ΔVB vs TB
    for N in [2, 7, 13, 55, 83]:
        TB, dVB = data[N]
        if not len(TB):
            continue
        ax1.scatter(TB, dVB, s=55, color=COLORES[N], marker=MARKERS[N],
                    edgecolors="k", lw=0.5, label=f"$N={N}$", zorder=4)

    ax1.plot(TB_range, quad(TB_range, *popt_g), "k-", lw=2.2, zorder=6,
             label=f"$a={a_fit:.3f},\\;b={b_fit:.4f}$")
    ax1.plot(TB_range, TB_range, "k:", lw=1.2, alpha=0.5,
             label="armónico puro  ($a=1, b=0$)")

    ax1.set_xlabel("$T_{B1}$  (KE interna en 1.º apogeo, $\\varepsilon$)")
    ax1.set_ylabel("$\\Delta V_B$  (cambio PE intra-bola, $\\varepsilon$)")
    ax1.set_title("$\\Delta V_B$ vs $T_{B1}$ — todos los $N$ y $h_0$")
    ax1.legend(fontsize=8.5, loc="upper left")
    ax1.set_xlim(-0.5, max(data[83][0]) * 1.05)

    # Panel der: ΔVB/TB vs TB — razón supraunitaria
    for N in [2, 7, 13, 55, 83]:
        TB, dVB = data[N]
        mask = TB > 0.5
        if not mask.any():
            continue
        ax2.scatter(TB[mask], dVB[mask] / TB[mask], s=55,
                    color=COLORES[N], marker=MARKERS[N],
                    edgecolors="k", lw=0.5, label=f"$N={N}$", zorder=4)

    ax2.plot(TB_range[TB_range > 0.1],
             quad(TB_range[TB_range > 0.1], *popt_g) / TB_range[TB_range > 0.1],
             "k-", lw=2.2, label=f"$\\Delta V_B/T_B = {a_fit:.3f} + {b_fit:.4f}\\,T_B$")
    ax2.axhline(1.0, color="k", ls=":", lw=1.2, alpha=0.5,
                label="armónico puro ($\\Delta V_B = T_B$)")

    ax2.set_xlabel("$T_{B1}$  ($\\varepsilon$)")
    ax2.set_ylabel("$\\Delta V_B / T_{B1}$")
    ax2.set_title("Razón $\\Delta V_B/T_B$ — supraunitaria  (LJ anarmónico)")
    ax2.legend(fontsize=8.5, loc="upper left")
    ax2.set_xlim(-0.5, max(data[83][0]) * 1.05)
    ax2.set_ylim(0.5, 2.2)

    fig.tight_layout()
    out = SWEEP / "fig_C1_dVB_vs_TB.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"\nFigura guardada: {out}")


if __name__ == "__main__":
    main()
