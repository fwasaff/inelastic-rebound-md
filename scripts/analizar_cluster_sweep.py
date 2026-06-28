"""
Análisis del sweep N ∈ {1, 2, 7, 13, 55} + N=83 (reproducción).

Genera dos figuras:
  Fig B1 – α_TB(N) medido vs predicción de Gutiérrez: (3N-3)/(3N-3+3Nc)
  Fig B2 – e₁(N) y α_surf(N) vs número de átomos de la bola

Partición energética directa (3 canales, suman 1):
  α_CM   = M·g·h₁ / E₀
  α_surf = ΔE_surf / E₀  (medido de KS + PE_surf)
  α_ball = 1 − α_CM − α_surf  (total interno bola = α_TB + α_VB)
  α_TB   = TB₁ / E₀  (canal cinético — comparar con fórmula Gutiérrez)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from pathlib import Path

G      = 0.026
NS     = [1, 2, 7, 13, 55, 83]
H0S    = [11.03, 18.64, 28.15, 41.47, 58.6, 79.53]
SWEEP  = Path("data/outputs/cluster_sweep")
REPRO  = Path("data/outputs/reprod_2012")
FIGDIR = SWEEP

COLORES = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b"]

plt.rcParams.update({
    "font.size": 11, "axes.labelsize": 12, "legend.fontsize": 9,
    "figure.dpi": 150, "axes.grid": True, "grid.alpha": 0.3,
})


# ---------------------------------------------------------------------------
# Extracción
# ---------------------------------------------------------------------------

def get_particion(path, N):
    """Devuelve dict con α_CM, α_TB, α_ball, α_surf, v_imp, e1."""
    d = np.loadtxt(path, comments="#")
    Zcm=d[:,1]; TB=d[:,4]; KS=d[:,9]; PE_surf=d[:,10]
    M = float(N)
    mins,_=find_peaks(-Zcm, prominence=0.5, distance=50)
    maxs,_=find_peaks( Zcm, prominence=0.5, distance=50)
    if not len(mins) or not len(maxs): return None
    z0=Zcm[mins].mean(); h0_eff=float(Zcm[0]-z0); E0=M*G*h0_eff
    if E0 < 0.05: return None
    h1  = float(Zcm[maxs[0]]-z0)
    TB1 = float(TB[maxs[0]])
    dEs = float((KS[maxs[0]]+PE_surf[maxs[0]])-(KS[0]+PE_surf[0]))
    a_CM   = M*G*h1/E0
    a_surf = dEs/E0
    a_ball = 1.0 - a_CM - a_surf
    a_TB   = TB1/E0
    return {
        "v":      float(np.sqrt(2*G*h0_eff)),
        "e1":     float(np.sqrt(max(0, h1/h0_eff))),
        "a_CM":   a_CM,
        "a_surf": a_surf,
        "a_ball": a_ball,
        "a_TB":   a_TB,
        "a_VB":   a_ball - a_TB,
    }


def cargar_N(N):
    """Devuelve lista de dicts (uno por h0) para un N dado."""
    resultado = []
    for h0 in H0S:
        p = REPRO/f"h0_{h0}/bounce.dat" if N==83 else SWEEP/f"N{N}/h0_{h0}/bounce.dat"
        if not p.exists(): continue
        r = get_particion(p, N)
        if r: resultado.append(r)
    return resultado


# ---------------------------------------------------------------------------
# Fig B1 – α_TB(N) vs predicción de Gutiérrez
# ---------------------------------------------------------------------------

def fig_aTB_vs_N(datos_N: dict, outdir: Path):
    """
    datos_N: {N: [lista de dicts por h0]}
    Promedia solo sobre v > 0.9 σ/τ (régimen saturado).
    """
    # Ajustar Nc desde N=83
    hi83 = [r for r in datos_N[83] if r["v"] > 0.9]
    aTB_83 = np.mean([r["a_TB"] for r in hi83])
    Nc = (3*83-3) * (1/aTB_83 - 1) / 3

    N_arr=[]; aTB_m=[]; aTB_s=[]; aTB_p=[]; e1_m=[]
    for N, rows in sorted(datos_N.items()):
        hi = [r for r in rows if r["v"] > 0.9]
        if not hi: continue
        aTB = np.mean([r["a_TB"] for r in hi])
        N_arr.append(N); aTB_m.append(aTB)
        aTB_s.append(np.std([r["a_TB"] for r in hi]))
        aTB_p.append((3*N-3)/(3*N-3+3*Nc) if N>1 else 0.0)
        e1_m.append(np.mean([r["e1"] for r in hi]))

    fig, ax = plt.subplots(figsize=(7, 5))

    # Curva predicción
    N_cont = np.linspace(1, 90, 300)
    pred = (3*N_cont-3)/(3*N_cont-3+3*Nc)
    ax.plot(N_cont, pred, "k--", lw=1.5,
            label=f"$\\alpha_{{T_B}}(N) = \\frac{{3N-3}}{{3N-3+3N_c}}$,  $N_c = {Nc:.0f}$")

    # Puntos medidos con barras de error
    ax.errorbar(N_arr, aTB_m, yerr=aTB_s, fmt="none", ecolor="gray", capsize=3, lw=1.0)
    sc = ax.scatter(N_arr, aTB_m, s=90, c=COLORES[:len(N_arr)],
                    edgecolors="k", lw=0.7, zorder=5)
    for n, m, p in zip(N_arr, aTB_m, aTB_p):
        ax.annotate(f"$N={n}$", (n, m), xytext=(6, 4), textcoords="offset points",
                    fontsize=8.5)

    ax.set_xlabel("Número de átomos de la bola $N$")
    ax.set_ylabel("$\\alpha_{T_B} = T_{B1}/E_0$ (canal cinético interno)")
    ax.set_title("Dependencia con $N$ del canal cinético interno\n"
                 "promedio sobre $v_0 > 0.9\\,\\sigma/\\tau$ (régimen saturado)")
    ax.set_xlim(0, 92); ax.set_ylim(-0.01, 0.32)
    ax.legend(loc="upper left")

    out = outdir / "fig_B1_aTB_vs_N.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"  Guardado: {out}")
    return Nc, N_arr, aTB_m, aTB_p


# ---------------------------------------------------------------------------
# Fig B2 – e₁(N) y α_surf(N) — dos paneles
# ---------------------------------------------------------------------------

def fig_e1_asurf_vs_N(datos_N: dict, outdir: Path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    N_arr=[]; e1_m=[]; e1_s=[]; as_m=[]; as_s=[]
    for N, rows in sorted(datos_N.items()):
        hi = [r for r in rows if r["v"] > 0.9]
        if not hi: continue
        N_arr.append(N)
        e1_m.append(np.mean([r["e1"] for r in hi]))
        e1_s.append(np.std([r["e1"] for r in hi]))
        as_m.append(np.mean([r["a_surf"] for r in hi]))
        as_s.append(np.std([r["a_surf"] for r in hi]))

    for ax, ys, yer, ylabel, title in [
        (ax1, e1_m, e1_s,
         "$e_1 = \\sqrt{h_1/h_0}$ (COR primer rebote)",
         "Coeficiente de restitución vs $N$"),
        (ax2, as_m, as_s,
         "$\\alpha_\\mathrm{surf}$ (fonones superficie / $E_0$)",
         "Canal superficial vs $N$"),
    ]:
        ax.errorbar(N_arr, ys, yerr=yer, fmt="none", ecolor="gray", capsize=3, lw=1.0)
        ax.scatter(N_arr, ys, s=90, c=COLORES[:len(N_arr)],
                   edgecolors="k", lw=0.7, zorder=5)
        for n, y in zip(N_arr, ys):
            ax.annotate(f"$N={n}$", (n, y), xytext=(4, 5),
                        textcoords="offset points", fontsize=8.5)
        ax.set_xlabel("Número de átomos $N$")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xlim(0, 92)

    ax1.set_ylim(0.3, 1.05)
    ax1.axhline(1.0, color="k", ls=":", lw=0.8, alpha=0.5)
    ax2.set_ylim(0, 0.45)

    fig.tight_layout()
    out = outdir / "fig_B2_e1_asurf_vs_N.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"  Guardado: {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Cargando datos...")
    datos_N = {N: cargar_N(N) for N in NS}

    # Tabla resumen
    print("\n=== TABLA: α_TB(N) medido vs predicción Gutiérrez ===")
    hi83 = [r for r in datos_N[83] if r["v"] > 0.9]
    aTB_83 = np.mean([r["a_TB"] for r in hi83])
    Nc = (3*83-3) * (1/aTB_83 - 1) / 3
    print(f"  Nc = {Nc:.0f}  (ajustado desde N=83, α_TB={aTB_83:.3f})\n")
    print(f"{'N':>5}  {'α_TB_pred':>10}  {'α_TB_meas':>10}  {'α_surf':>8}  {'e1':>6}")
    print("-"*50)
    for N in NS:
        hi = [r for r in datos_N[N] if r["v"] > 0.9]
        if not hi: continue
        pred = (3*N-3)/(3*N-3+3*Nc) if N>1 else 0.0
        print(f"{N:>5}  {pred:>10.3f}  "
              f"{np.mean([r['a_TB'] for r in hi]):>10.3f}  "
              f"{np.mean([r['a_surf'] for r in hi]):>8.3f}  "
              f"{np.mean([r['e1'] for r in hi]):>6.3f}")

    print("\nGenerando figuras...")
    fig_aTB_vs_N(datos_N, FIGDIR)
    fig_e1_asurf_vs_N(datos_N, FIGDIR)
    print(f"\nFiguras guardadas en: {FIGDIR}")


if __name__ == "__main__":
    main()
