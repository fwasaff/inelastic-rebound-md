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

    # Curva predicción en escala log
    N_cont = np.logspace(np.log10(1.2), np.log10(90), 300)
    pred = (3*N_cont-3)/(3*N_cont-3+3*Nc)
    ax.plot(N_cont, pred, "k--", lw=1.5,
            label=f"$\\alpha_{{T_B}}(N) = \\frac{{3N-3}}{{3N-3+3N_c}}$,  $N_c = {Nc:.0f}$")

    # Puntos medidos con barras de error
    ax.errorbar(N_arr, aTB_m, yerr=aTB_s, fmt="none", ecolor="gray", capsize=3, lw=1.0)
    ax.scatter(N_arr, aTB_m, s=90, c=COLORES[:len(N_arr)],
               edgecolors="k", lw=0.7, zorder=5)

    # Etiquetas: desplazamiento diferenciado para N=1 y N=2
    offsets = {1: (-18, 6), 2: (6, -14), 7: (6, 4), 13: (6, 4), 55: (6, 4), 83: (6, 4)}
    for n, m in zip(N_arr, aTB_m):
        dx, dy = offsets.get(n, (6, 4))
        ax.annotate(f"$N={n}$", (n, m), xytext=(dx, dy),
                    textcoords="offset points", fontsize=8.5)

    ax.set_xscale("log")
    ax.set_xlabel("Número de átomos de la bola $N$")
    ax.set_ylabel("$\\alpha_{T_B} = T_{B1}/E_0$ (canal cinético interno)")
    ax.set_xlim(0.8, 110); ax.set_ylim(-0.01, 0.32)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}" if x >= 1 else ""))
    ax.set_xticks([1, 2, 7, 13, 55, 83])
    ax.legend(loc="upper left")
    # sin título: va en caption

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

    offsets_B2 = {1: (-20, 5), 2: (5, -14), 7: (5, 5), 13: (5, 5), 55: (5, 5), 83: (5, 5)}
    for ax, ys, yer, ylabel in [
        (ax1, e1_m, e1_s,  "$e_1 = \\sqrt{h_1/h_0}$ (COR primer rebote)"),
        (ax2, as_m, as_s,  "$\\alpha_\\mathrm{surf}$ (fonones superficie / $E_0$)"),
    ]:
        ax.errorbar(N_arr, ys, yerr=yer, fmt="none", ecolor="gray", capsize=3, lw=1.0)
        ax.scatter(N_arr, ys, s=90, c=COLORES[:len(N_arr)],
                   edgecolors="k", lw=0.7, zorder=5)
        for n, y in zip(N_arr, ys):
            dx, dy = offsets_B2.get(n, (5, 5))
            ax.annotate(f"$N={n}$", (n, y), xytext=(dx, dy),
                        textcoords="offset points", fontsize=8.5)
        ax.set_xscale("log")
        ax.set_xlabel("Número de átomos $N$")
        ax.set_ylabel(ylabel)
        ax.set_xlim(0.8, 110)
        ax.set_xticks([1, 2, 7, 13, 55, 83])
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}" if x >= 1 else ""))

    ax1.set_ylim(0.3, 1.05)
    ax1.axhline(1.0, color="k", ls=":", lw=0.8, alpha=0.5)
    ax2.set_ylim(0, 0.45)
    # Anotar máximo de α_surf en N=13
    idx13 = N_arr.index(13)
    ax2.annotate("máximo", (13, as_m[idx13]),
                 xytext=(20, -18), textcoords="offset points",
                 fontsize=8, color="#d62728",
                 arrowprops=dict(arrowstyle="->", color="#d62728", lw=0.8))

    fig.tight_layout()
    out = outdir / "fig_B2_e1_asurf_vs_N.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"  Guardado: {out}")


# ---------------------------------------------------------------------------
# Fig B3 – comparación α_TB vs α_ball con dos Nc: opción 2
# ---------------------------------------------------------------------------

def fig_aTB_vs_aball(datos_N: dict, outdir: Path):
    """
    Dos paneles lado a lado:
      Izq: α_TB(N) con Nc ajustado desde TB (Nc≈304)
      Der: α_ball(N) = α_TB+α_VB con Nc ajustado desde α_ball (Nc≈82≈N_ball)

    N=7 es problemático (E0 muy pequeño, relación señal/ruido baja);
    se marca con marcador abierto y nota.
    """
    V_CUT = 0.9

    # Nc desde α_TB(83)
    hi83 = [r for r in datos_N[83] if r["v"] > V_CUT]
    aTB_83   = np.mean([r["a_TB"]   for r in hi83])
    aball_83 = np.mean([r["a_ball"] for r in hi83])
    Nc_TB   = (3*83-3) * (1/aTB_83   - 1) / 3
    Nc_ball = (3*83-3) * (1/aball_83 - 1) / 3

    print(f"\n  Nc(α_TB)   = {Nc_TB:.0f}  (ajustado desde α_TB={aTB_83:.3f})")
    print(f"  Nc(α_ball) = {Nc_ball:.0f}  (ajustado desde α_ball={aball_83:.3f})")
    print(f"\n{'N':>5}  {'α_TB_pred':>10}  {'α_TB':>8}  {'α_ball_pred':>12}  {'α_ball':>8}  {'fiable':>7}")
    print("-"*65)

    N_arr=[]; TB_m=[]; TB_s=[]; ball_m=[]; ball_s=[]; fiable=[]
    for N, rows in sorted(datos_N.items()):
        hi = [r for r in rows if r["v"] > V_CUT]
        if not hi: continue
        # Para N=7, solo usar runs donde α_ball > -0.05 (filtrar ruido)
        if N == 7:
            hi_ok = [r for r in hi if r["a_ball"] > -0.05]
            ok = len(hi_ok) >= 2
            hi_use = hi_ok if ok else hi
        else:
            hi_use = hi; ok = True

        aTB_vals   = [r["a_TB"]   for r in hi_use]
        aball_vals = [r["a_ball"] for r in hi_use]
        N_arr.append(N)
        TB_m.append(np.mean(aTB_vals));   TB_s.append(np.std(aTB_vals))
        ball_m.append(np.mean(aball_vals)); ball_s.append(np.std(aball_vals))
        fiable.append(ok)

        pred_TB   = (3*N-3)/(3*N-3+3*Nc_TB)   if N>1 else 0.0
        pred_ball = (3*N-3)/(3*N-3+3*Nc_ball) if N>1 else 0.0
        print(f"{N:>5}  {pred_TB:>10.3f}  {np.mean(aTB_vals):>8.3f}  "
              f"{pred_ball:>12.3f}  {np.mean(aball_vals):>8.3f}  "
              f"{'✓' if ok else '⚠ ruido':>7}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

    N_cont = np.logspace(np.log10(1.2), np.log10(90), 300)

    for ax, ys, ys_err, Nc, ylabel, title_suffix in [
        (ax1, TB_m,   TB_s,   Nc_TB,
         "$\\alpha_{T_B} = T_{B1}/E_0$  (canal cinético interno)",
         f"solo $T_B$,  $N_c = {Nc_TB:.0f}$"),
        (ax2, ball_m, ball_s, Nc_ball,
         "$\\alpha_\\mathrm{ball} = \\alpha_{{T_B}} + \\alpha_{{\\Delta V_B}}$  (energía interna total)",
         f"$T_B + \\Delta V_B$,  $N_c = {Nc_ball:.0f} \\approx N_\\mathrm{{ball}}$"),
    ]:
        pred = (3*N_cont-3)/(3*N_cont-3+3*Nc)
        ax.plot(N_cont, pred, "k--", lw=1.5,
                label=f"$(3N-3)/(3N-3+3N_c)$,  $N_c={Nc:.0f}$")

        for i, (n, y, ye, ok) in enumerate(zip(N_arr, ys, ys_err, fiable)):
            ms  = 90 if ok else 70
            ec  = "k" if ok else "gray"
            mk  = "o" if ok else "^"
            lw  = 0.7
            ax.errorbar(n, y, yerr=ye, fmt="none", ecolor="gray", capsize=3, lw=1.0)
            ax.scatter(n, y, s=ms, c=[COLORES[i]], marker=mk,
                       edgecolors=ec, linewidths=lw, zorder=5)
            dx, dy = {1: (-20,5), 2:(5,-14), 7:(5,6), 13:(5,5), 55:(5,5), 83:(5,5)}.get(n,(5,5))
            ax.annotate(f"$N={n}${'*' if not ok else ''}", (n, y),
                        xytext=(dx, dy), textcoords="offset points", fontsize=8.5,
                        color="gray" if not ok else "black")

        ax.set_xscale("log")
        ax.set_xlim(0.8, 110)
        ax.set_xticks([1, 2, 7, 13, 55, 83])
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}" if x >= 1 else ""))
        ax.set_xlabel("Número de átomos $N$")
        ax.set_ylabel(ylabel)
        ax.set_ylim(-0.02, 0.62)
        ax.legend(loc="upper left", fontsize=9)

    if any(not ok for ok in fiable):
        ax2.text(0.02, 0.02, "* N=7: E₀ pequeño, relación S/N baja",
                 transform=ax2.transAxes, fontsize=8, color="gray",
                 va="bottom")

    fig.tight_layout()
    out = outdir / "fig_B3_aTB_vs_aball.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"\n  Guardado: {out}")


# ---------------------------------------------------------------------------
# Fig B_consolidada – 2×2: e₁ | α_surf | α_TB | α_ball
# ---------------------------------------------------------------------------

def fig_B_consolidada(datos_N: dict, outdir: Path):
    """
    Figura única de publicación con 4 paneles (2×2).
      (a) e₁(N)        — COR decrece con N
      (b) α_surf(N)    — máximo no-monótono en N=13
      (c) α_TB(N)      — fórmula Gutiérrez con Nc=304 (solo KE interna)
      (d) α_ball(N)    — fórmula Gutiérrez con Nc≈84≈N_ball (KE+PE interna)
    """
    V_CUT = 0.9
    OFFSETS = {1:(-22,5), 2:(5,-14), 7:(5,6), 13:(5,5), 55:(5,5), 83:(5,5)}

    # Calibrar Nc desde N=83
    hi83     = [r for r in datos_N[83] if r["v"] > V_CUT]
    Nc_TB    = (3*83-3) * (1/np.mean([r["a_TB"]   for r in hi83]) - 1) / 3
    Nc_ball  = (3*83-3) * (1/np.mean([r["a_ball"] for r in hi83]) - 1) / 3

    # Recopilar estadísticas por N
    N_arr=[]; e1_m=[]; e1_s=[]; as_m=[]; as_s=[]
    TB_m=[]; TB_s=[]; ball_m=[]; ball_s=[]
    for N, rows in sorted(datos_N.items()):
        hi = [r for r in rows if r["v"] > V_CUT]
        if not hi: continue
        N_arr.append(N)
        e1_m.append(np.mean([r["e1"]     for r in hi]));  e1_s.append(np.std([r["e1"]     for r in hi]))
        as_m.append(np.mean([r["a_surf"] for r in hi]));  as_s.append(np.std([r["a_surf"] for r in hi]))
        TB_m.append(np.mean([r["a_TB"]   for r in hi]));  TB_s.append(np.std([r["a_TB"]   for r in hi]))
        ball_m.append(np.mean([r["a_ball"] for r in hi])); ball_s.append(np.std([r["a_ball"] for r in hi]))

    N_cont = np.logspace(np.log10(1.2), np.log10(90), 300)

    def _scatter(ax, ys, ys_err, ylim, ylabel, panel_label):
        ax.errorbar(N_arr, ys, yerr=ys_err, fmt="none", ecolor="gray", capsize=3, lw=1.0)
        ax.scatter(N_arr, ys, s=80, c=COLORES[:len(N_arr)],
                   edgecolors="k", lw=0.7, zorder=5)
        for n, y in zip(N_arr, ys):
            dx, dy = OFFSETS.get(n, (5,5))
            ax.annotate(f"$N={n}$", (n, y), xytext=(dx, dy),
                        textcoords="offset points", fontsize=8)
        ax.text(0.03, 0.96, panel_label, transform=ax.transAxes,
                fontsize=12, fontweight="bold", va="top")
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_ylim(*ylim)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    (ax_e1, ax_surf), (ax_TB, ax_ball) = axes

    # (a) e₁
    _scatter(ax_e1, e1_m, e1_s, (0.28, 1.08),
             "$e_1 = \\sqrt{h_1/h_0}$  (COR primer rebote)", "(a)")
    ax_e1.axhline(1.0, color="k", ls=":", lw=0.8, alpha=0.5)

    # (b) α_surf  — anotar máximo
    _scatter(ax_surf, as_m, as_s, (-0.01, 0.48),
             "$\\alpha_\\mathrm{surf}$  (fonones superficie / $E_0$)", "(b)")
    idx13 = N_arr.index(13)
    ax_surf.annotate("máximo", (13, as_m[idx13]),
                     xytext=(22, -16), textcoords="offset points",
                     fontsize=8, color="#d62728",
                     arrowprops=dict(arrowstyle="->", color="#d62728", lw=0.9))

    # (c) α_TB + fórmula
    ax_TB.plot(N_cont, (3*N_cont-3)/(3*N_cont-3+3*Nc_TB), "k--", lw=1.5,
               label=f"$(3N-3)/(3N-3+3N_c)$\n$N_c={Nc_TB:.0f}$")
    _scatter(ax_TB, TB_m, TB_s, (-0.02, 0.32),
             "$\\alpha_{{T_B}} = T_{{B1}}/E_0$  (KE interna)", "(c)")
    ax_TB.legend(loc="upper left", fontsize=8, framealpha=0.9)
    ax_TB.set_xlabel("Número de átomos $N$", fontsize=10)

    # (d) α_ball + fórmula  — resultado principal
    ax_ball.plot(N_cont, (3*N_cont-3)/(3*N_cont-3+3*Nc_ball), "k--", lw=1.5,
                 label=f"$(3N-3)/(3N-3+3N_c)$\n$N_c={Nc_ball:.0f}\\approx N_\\mathrm{{ball}}$")
    _scatter(ax_ball, ball_m, ball_s, (-0.02, 0.62),
             "$\\alpha_\\mathrm{ball} = \\alpha_{{T_B}}+\\alpha_{{\\Delta V_B}}$  (KE+PE interna)", "(d)")
    ax_ball.legend(loc="upper left", fontsize=8, framealpha=0.9)
    ax_ball.set_xlabel("Número de átomos $N$", fontsize=10)
    # N=7 tiene α_ball < 0 (S/R bajo): marcar con triángulo en el borde inferior
    idx7 = N_arr.index(7)
    ax_ball.annotate("", xy=(7, -0.018), xytext=(7, 0.03),
                     arrowprops=dict(arrowstyle="-|>", color="gray", lw=1.2))
    ax_ball.text(7, 0.032, "$N=7$\n(S/R$<$0)", ha="left", va="bottom",
                 fontsize=7.5, color="gray")

    # Eje x: log scale, ticks en los N medidos
    fmt = plt.FuncFormatter(lambda x, _: f"{int(x)}" if x >= 1 else "")
    for ax in axes.flat:
        ax.set_xscale("log")
        ax.set_xlim(0.8, 130)
        ax.set_xticks([1, 2, 7, 13, 55, 83])
        ax.xaxis.set_major_formatter(fmt)

    fig.tight_layout(h_pad=0.8, w_pad=1.2)
    out = outdir / "fig_B_consolidada.pdf"
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
    fig_aTB_vs_aball(datos_N, FIGDIR)
    fig_B_consolidada(datos_N, FIGDIR)
    print(f"\nFiguras guardadas en: {FIGDIR}")


if __name__ == "__main__":
    main()
