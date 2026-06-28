"""
Análisis físico detallado de la reproducción de Gutiérrez et al. 2012.

Genera cinco figuras:
  Fig 1 – Cascada de alturas normalizadas h_n/h0 (escala log)
  Fig 2 – Partición CORRECTA de energía en el primer rebote (4 contribuciones)
  Fig 3 – e₁ vs velocidad de impacto
  Fig 4 – Acumulación de T_B/E₀ y evolución de α_ball(n) por rebote
  Fig 5 – Balance directo de energía: E_ball(t) + E_surf(t) vs tiempo

Partición correcta (3 canales, suman 1.000):
  α_CM   = M*g*h1 / E0               → fracción en CM gravitatorio
  α_surf = ΔE_surf / E0              → DIRECTA desde KS + PE_surf (fonones)
  α_ball = 1 − α_CM − α_surf        → energía interna de la bola
           = α_TB + α_VB            (TB cinético + ΔVB potencial anarmónico)

Nota sobre ΔVB: al calentarse, la bola LJ se expande (anarmónico). La PE media
por par sube (menos negativa), por lo que ΔVB > 0. Para v_imp ≳ 0.9 σ/τ,
ΔVB > TB; el canal potencial domina el calentamiento interno.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.signal import find_peaks
from pathlib import Path

H0S   = [11.03, 18.64, 28.15, 41.47, 58.6, 79.53]
M     = 83
G     = 0.026
OUTDIR = Path("data/outputs/reprod_2012")
FIGDIR = OUTDIR

COLS  = "time Zcm Vcmz KB TB VB Uc UF KB3d KS PE_surf TS".split()
COLORES = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
MARKERS = ["o", "s", "^", "D", "v", "P"]

plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "legend.fontsize": 9,
    "figure.dpi": 150,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "lines.linewidth": 1.5,
    "lines.markersize": 6,
})


# ---------------------------------------------------------------------------
# Extracción de datos
# ---------------------------------------------------------------------------

def cargar_rebotes(h0: float) -> dict:
    """
    Devuelve dict con energías, alturas y partición completa para una altura h0.

    Partición en el primer rebote (3 canales):
      a_CM   = M*g*h1 / E0
      a_surf = ΔE_surf / E0  (medido directamente)
      a_ball = 1 − a_CM − a_surf
      a_TB   = TB(ap1) / E0  (sub-canal cinético de a_ball)
      a_VB   = a_ball − a_TB (sub-canal potencial de a_ball = ΔVB/E0)
    """
    path = OUTDIR / f"h0_{h0}" / "bounce.dat"
    d = np.loadtxt(path, comments="#")
    cols = COLS[:d.shape[1]]
    arr = {c: d[:, i] for i, c in enumerate(cols)}

    Zcm = arr["Zcm"];  TB = arr["TB"];  KB3d = arr["KB3d"]
    KS = arr.get("KS", np.zeros_like(Zcm))
    PE_surf = arr.get("PE_surf", np.zeros_like(Zcm))

    mins, _ = find_peaks(-Zcm, prominence=0.5, distance=50)
    maxs, _ = find_peaks( Zcm, prominence=0.5, distance=50)

    z0     = float(Zcm[mins].mean())
    h0_eff = float(Zcm[0] - z0)
    E0     = M * G * h0_eff

    # Energía superficial inicial y en el primer apogeo
    E_surf0  = KS[0] + PE_surf[0]
    E_surf_1 = KS[maxs[0]] + PE_surf[maxs[0]]
    dE_surf  = float(E_surf_1 - E_surf0)

    h1   = float(Zcm[maxs[0]] - z0)
    TB1  = float(TB[maxs[0]])

    a_CM   = h1 / h0_eff
    a_surf = dE_surf / E0
    a_ball = 1.0 - a_CM - a_surf
    a_TB   = TB1 / E0
    a_VB   = a_ball - a_TB

    apogeos = []
    for k, m in enumerate(maxs):
        apogeos.append({
            "n":    k + 1,
            "time": float(arr["time"][m]),
            "h_n":  float(Zcm[m] - z0),
            "TB":   float(TB[m]),
            "KB3d": float(KB3d[m]),
        })

    return {
        "h0":    h0,
        "z0":    z0,
        "h0_eff": h0_eff,
        "E0":    E0,
        "v_imp": float(np.sqrt(2 * G * h0_eff)),
        "h1":    h1,
        "e1":    float(np.sqrt(h1 / h0_eff)),
        "a_CM":  a_CM,
        "a_surf": a_surf,
        "a_ball": a_ball,
        "a_TB":  a_TB,
        "a_VB":  a_VB,
        "dE_surf": dE_surf,
        "apogeos": apogeos,
        "arr":   arr,
        "mins":  mins,
        "maxs":  maxs,
    }


# ---------------------------------------------------------------------------
# Fig 1 – Cascada de alturas normalizadas
# ---------------------------------------------------------------------------

def fig_cascada(datos: list[dict], outdir: Path):
    fig, ax = plt.subplots(figsize=(7, 5))

    for i, dat in enumerate(datos):
        h0_eff = dat["h0_eff"]
        ns  = [0] + [a["n"] for a in dat["apogeos"]]
        hs  = [1.0] + [a["h_n"] / h0_eff for a in dat["apogeos"]]
        ax.semilogy(ns, hs, color=COLORES[i], marker=MARKERS[i],
                    label=f"$h_0^\\mathrm{{eff}} = {h0_eff:.1f}\\,\\sigma$")

    n_ref = np.linspace(0, 12, 200)
    for e_ref, ls in [(0.9, ":"), (0.7, "--")]:
        ax.semilogy(n_ref, e_ref**(2*n_ref), color="gray", ls=ls, lw=0.8, alpha=0.6,
                    label=f"$e={e_ref}$ (geométrico)")

    ax.set_xlabel("Número de rebote $n$")
    ax.set_ylabel("Altura normalizada $h_n / h_0$")
    ax.set_xlim(-0.3, 25); ax.set_ylim(5e-3, 2)
    ax.legend(loc="upper right", ncol=2)
    out = outdir / "fig_A1_cascada.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"  Guardado: {out}")


# ---------------------------------------------------------------------------
# Fig 2 – Partición CORRECTA de energía en el primer rebote
# ---------------------------------------------------------------------------

def fig_particion(datos: list[dict], outdir: Path):
    """
    Barras apiladas con 4 capas:
      α_CM (azul) | α_TB (rojo) | α_VB (naranja) | α_surf (gris)

    α_ball = α_TB + α_VB es la energía interna total de la bola.
    α_VB (anarmónico) domina para v > 0.9 σ/τ.
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(datos)); w = 0.55

    a_CMs   = np.array([d["a_CM"]   for d in datos])
    a_TBs   = np.array([d["a_TB"]   for d in datos])
    a_VBs   = np.array([d["a_VB"]   for d in datos])
    a_surfs = np.array([d["a_surf"] for d in datos])

    ax.bar(x, a_CMs,  w, label="$\\alpha_\\mathrm{CM}$  KE traslacional conservada",
           color="#2196F3")
    ax.bar(x, a_TBs,  w, bottom=a_CMs,
           label="$\\alpha_{T_B}$  KE interna (fonones bola)",  color="#F44336")
    ax.bar(x, a_VBs,  w, bottom=a_CMs+a_TBs,
           label="$\\alpha_{\\Delta V_B}$  PE anarmónico (expansión bola)", color="#FF9800")
    ax.bar(x, a_surfs, w, bottom=a_CMs+a_TBs+a_VBs,
           label="$\\alpha_\\mathrm{surf}$  Fonones superficie (directo)", color="#9E9E9E")

    # Etiquetas de porcentaje dentro de cada segmento
    # Segmentos: (array de valores, array de bases, color de texto)
    MIN_LABEL = 0.04   # no etiquetar segmentos más delgados que 4%
    layers = [
        (a_CMs,              np.zeros(len(datos)), "white"),
        (a_TBs,              a_CMs,                "white"),
        (a_VBs,              a_CMs + a_TBs,        "black"),
        (a_surfs,            a_CMs + a_TBs + a_VBs, "black"),
    ]
    for vals, bottoms, tcol in layers:
        for j, (v, b) in enumerate(zip(vals, bottoms)):
            if v >= MIN_LABEL:
                ax.text(j, b + v / 2, f"{v:.0%}",
                        ha="center", va="center", fontsize=7.5,
                        color=tcol, fontweight="bold")

    for j, d in enumerate(datos):
        ax.text(j, 1.02, f"$e_1={d['e1']:.2f}$",
                ha="center", va="bottom", fontsize=8.5)

    ax.set_xticks(x)
    ax.set_xticklabels([f"$h_0={d['h0_eff']:.1f}$\n$v={d['v_imp']:.2f}\\,\\sigma/\\tau$"
                        for d in datos], fontsize=8)
    ax.set_ylabel("Fracción de $E_0 = Mgh_0$")
    ax.set_ylim(0, 1.14)
    ax.axhline(1.0, color="k", lw=0.7, ls="--", alpha=0.4)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0,
              ncol=1, fontsize=8.5, framealpha=0.95, edgecolor="0.7")
    fig.tight_layout()
    out = outdir / "fig_A2_particion.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"  Guardado: {out}")


# ---------------------------------------------------------------------------
# Fig 3 – e₁ vs velocidad de impacto
# ---------------------------------------------------------------------------

def fig_e1_vs_v(datos: list[dict], outdir: Path):
    fig, ax = plt.subplots(figsize=(6, 4.5))

    v_imp = np.array([d["v_imp"] for d in datos])
    e1    = np.array([d["e1"]    for d in datos])

    ax.scatter(v_imp, e1, s=80, zorder=5, c=COLORES[:len(datos)], marker="o",
               edgecolors="k", linewidths=0.5)

    for i, dat in enumerate(datos):
        aps = dat["apogeos"]
        if len(aps) >= 2:
            e2 = float(np.sqrt(abs(aps[1]["h_n"] / aps[0]["h_n"])))
            ax.scatter(v_imp[i], e2, s=55, marker="s", color=COLORES[i], alpha=0.85,
                       edgecolors="k", linewidths=0.4)
        if len(aps) >= 3:
            e3 = float(np.sqrt(abs(aps[2]["h_n"] / aps[1]["h_n"])))
            ax.scatter(v_imp[i], e3, s=35, marker="^", color=COLORES[i], alpha=0.85,
                       edgecolors="k", linewidths=0.4)

    mask = (v_imp > 0.8) & (v_imp < 2.1)
    if mask.sum() >= 3:
        coeffs = np.polyfit(np.log(v_imp[mask]), np.log(e1[mask]), 1)
        v_fit  = np.linspace(0.6, 2.2, 100)
        e_fit  = np.exp(coeffs[1]) * v_fit**coeffs[0]
        ax.plot(v_fit, e_fit, "k--", lw=1.2, alpha=0.6,
                label=f"$e_1 \\propto v^{{{coeffs[0]:.2f}}}$")

    # Leyenda de forma: círculo=e1, cuadrado=e2, triángulo=e3
    ax.scatter([], [], s=80, marker="o", color="gray", edgecolors="k", lw=0.4,
               label="$e_1$ ($n=1$, primer rebote)")
    ax.scatter([], [], s=55, marker="s", color="gray", edgecolors="k", lw=0.4,
               label="$e_2$ ($n=2$)")
    ax.scatter([], [], s=35, marker="^", color="gray", edgecolors="k", lw=0.4,
               label="$e_3$ ($n=3$)")
    ax.axhline(1.0, color="k", lw=0.8, ls=":", alpha=0.5)
    ax.set_xlabel("Velocidad de impacto $v_0 = \\sqrt{2gh_0}$ ($\\sigma/\\tau$)")
    ax.set_ylabel("Coeficiente de restitución $e_n$")
    ax.set_xlim(0.5, 2.3); ax.set_ylim(0.0, 1.15)
    ax.legend(loc="lower left", fontsize=8.5)
    out = outdir / "fig_A3_e1_vs_v.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"  Guardado: {out}")


# ---------------------------------------------------------------------------
# Fig 4 – Acumulación de energía interna por rebote (TB y E_ball total)
# ---------------------------------------------------------------------------

def fig_TB_acumulacion(datos: list[dict], outdir: Path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    for i, dat in enumerate(datos):
        ns   = [a["n"]              for a in dat["apogeos"]]
        tbno = [a["TB"] / dat["E0"] for a in dat["apogeos"]]
        ax1.plot(ns, tbno, color=COLORES[i], marker=MARKERS[i],
                 label=f"$h_0={dat['h0_eff']:.1f}\\,\\sigma$")

    ax1.set_xlabel("Número de rebote $n$")
    ax1.set_ylabel("$T_B / E_0$ (KE interna cinética)")
    ax1.set_title("Calentamiento cinético interno de la bola")
    ax1.set_xlim(0.5, 13); ax1.set_ylim(0, 0.55)
    ax1.legend(fontsize=8)

    # Panel derecho: convergencia al plateau h* — solo h0=28.15 (run largo)
    dat_largo = next((d for d in datos if abs(d["h0"] - 28.15) < 0.1), None)
    if dat_largo and len(dat_largo["apogeos"]) > 5:
        ns_l  = [a["n"]   for a in dat_largo["apogeos"]]
        hn_l  = [a["h_n"] for a in dat_largo["apogeos"]]
        # Plateau: media de los últimos 6 apogeos
        h_star = float(np.mean(hn_l[-6:]))
        TB_plateau = float(np.mean([a["TB"] for a in dat_largo["apogeos"][-6:]]))
        # T_B del script es energía cinética interna total. La fórmula Gutiérrez usa
        # temperatura efectiva k_B*T_ball = T_B / (3N-3)/2, con N = N_ball real (79)
        N_ball_real = int(round(dat_largo["E0"] / (G * dat_largo["h0_eff"])))  # ≈ 83
        n_int_modes = 3 * N_ball_real - 3          # grados internos
        kT_ball     = TB_plateau / (n_int_modes / 2)  # temperatura efectiva [ε]
        h_pred      = kT_ball / G                  # h* = k_B T / (m*g), m=1 por átomo

        ax2.plot(ns_l, hn_l, color=COLORES[2], marker=MARKERS[2], ms=4, lw=1.0,
                 label=f"$h_0={dat_largo['h0_eff']:.1f}\\,\\sigma$")
        ax2.axhline(h_star, color="k",       lw=1.2, ls="--",
                    label=f"$h^* = {h_star:.2f}\\,\\sigma$ (simulación)")
        ax2.axhline(h_pred,  color="#E53935", lw=1.2, ls=":",
                    label=f"$k_BT_\\mathrm{{ball}}/mg = {h_pred:.2f}\\,\\sigma$ (Gutiérrez)")
        ax2.set_xlabel("Número de rebote $n$")
        ax2.set_ylabel("Altura de apogeo $h_n$ ($\\sigma$)")
        ax2.set_title("Convergencia al plateau $h^*$ — run 5000$\\,\\tau$")
        ax2.set_ylim(bottom=0)
        ax2.legend(fontsize=8, loc="upper right")
    else:
        ax2.text(0.5, 0.5, "Sin datos de run largo\n(h0=28.15, 5000τ)",
                 ha="center", va="center", transform=ax2.transAxes, fontsize=10)
    fig.tight_layout()
    out = outdir / "fig_A4_TB_acumulacion.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"  Guardado: {out}")


# ---------------------------------------------------------------------------
# Fig 5 – Balance directo de energía E_ball(t) + E_surf(t) = cte
# ---------------------------------------------------------------------------

def fig_balance_energia(datos: list[dict], outdir: Path):
    """
    Muestra E_ball(t) y E_surf(t)−E_surf(0) para verificar conservación y
    cuantificar el flujo de energía bola → superficie.
    Solo para el caso h0=28.15 (v≈1.16, más ilustrativo).
    """
    # Elegir h0=28.15 como caso representativo
    dat = next((d for d in datos if abs(d["h0"] - 28.15) < 0.1), datos[2])
    arr  = dat["arr"]
    maxs = dat["maxs"]

    time    = arr["time"]
    KB3d    = arr["KB3d"]
    TB      = arr["TB"]
    VB      = arr["VB"]
    Uc      = arr["Uc"]
    UF      = arr["UF"]
    KS      = arr.get("KS", np.zeros_like(time))
    PE_surf = arr.get("PE_surf", np.zeros_like(time))

    # Energía total bola (con gravedad)
    E_ball = KB3d + TB + VB + 0.5*Uc + UF
    # Energía total superficie
    E_surf = KS + PE_surf

    # Normalizar: mostrar cambio respecto a t=0
    dE_ball = E_ball - E_ball[0]
    dE_surf = E_surf - E_surf[0]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    # Panel superior: componentes de E_ball
    ax1.plot(time, KB3d, label="$K_B$ (CM traslacional)", color="black", lw=1.2)
    ax1.plot(time, TB,   label="$T_B$ (KE interna bola)", color="#F44336", lw=1.2)
    ax1.plot(time, VB - VB[0], label="$\\Delta V_B$ (PE interna bola)", color="#FF9800",
             lw=1.2)
    ax1.plot(time, UF - UF[0], label="$\\Delta U_F$ (gravitatorio)", color="purple",
             lw=1.2, ls="--")

    for m in maxs[:4]:
        ax1.axvline(time[m], color="gray", lw=0.6, ls=":")
    ax1.set_ylabel("Energía ($\\varepsilon$)")
    ax1.set_title(f"Balance energético — $h_0={dat['h0_eff']:.1f}\\,\\sigma$, "
                  f"$v_0={dat['v_imp']:.2f}\\,\\sigma/\\tau$")
    ax1.legend(fontsize=8.5, ncol=2, loc="upper right")

    # Panel inferior: flujo energético y verificación NVE
    ax2.plot(time, dE_surf, color="#4CAF50", lw=1.5,
             label="$\\Delta E_\\mathrm{surf}$ (fonones superficie)")
    # Verificación NVE: dE_ball + dE_surf ≈ 0 si la energía se conserva
    nve_check = dE_ball + dE_surf
    ax2.plot(time, nve_check, color="k", lw=0.9, ls="--", alpha=0.7,
             label="$\\delta E_\\mathrm{total} = \\Delta E_\\mathrm{ball} + \\Delta E_\\mathrm{surf}$"
                   " (verif. NVE, $\\approx 0$)")
    for m in maxs[:4]:
        ax2.axvline(time[m], color="gray", lw=0.6, ls=":")

    ax2.axhline(0, color="k", lw=0.6, alpha=0.5)
    ax2.set_xlabel("Tiempo ($\\tau$)")
    ax2.set_ylabel("Cambio de energía ($\\varepsilon$)")
    ax2.set_title(f"Flujo bola → superficie "
                  f"($\\alpha_\\mathrm{{surf}} = {dat['a_surf']:.3f}$ al primer apogeo)"
                  f" — verificación NVE")
    ax2.legend(fontsize=8.5)

    fig.tight_layout()
    out = outdir / "fig_A5_balance_energia.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"  Guardado: {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Cargando datos...")
    datos = [cargar_rebotes(h0) for h0 in H0S]

    print("\n=== PARTICIÓN DE ENERGÍA CORRECTA — PRIMER REBOTE ===")
    print(f"{'h0_eff':>8} {'v_imp':>6} {'e1':>6} | "
          f"{'α_CM':>6} {'α_TB':>6} {'α_VB':>6} {'α_surf':>7} | suma")
    for d in datos:
        s = d["a_CM"] + d["a_TB"] + d["a_VB"] + d["a_surf"]
        print(f"{d['h0_eff']:8.2f} {d['v_imp']:6.3f} {d['e1']:6.3f} | "
              f"{d['a_CM']:6.3f} {d['a_TB']:6.3f} {d['a_VB']:6.3f} {d['a_surf']:7.3f} | "
              f"{s:.3f}")

    print("\n  Nota: α_ball = α_TB + α_VB = energía interna TOTAL de la bola")
    print("  α_VB domina para v > 0.9 σ/τ: expansión anarmónica del cluster LJ")

    print("\nGenerando figuras...")
    fig_cascada(datos, FIGDIR)
    fig_particion(datos, FIGDIR)
    fig_e1_vs_v(datos, FIGDIR)
    fig_TB_acumulacion(datos, FIGDIR)
    fig_balance_energia(datos, FIGDIR)
    print(f"\nFiguras guardadas en: {FIGDIR}")


if __name__ == "__main__":
    main()
