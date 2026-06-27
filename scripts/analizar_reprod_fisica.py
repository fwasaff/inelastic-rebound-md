"""
Análisis físico detallado de la reproducción de Gutiérrez et al. 2012.

Genera cuatro figuras centradas en la física del impacto:
  Fig 1 – Cascada de alturas normalizadas h_n/h0 (escala log)
  Fig 2 – Partición de energía en el primer rebote (3 canales)
  Fig 3 – e₁ vs velocidad de impacto (ley empírica central)
  Fig 4 – Acumulación de energía interna T_B/E₀ por rebote

Las figuras están diseñadas para comparación directa con el estudio extendido
(N variable, T*>0, c_ps>0) que viene a continuación.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.signal import find_peaks
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuración
# ---------------------------------------------------------------------------
H0S   = [11.03, 18.64, 28.15, 41.47, 58.6, 79.53]
M     = 83        # átomos en la bola
G     = 0.026     # gravedad [sigma/tau^2]
OUTDIR = Path("data/outputs/reprod_2012")
FIGDIR = OUTDIR

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
# Extracción de datos por simulación
# ---------------------------------------------------------------------------

def cargar_rebotes(h0: float) -> dict:
    """
    Devuelve dict con:
      z0, h0_eff, v_imp0, apogeos:[{n, h_n, TB, KB3d, time}], ...
    """
    path = OUTDIR / f"h0_{h0}" / "bounce.dat"
    d    = np.loadtxt(path, comments="#")
    time, Zcm, Vcmz, KB, TB, VB, Uc, UF, KB3d = d.T

    mins, _ = find_peaks(-Zcm, prominence=0.5, distance=50)
    maxs, _ = find_peaks( Zcm, prominence=0.5, distance=50)

    z0     = float(Zcm[mins].mean())
    h0_eff = float(Zcm[0] - z0)
    E0     = M * G * h0_eff        # energía gravitatoria inicial

    apogeos = []
    for k, m in enumerate(maxs):
        h_n = float(Zcm[m] - z0)
        apogeos.append({
            "n":    k + 1,
            "time": float(time[m]),
            "h_n":  h_n,
            "TB":   float(TB[m]),
            "KB3d": float(KB3d[m]),
        })

    # Partición de energía en el PRIMER rebote
    TB1   = TB[maxs[0]] if len(maxs) else 0.0
    h1    = float(Zcm[maxs[0]] - z0) if len(maxs) else 0.0
    a_CM  = h1 / h0_eff
    a_int = TB1 / E0
    a_surf = max(0.0, 1.0 - a_CM - a_int)

    return {
        "h0":    h0,
        "z0":    z0,
        "h0_eff": h0_eff,
        "E0":    E0,
        "v_imp": float(np.sqrt(2 * G * h0_eff)),
        "e1":    float(np.sqrt(h1 / h0_eff)) if h0_eff > 0 else 0.0,
        "a_CM":  a_CM,
        "a_int": a_int,
        "a_surf": a_surf,
        "apogeos": apogeos,
    }


# ---------------------------------------------------------------------------
# Fig 1 – Cascada de alturas normalizadas h_n/h0  (escala log)
# ---------------------------------------------------------------------------

def fig_cascada(datos: list[dict], outdir: Path):
    fig, ax = plt.subplots(figsize=(7, 5))

    for i, dat in enumerate(datos):
        h0_eff = dat["h0_eff"]
        ns     = [0] + [a["n"] for a in dat["apogeos"]]
        hs_norm = [1.0] + [a["h_n"] / h0_eff for a in dat["apogeos"]]

        ax.semilogy(
            ns, hs_norm,
            color=COLORES[i], marker=MARKERS[i],
            label=f"$h_0^\\mathrm{{eff}} = {h0_eff:.1f}\\,\\sigma$",
        )

    # Líneas de referencia: e=const (decaimiento geométrico puro)
    n_ref = np.linspace(0, 12, 100)
    for e_ref, ls in [(0.9, ":"), (0.7, "--")]:
        ax.semilogy(n_ref, e_ref ** (2 * n_ref),
                    color="gray", ls=ls, lw=0.8, alpha=0.6,
                    label=f"$e_n = {e_ref}$ (geométrico)")

    ax.set_xlabel("Número de rebote $n$")
    ax.set_ylabel("Altura normalizada $h_n / h_0$")
    ax.set_title("Cascada de alturas — reproducción Gutiérrez et al. 2012")
    ax.set_xlim(-0.3, 13)
    ax.set_ylim(5e-3, 2)
    ax.legend(loc="upper right", ncol=2)

    out = outdir / "fig_A1_cascada.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Guardado: {out}")


# ---------------------------------------------------------------------------
# Fig 2 – Partición de energía en el primer rebote (stacked bar)
# ---------------------------------------------------------------------------

def fig_particion(datos: list[dict], outdir: Path):
    fig, ax = plt.subplots(figsize=(7, 4.5))

    v_imps  = [d["v_imp"]  for d in datos]
    a_CMs   = [d["a_CM"]   for d in datos]
    a_ints  = [d["a_int"]  for d in datos]
    a_surfs = [d["a_surf"] for d in datos]
    h0_effs = [d["h0_eff"] for d in datos]

    x = np.arange(len(datos))
    w = 0.55

    b0 = ax.bar(x, a_CMs,   w, label="$K_B$ (CM conservado)",      color="#2196F3")
    b1 = ax.bar(x, a_ints,  w, bottom=a_CMs,
                label="$\\Delta T_B$ (modos internos bola)",       color="#F44336")
    b2 = ax.bar(x, a_surfs, w, bottom=[c+i for c,i in zip(a_CMs, a_ints)],
                label="Fonones superficie",                        color="#9E9E9E")

    # Anotar el valor de e1 sobre cada barra
    for j, d in enumerate(datos):
        ax.text(j, 1.02, f"$e_1={d['e1']:.2f}$",
                ha="center", va="bottom", fontsize=8.5)

    ax.set_xticks(x)
    ax.set_xticklabels([f"$h_0={d['h0_eff']:.1f}$\n$v={d['v_imp']:.2f}$"
                        for d in datos], fontsize=8.5)
    ax.set_xlabel("Altura inicial $h_0^\\mathrm{eff}$ (σ)  /  velocidad $v_0$ (σ/τ)")
    ax.set_ylabel("Fracción de $E_\\mathrm{grav,0} = Mgh_0$")
    ax.set_title("Partición de energía en el primer rebote")
    ax.set_ylim(0, 1.12)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22),
              ncol=3, frameon=True)

    # Línea horizontal en α_int ≈ 0.22
    ax.axhline(1.0, color="k", lw=0.7, ls="--", alpha=0.4)

    fig.tight_layout()
    out = outdir / "fig_A2_particion.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Guardado: {out}")


# ---------------------------------------------------------------------------
# Fig 3 – e₁ vs velocidad de impacto  (ley empírica)
# ---------------------------------------------------------------------------

def fig_e1_vs_v(datos: list[dict], outdir: Path):
    fig, ax = plt.subplots(figsize=(6, 4.5))

    v_imp = np.array([d["v_imp"] for d in datos])
    e1    = np.array([d["e1"]    for d in datos])

    # Scatter e1 (primer rebote)
    ax.scatter(v_imp, e1, s=80, zorder=5,
               c=COLORES[:len(datos)], marker="o",
               label="$e_1$ (primer rebote)", edgecolors="k", linewidths=0.5)

    # También mostrar e2, e3 para contexto
    for i, dat in enumerate(datos):
        aps = dat["apogeos"]
        h0e = dat["h0_eff"]
        if len(aps) >= 2:
            e2 = float(np.sqrt(abs(aps[1]["h_n"] / aps[0]["h_n"])))
            ax.scatter(v_imp[i], e2, s=40, marker="s",
                       color=COLORES[i], alpha=0.5, zorder=4)
        if len(aps) >= 3:
            e3 = float(np.sqrt(abs(aps[2]["h_n"] / aps[1]["h_n"])))
            ax.scatter(v_imp[i], e3, s=20, marker="^",
                       color=COLORES[i], alpha=0.35, zorder=3)

    # Ajuste power-law solo para los 4 puntos intermedios (excluir saturación)
    mask = (v_imp > 0.8) & (v_imp < 2.1)
    if mask.sum() >= 3:
        coeffs = np.polyfit(np.log(v_imp[mask]), np.log(e1[mask]), 1)
        v_fit  = np.linspace(0.6, 2.1, 100)
        e_fit  = np.exp(coeffs[1]) * v_fit ** coeffs[0]
        ax.plot(v_fit, e_fit, "k--", lw=1.2, alpha=0.6,
                label=f"$e_1 \\propto v^{{{coeffs[0]:.2f}}}$")

    # Dummy markers para leyenda
    ax.scatter([], [], s=80, marker="o", color="gray", label="$e_1$ (n=1)")
    ax.scatter([], [], s=40, marker="s", color="gray", alpha=0.5, label="$e_2$ (n=2)")
    ax.scatter([], [], s=20, marker="^", color="gray", alpha=0.35, label="$e_3$ (n=3)")

    ax.axhline(1.0, color="k", lw=0.8, ls=":", alpha=0.5)
    ax.set_xlabel("Velocidad de impacto $v_0 = \\sqrt{2gh_0}$ (σ/τ)")
    ax.set_ylabel("Coeficiente de restitución $e_n$")
    ax.set_title("COR vs velocidad de impacto")
    ax.set_xlim(0.5, 2.2)
    ax.set_ylim(0.0, 1.05)
    ax.legend(loc="lower left")

    out = outdir / "fig_A3_e1_vs_v.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Guardado: {out}")


# ---------------------------------------------------------------------------
# Fig 4 – Acumulación de T_B / E₀ por rebote
# ---------------------------------------------------------------------------

def fig_TB_acumulacion(datos: list[dict], outdir: Path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Panel izquierdo: T_B/E0 vs rebote n
    for i, dat in enumerate(datos):
        ns   = [a["n"]             for a in dat["apogeos"]]
        tbno = [a["TB"] / dat["E0"] for a in dat["apogeos"]]
        ax1.plot(ns, tbno, color=COLORES[i], marker=MARKERS[i],
                 label=f"$h_0={dat['h0_eff']:.1f}\\,\\sigma$")

    ax1.axhline(0.22, color="k", ls="--", lw=1.0, alpha=0.5,
                label="$\\alpha_\\mathrm{int} \\approx 0.22$")
    ax1.set_xlabel("Número de rebote $n$")
    ax1.set_ylabel("Energía interna bola $T_B / E_0$")
    ax1.set_title("Calentamiento interno de la bola")
    ax1.set_xlim(0.5, 13)
    ax1.set_ylim(0, 0.55)
    ax1.legend(fontsize=8, loc="upper left")

    # Panel derecho: h_n/h0 vs TB/E0  (espacio de fase energético)
    for i, dat in enumerate(datos):
        h0e  = dat["h0_eff"]
        hn   = [a["h_n"] / h0e   for a in dat["apogeos"]]
        tbno = [a["TB"] / dat["E0"] for a in dat["apogeos"]]
        # Trayectoria: flecha del primero al último
        ax2.plot(tbno, hn, color=COLORES[i], marker=MARKERS[i],
                 label=f"$h_0={h0e:.1f}$", alpha=0.8)
        # Punto inicial (n=1)
        if hn:
            ax2.annotate("", xy=(tbno[0], hn[0]),
                         xytext=(0, 1.0),
                         arrowprops=dict(arrowstyle="->", color=COLORES[i],
                                         lw=0.8, alpha=0.4))

    ax2.scatter([0], [1.0], s=60, c="black", zorder=5, label="n=0 (todos)")
    ax2.set_xlabel("Energía interna $T_B / E_0$")
    ax2.set_ylabel("Altura normalizada $h_n / h_0$")
    ax2.set_title("Espacio de fase $(T_B, h_n)$ — rebotes sucesivos")
    ax2.set_xlim(-0.01, 0.55)
    ax2.set_ylim(-0.02, 1.05)
    ax2.legend(fontsize=8, loc="upper right")

    fig.tight_layout()
    out = outdir / "fig_A4_TB_acumulacion.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Guardado: {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Cargando datos...")
    datos = [cargar_rebotes(h0) for h0 in H0S]

    # Tabla de partición de energía
    print("\n=== PARTICIÓN DE ENERGÍA — PRIMER REBOTE ===")
    print(f"{'h0_eff':>8} {'v_imp':>6} {'e1':>6} {'α_CM':>7} {'α_int':>7} {'α_surf':>7}")
    for d in datos:
        print(f"{d['h0_eff']:8.2f} {d['v_imp']:6.3f} {d['e1']:6.3f} "
              f"{d['a_CM']:7.3f} {d['a_int']:7.3f} {d['a_surf']:7.3f}")

    print("\nGenerando figuras...")
    fig_cascada(datos, FIGDIR)
    fig_particion(datos, FIGDIR)
    fig_e1_vs_v(datos, FIGDIR)
    fig_TB_acumulacion(datos, FIGDIR)
    print(f"\nFiguras guardadas en: {FIGDIR}")


if __name__ == "__main__":
    main()
