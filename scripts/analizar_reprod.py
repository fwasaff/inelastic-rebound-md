"""
Análisis de la reproducción de Gutiérrez et al. 2012.

Genera:
  Fig. R1 – Alturas por rebote h_n vs n (reproduce Fig. 7 del 2012)
  Fig. R2 – Evolución energética K_B, T_B, V_B, U_c, U_F (reproduce Fig. 9)
  Fig. R3 – Secuencia e_n = COR por rebote (observable ausente en 2012)

Uso:
  python3 scripts/analizar_reprod.py --outdir data/outputs/reprod_2012
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
from scipy.signal import find_peaks


# ---------------------------------------------------------------------------
# Carga de datos
# ---------------------------------------------------------------------------

COLS = ["time", "Zcm", "Vcmz", "KB", "TB", "VB", "Uc", "UF", "KB3d",
        "KS", "PE_surf", "TS"]

def cargar_datos(path: Path) -> dict:
    """Lee bounce.dat y devuelve dict con arrays por columna."""
    data = np.loadtxt(path, comments="#")
    n_cols = data.shape[1]
    cols = COLS[:n_cols]
    return {col: data[:, i] for i, col in enumerate(cols)}


def encontrar_rebotes(d: dict) -> tuple[list[dict], float]:
    """
    Detecta eventos de rebote a partir de mínimos locales de Zcm.

    Las alturas h_n se miden sobre el nivel de contacto efectivo z0, definido
    como la media de los Zcm en los mínimos (profundidad de penetración promedio).
    El COR se calcula a partir de alturas: e_n = sqrt(h_{n+1} / h_n).

    Returns
    -------
    rebotes : list[dict]
        Cada dict tiene: idx, time, h_n (altura apogeo post-rebote sobre z0), e_n
    z0 : float
        Nivel de contacto efectivo en coordenadas de caja [sigma]
    """
    Zcm  = d["Zcm"]
    time = d["time"]

    min_idx, _ = find_peaks(-Zcm, prominence=0.5, distance=50)
    max_idx, _ = find_peaks( Zcm, prominence=0.5, distance=50)

    if len(min_idx) == 0:
        return [], float(Zcm.min())

    # Nivel de contacto efectivo: media de los mínimos (ajusta por deformación surface)
    z0 = float(np.mean(Zcm[min_idx]))

    # Construir lista de apogeos ENTRE contactos consecutivos
    # apogeo[n] = max(Zcm) entre contacto n y n+1
    apogeos = []
    for k, idx in enumerate(min_idx):
        next_idx = min_idx[k + 1] if k + 1 < len(min_idx) else len(Zcm) - 1
        apex_idx = idx + np.argmax(Zcm[idx:next_idx])
        apogeos.append({
            "idx":  int(apex_idx),
            "time": float(time[apex_idx]),
            "Zcm_apex": float(Zcm[apex_idx]),
            "h_n":  float(Zcm[apex_idx] - z0),   # altura sobre nivel contacto
        })

    # Calcular e_n = sqrt(h_{n+1} / h_n) entre apogeos consecutivos
    rebotes = []
    for k, ap in enumerate(apogeos):
        h_n      = ap["h_n"]
        h_n_prev = float(Zcm[0] - z0) if k == 0 else apogeos[k - 1]["h_n"]
        # COR basado en alturas; puede ser >1 si energía interna → CM
        e_n = float(np.sqrt(abs(h_n / h_n_prev))) if h_n_prev > 0.01 else float("nan")
        rebotes.append({
            "idx":    int(min_idx[k]),
            "time":   float(time[min_idx[k]]),
            "h_n":    h_n,
            "e_n":    e_n,
        })

    return rebotes, z0


# ---------------------------------------------------------------------------
# Figura R1: alturas por rebote (reproduce Fig. 7 de Gutiérrez 2012)
# ---------------------------------------------------------------------------

def fig_alturas(resultados: dict, outdir: Path):
    """
    resultados: {h0: {"data": d, "rebotes": rebotes, "z0": z0}}
    """
    fig, ax = plt.subplots(figsize=(7, 5))

    colores = list(mcolors.TABLEAU_COLORS.values())
    markers = ["o", "s", "^", "D", "v", "P"]

    for i, (h0, res) in enumerate(sorted(resultados.items())):
        rebotes = res["rebotes"]
        z0      = res["z0"]
        d       = res["data"]
        if not rebotes:
            continue

        # n=0: altura inicial sobre nivel de contacto
        h_init  = float(d["Zcm"][0] - z0)
        alturas = [h_init] + [r["h_n"] for r in rebotes]
        bounces = list(range(len(alturas)))

        ax.plot(bounces, alturas,
                color=colores[i % len(colores)],
                marker=markers[i % len(markers)],
                markersize=6, linewidth=1.5,
                label=f"$h_0 = {h0}\\,\\sigma$")

    ax.set_xlabel("Número de rebote $n$", fontsize=12)
    ax.set_ylabel("Altura $h_n$ sobre superficie ($\\sigma$)", fontsize=12)
    ax.set_title("Reproducción Fig. 7 — Gutiérrez et al. 2012", fontsize=12)
    ax.legend(fontsize=9, loc="upper right")
    ax.set_xlim(-0.2, None)
    ax.set_ylim(0, None)
    ax.grid(True, alpha=0.3)

    out = outdir / "fig_R1_alturas.pdf"
    fig.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Guardado: {out}")


# ---------------------------------------------------------------------------
# Figura R2: evolución energética (reproduce Fig. 9 de Gutiérrez 2012)
# ---------------------------------------------------------------------------

def fig_energia(h0: float, d: dict, rebotes: list, z0: float, outdir: Path):
    """
    Evolución temporal de K_B, T_B, ΔV_B, U_c y ΔU_F para un h0.
    Dos paneles con escalas independientes:
      - Superior: canales dinámicos (KB, TB, ΔVB, Uc) — escala pequeña
      - Inferior: ΔUF gravitatorio — escala grande
    """
    time = d["time"]
    KB   = d["KB3d"]
    TB   = d["TB"]
    VB   = d["VB"]
    Uc   = d["Uc"]
    UF   = d["UF"]

    dVB = VB - VB[0]    # cambio PE intra-bola (anarmónico)
    dUF = UF - UF[0]    # cambio PE gravitatoria

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True,
                                    gridspec_kw={"height_ratios": [2, 1]})

    # --- Panel superior: canales cinéticos e internos ---
    ax1.plot(time, KB,  label="$K_B$ (traslacional CM)", color="black",  lw=1.5)
    ax1.plot(time, TB,  label="$T_B$ (KE interna bola)", color="green",  lw=1.0)
    ax1.plot(time, dVB, label="$\\Delta V_B$ (PE intra-bola)", color="red",   lw=1.0)
    ax1.plot(time, Uc,  label="$U_c$ (colisional)",      color="blue",   lw=0.8, alpha=0.8)

    # Marcadores de rebote — evitar solapamiento para runs largos
    n_rebotes  = len(rebotes)
    label_step = max(1, n_rebotes // 8)
    for n, rb in enumerate(rebotes, start=1):
        ax1.axvline(rb["time"], color="gray", lw=0.5, ls=":")
        ax2.axvline(rb["time"], color="gray", lw=0.5, ls=":")

    # Etiquetas numéricas de rebote: solo las que no se superponen
    ymax1 = max(KB.max(), TB.max(), max(dVB.max(), 0), Uc.max()) * 1.08
    ymin1 = min(KB.min(), TB.min(), dVB.min(), Uc.min(), 0) * 1.05
    ax1.set_ylim(ymin1, ymax1 * 1.12)
    for n, rb in enumerate(rebotes, start=1):
        if n == 1 or n % label_step == 0 or n == n_rebotes:
            ax1.text(rb["time"], ymax1 * 1.05,
                     f"{n}", ha="center", fontsize=7, color="gray")

    ax1.set_ylabel("Energía ($\\varepsilon$)", fontsize=11)
    ax1.set_title(f"Reproducción Fig. 9 — $h_0 = {h0}\\,\\sigma$", fontsize=12)
    ax1.legend(fontsize=8, loc="upper right", ncol=2)
    ax1.grid(True, alpha=0.3)

    # --- Panel inferior: energía gravitatoria ---
    ax2.plot(time, dUF, label="$\\Delta U_F$ (gravitatorio)", color="purple", lw=1.0, ls="--")
    ax2.axhline(0, color="k", lw=0.5, alpha=0.4)
    ax2.set_xlabel("Tiempo ($\\tau$)", fontsize=11)
    ax2.set_ylabel("$\\Delta U_F$ ($\\varepsilon$)", fontsize=11)
    ax2.legend(fontsize=8, loc="lower right")
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    out = outdir / f"fig_R2_energia_h0_{h0}.pdf"
    fig.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Guardado: {out}")


# ---------------------------------------------------------------------------
# Figura R3: secuencia e_n (observable nuevo, ausente en el 2012)
# ---------------------------------------------------------------------------

def fig_cor_secuencia(resultados: dict, outdir: Path):
    """Secuencia e_n vs n para cada h0."""
    fig, ax = plt.subplots(figsize=(7, 5))

    colores = list(mcolors.TABLEAU_COLORS.values())
    markers = ["o", "s", "^", "D", "v", "P"]

    for i, (h0, res) in enumerate(sorted(resultados.items())):
        rebotes = res["rebotes"]
        if not rebotes:
            continue

        en = [r["e_n"] for r in rebotes if not np.isnan(r["e_n"])]
        ns = list(range(1, len(en) + 1))

        ax.plot(ns, en,
                color=colores[i % len(colores)],
                marker=markers[i % len(markers)],
                markersize=7, linewidth=1.5,
                label=f"$h_0 = {h0}\\,\\sigma$")

    ax.set_xlabel("Número de rebote $n$", fontsize=12)
    ax.set_ylabel("Coeficiente de restitución $e_n = \\sqrt{h_{n}/h_{n-1}}$", fontsize=12)
    ax.set_title("Secuencia $\\{e_n\\}$ — observable no reportado en 2012", fontsize=12)
    ax.legend(fontsize=9, loc="upper right")
    ax.set_xlim(0.5, None)
    ax.set_ylim(0, None)
    ax.axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.5, label="$e_n=1$")
    ax.grid(True, alpha=0.3)

    # e_n > 1 es físicamente válido: energía interna → KE traslacional
    ax.text(0.02, 0.05,
            "$e_n>1$: energía interna de la bola → KE traslacional CM",
            transform=ax.transAxes, fontsize=8, color="gray")

    out = outdir / "fig_R3_COR_secuencia.pdf"
    fig.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Guardado: {out}")


# ---------------------------------------------------------------------------
# Tabla resumen de rebotes
# ---------------------------------------------------------------------------

def imprimir_tabla_rebotes(resultados: dict):
    print("\n" + "=" * 65)
    print(f"{'h0_param':>9} {'h0_eff':>7} {'n':>4} {'h_n (σ)':>10} {'e_n':>8}")
    print("-" * 65)
    for h0, res in sorted(resultados.items()):
        z0     = res["z0"]
        h_init = float(res["data"]["Zcm"][0] - z0)
        for n, rb in enumerate(res["rebotes"], start=1):
            print(f"{h0:>9.2f} {h_init:>7.2f} {n:>4d} {rb['h_n']:>10.4f} {rb['e_n']:>8.4f}")
        print()
    print("=" * 65)
    print("h0_param: parámetro LAMMPS (desde cima box)")
    print("h0_eff  : altura efectiva inicial sobre nivel de contacto z0")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Análisis de la reproducción de Gutiérrez et al. 2012"
    )
    parser.add_argument("--outdir", default="data/outputs/reprod_2012",
                        help="Directorio con los resultados de las simulaciones")
    parser.add_argument("--h0-fig9", type=float, default=11.03,
                        help="Altura h0 para la figura de energía (Fig. R2)")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    if not outdir.exists():
        print(f"ERROR: directorio no encontrado: {outdir}")
        return

    # --- Cargar todos los resultados -----------------------------------------
    H0_VALUES = [11.03, 18.64, 28.15, 41.47, 58.6, 79.53]
    resultados = {}

    for h0 in H0_VALUES:
        dat_file = outdir / f"h0_{h0}" / "bounce.dat"
        if not dat_file.exists():
            print(f"  [AVISO] No encontrado: {dat_file}")
            continue

        print(f"  Cargando h0={h0}...")
        d = cargar_datos(dat_file)
        rebotes, z0 = encontrar_rebotes(d)
        resultados[h0] = {"data": d, "rebotes": rebotes, "z0": z0}
        h_init = float(d["Zcm"][0] - z0)
        print(f"    → {len(rebotes)} rebotes detectados, z0={z0:.3f}σ, h_init={h_init:.2f}σ")

    if not resultados:
        print("No se encontraron datos. Correr primero: bash scripts/run_reprod.sh")
        return

    # --- Generar figuras -----------------------------------------------------
    print("\nGenerando figuras...")

    fig_alturas(resultados, outdir)

    for h0, res in resultados.items():
        fig_energia(h0, res["data"], res["rebotes"], res["z0"], outdir)

    fig_cor_secuencia(resultados, outdir)

    imprimir_tabla_rebotes(resultados)

    print(f"\nFiguras guardadas en: {outdir}")


if __name__ == "__main__":
    main()
