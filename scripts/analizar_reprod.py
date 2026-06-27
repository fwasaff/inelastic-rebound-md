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

COLS = ["time", "Zcm", "Vcmz", "KB", "TB", "VB", "Uc", "UF", "KB3d"]

def cargar_datos(path: Path) -> dict:
    """Lee bounce_h0_XX.dat y devuelve dict con arrays por columna."""
    data = np.loadtxt(path, comments="#")
    return {col: data[:, i] for i, col in enumerate(COLS)}


def encontrar_rebotes(d: dict, z_surf_top: float = 9.3,
                      ventana: int = 5) -> list[dict]:
    """
    Detecta eventos de rebote como mínimos locales de Zcm.
    Devuelve lista de dicts con: idx, time, h_antes, h_despues, Vcmz_antes, Vcmz_despues
    """
    Zcm  = d["Zcm"]
    Vcmz = d["Vcmz"]
    time = d["time"]

    # Mínimos de Zcm: el CM de la bola toca la superficie
    # La distancia de contacto es aproximadamente z_surf_top + rball
    # Usamos find_peaks sobre -Zcm para encontrar mínimos
    min_idx, _ = find_peaks(-Zcm, prominence=0.5, distance=50)

    rebotes = []
    for idx in min_idx:
        # Ventana antes y después del mínimo
        i0 = max(0, idx - ventana)
        i1 = min(len(Zcm) - 1, idx + ventana)

        # Velocidad antes (máxima en módulo, descendente) y después (ascendente)
        vz_antes  = Vcmz[i0]   # negativo (bajando)
        vz_despues = Vcmz[i1]  # positivo (subiendo)

        if vz_antes >= 0 or vz_despues <= 0:
            continue  # no es un rebote real

        # Altura antes: máximo de Zcm en la fase de caída previa
        if rebotes:
            last_idx = rebotes[-1]["idx"]
            h_antes = np.max(Zcm[last_idx:idx])
        else:
            h_antes = np.max(Zcm[:idx])

        # Altura después: máximo de Zcm en la fase de subida siguiente
        next_idx = min_idx[list(min_idx).index(idx) + 1] \
                   if idx != min_idx[-1] else len(Zcm) - 1
        h_despues = np.max(Zcm[idx:next_idx])

        cor = abs(vz_despues) / abs(vz_antes)

        rebotes.append({
            "idx":       idx,
            "time":      time[idx],
            "h_antes":   h_antes,
            "h_despues": h_despues,
            "Vcmz_antes":  vz_antes,
            "Vcmz_despues": vz_despues,
            "e_n":       cor,
        })

    return rebotes


# ---------------------------------------------------------------------------
# Figura R1: alturas por rebote (reproduce Fig. 7 de Gutiérrez 2012)
# ---------------------------------------------------------------------------

def fig_alturas(resultados: dict, outdir: Path):
    """
    resultados: {h0: {"data": d, "rebotes": rebotes}}
    """
    fig, ax = plt.subplots(figsize=(7, 5))

    colores = list(mcolors.TABLEAU_COLORS.values())
    markers = ["o", "s", "^", "D", "v", "P"]

    for i, (h0, res) in enumerate(sorted(resultados.items())):
        rebotes = res["rebotes"]
        if not rebotes:
            continue

        # Bounce 0 = altura inicial
        alturas = [h0] + [r["h_despues"] for r in rebotes]
        bounces = list(range(len(alturas)))

        ax.plot(bounces, alturas,
                color=colores[i % len(colores)],
                marker=markers[i % len(markers)],
                markersize=6, linewidth=1.5,
                label=f"$h_0 = {h0}\\,\\sigma$")

    ax.set_xlabel("Número de rebote $n$", fontsize=12)
    ax.set_ylabel("Altura máxima $h_n$ ($\\sigma$)", fontsize=12)
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

def fig_energia(h0: float, d: dict, rebotes: list, outdir: Path):
    """Evolución temporal de K_B, T_B, V_B, U_c, U_F para un h0."""
    time = d["time"]
    KB   = d["KB3d"]   # KE traslacional 3D del CM
    TB   = d["TB"]     # KE interna
    VB   = d["VB"]     # PE intra-bola
    Uc   = d["Uc"]     # PE colisional
    UF   = d["UF"]     # PE gravitatoria

    fig, ax = plt.subplots(figsize=(10, 5))

    # Escalar UF para que quepa en el mismo gráfico (como en Fig. 9 del 2012)
    ax.plot(time, KB,      label="$K_B$ (traslacional CM)", color="black",   lw=1.5)
    ax.plot(time, TB,      label="$T_B$ (KE interna)",      color="green",   lw=1.0)
    ax.plot(time, VB,      label="$V_B$ (PE intra-bola)",   color="red",     lw=1.0)
    ax.plot(time, Uc,      label="$U_c$ (colisional)",      color="blue",    lw=1.0)
    ax.plot(time, UF - UF[0], label="$\\Delta U_F$ (gravitatorio)", color="purple", lw=1.0, ls="--")

    # Marcar instantes de rebote
    for n, rb in enumerate(rebotes, start=1):
        ax.axvline(rb["time"], color="gray", lw=0.8, ls=":")
        ax.text(rb["time"], ax.get_ylim()[1] * 0.9,
                f"{n}°", ha="center", fontsize=8, color="gray")

    ax.set_xlabel("Tiempo ($\\tau$)", fontsize=12)
    ax.set_ylabel("Energía ($\\varepsilon$)", fontsize=12)
    ax.set_title(f"Reproducción Fig. 9 — $h_0 = {h0}\\,\\sigma$", fontsize=12)
    ax.legend(fontsize=8, loc="upper right", ncol=2)
    ax.grid(True, alpha=0.3)

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

        en = [r["e_n"] for r in rebotes]
        ns = list(range(1, len(en) + 1))

        ax.plot(ns, en,
                color=colores[i % len(colores)],
                marker=markers[i % len(markers)],
                markersize=7, linewidth=1.5,
                label=f"$h_0 = {h0}\\,\\sigma$")

        # Predicción geométrica (e constante = e_1)
        if en:
            e_geom = [en[0]] * len(ns)
            ax.plot(ns, e_geom, ls="--",
                    color=colores[i % len(colores)], alpha=0.4, lw=1.0)

    ax.set_xlabel("Número de rebote $n$", fontsize=12)
    ax.set_ylabel("Coeficiente de restitución $e_n$", fontsize=12)
    ax.set_title("Secuencia $\\{e_n\\}$ — observable no reportado en 2012", fontsize=12)
    ax.legend(fontsize=9, loc="upper right")
    ax.set_xlim(0.5, None)
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="k", lw=0.5)

    # Nota aclaratoria
    ax.text(0.02, 0.05,
            "Líneas discontinuas: predicción geométrica con $e_n = e_1 = $ const.",
            transform=ax.transAxes, fontsize=8, color="gray")

    out = outdir / "fig_R3_COR_secuencia.pdf"
    fig.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Guardado: {out}")


# ---------------------------------------------------------------------------
# Tabla resumen de rebotes
# ---------------------------------------------------------------------------

def imprimir_tabla_rebotes(resultados: dict):
    print("\n" + "=" * 60)
    print(f"{'h0':>8} {'n_rebote':>9} {'h_n (σ)':>10} {'e_n':>8}")
    print("-" * 60)
    for h0, res in sorted(resultados.items()):
        for n, rb in enumerate(res["rebotes"], start=1):
            print(f"{h0:>8.2f} {n:>9d} {rb['h_despues']:>10.4f} {rb['e_n']:>8.4f}")
        print()
    print("=" * 60)


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
        rebotes = encontrar_rebotes(d)
        resultados[h0] = {"data": d, "rebotes": rebotes}
        print(f"    → {len(rebotes)} rebotes detectados")

    if not resultados:
        print("No se encontraron datos. Correr primero: bash scripts/run_reprod.sh")
        return

    # --- Generar figuras -----------------------------------------------------
    print("\nGenerando figuras...")

    fig_alturas(resultados, outdir)

    for h0, res in resultados.items():
        fig_energia(h0, res["data"], res["rebotes"], outdir)

    fig_cor_secuencia(resultados, outdir)

    imprimir_tabla_rebotes(resultados)

    print(f"\nFiguras guardadas en: {outdir}")


if __name__ == "__main__":
    main()
