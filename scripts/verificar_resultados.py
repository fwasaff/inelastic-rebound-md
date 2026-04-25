"""
verificar_resultados.py
========================
Script de verificación y análisis exploratorio de los resultados del sweep
de clústeres. Produce 7 figuras de diagnóstico y 2 archivos CSV en
data/outputs/verificacion/.

Ejecución:
    python3 scripts/verificar_resultados.py

Archivos generados:
    fig1_sticking_map.pdf      — mapa de fracción de sticking por (N, v0, eps)
    fig2_balance_energia.pdf   — verificación COR² + α_int + α_surf = 1
    fig3_COR_vs_v0.pdf         — COR(v0) para cada N, T*, ε_ps (ang=0°)
    fig4_alpha_int_vs_N.pdf    — α_E_int(N) para cada v0 y ε_ps (T*=0.1, ang=0°)
    fig5_eps_comparacion.pdf   — comparación ε_ps=0.30 vs 1.00 para N=1 y N=13
    fig6_T_comparacion.pdf     — comparación T*=0.1 vs 0.3 para N=1 y N=13
    fig7_power_law_COR.pdf     — ley de potencia COR ~ A·v₀^{-β} en escala log-log
    resultados_raw.csv         — todos los runs, una fila por run (incluye We)
    resultados_resumen.csv     — medias y std por condición (ang=0°)

Cambios metodológicos implementados (basados en estado del arte, 2026-04):
    - Número de Weber: We = E_inc / (N·ε_ps)  — variable de escala para
      la transición rebote↔adhesión (Dahneke 1995, Hou et al. 2000)
    - Ley de potencia COR(v₀) ~ A·v₀^{-β}: ajuste log-log sobre N=1 (fig7)
    - Canal rotacional documentado: ke_internal incluye vibración+rotación;
      vx_cm registrado en cluster_rebound_data.dat para análisis futuro

Conceptos de Programación Avanzada aplicados (MSC 005):
    - pandas.DataFrame como estructura principal de datos (Unidad 2)
    - groupby / agg / pivot_table para análisis estadístico (Unidad 2)
    - query() para filtrado declarativo (Unidad 2)
    - exportación CSV con to_csv() (Unidad 2)
    - describe() para estadísticas descriptivas (Unidad 2)
"""

import sys
import re
import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# suprimir warnings de paquetes del sistema compilados con numpy<2
warnings.filterwarnings("ignore", message=".*NumPy.*")

# ── Path setup ────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
DATA_DIR = ROOT / "data" / "outputs"
OUT_DIR  = DATA_DIR / "verificacion"
OUT_DIR.mkdir(exist_ok=True)

from python.analysis import ClusterReboundAnalyzer

logging.basicConfig(level=logging.WARNING)

# ── Colores y estilo ──────────────────────────────────────────────────────────
COLORS_N   = {1: "#1f77b4", 2: "#ff7f0e", 7: "#2ca02c", 13: "#d62728", 55: "#9467bd"}
COLORS_EPS = {"0.30": "#e377c2", "1.00": "#17becf"}
COLORS_T   = {0.1: "#1f77b4", 0.3: "#ff7f0e", 0.5: "#2ca02c"}
MARKERS    = {1: "o", 2: "s", 7: "^", 13: "D", 55: "P"}

plt.rcParams.update({"font.size": 10, "axes.grid": True,
                      "grid.alpha": 0.3, "figure.dpi": 150})

# ═════════════════════════════════════════════════════════════════════════════
# 1. Carga de resultados → DataFrame
# ═════════════════════════════════════════════════════════════════════════════
_pattern = re.compile(
    r"cluster_N(\d+)_T([\d.]+)_v0_([\d.]+)_ang_([\d.]+)_eps([\d.]+)_rep(\d+)"
)


def load_all_results() -> pd.DataFrame:
    """Carga todos los runs eps-tagged en un DataFrame.

    Cada fila es un run individual (un directorio). Las columnas son:
        N, T, v0, ang, eps, rep          — parámetros de la simulación
        COR, alpha_int, alpha_surf, dT   — resultados físicos
        E_inc, balance, sticking         — diagnóstico energético

    Returns
    -------
    pd.DataFrame
        DataFrame con dtype explícito por columna. Vacío si no hay datos.
    """
    rows = []
    dirs = sorted(DATA_DIR.glob("cluster_N*_T*_v0*_ang*_eps*_rep*"))
    total = len(dirs)
    ok = fail = 0

    for d in dirs:
        m = _pattern.match(d.name)
        if not m:
            continue
        N, T, v0, ang, eps, rep = (
            int(m[1]), float(m[2]), float(m[3]),
            float(m[4]), m[5], int(m[6])
        )
        dat = d / "cluster_rebound_data.dat"
        if not dat.exists():
            fail += 1
            continue

        try:
            az = (ClusterReboundAnalyzer(
                    N_cluster=N, v0=v0, angle=ang,
                    T_surface_init=T, surface_top_z=9.3,
                    eps_ps=float(eps))
                  .load(dat))
            res = az.compute_all()
            rows.append({
                "N":          N,
                "T":          T,
                "v0":         v0,
                "ang":        ang,
                "eps":        eps,
                "rep":        rep,
                "COR":        res.COR,
                "alpha_int":  res.alpha_E_int,
                "alpha_surf": res.alpha_E_surf,
                "dT":         res.dT_surface,
                "E_inc":      res.E_incident,
                "We":         res.Weber_number,
                "balance":    az.energy_balance_residual(),
                "sticking":   res.COR < 0.02,
            })
            ok += 1
        except Exception as exc:
            fail += 1
            print(f"  FALLO {d.name}: {exc}")

    print(f"Cargados {ok}/{total} runs  ({fail} fallidos)")

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df = df.astype({"N": int, "rep": int, "sticking": bool})
    return df


# ═════════════════════════════════════════════════════════════════════════════
# 2. Figuras  (pandas para filtrado y agregación)
# ═════════════════════════════════════════════════════════════════════════════

def fig1_sticking_map(df: pd.DataFrame) -> None:
    """Mapa de fracción de sticking por (N, v0) para cada (T*, eps)."""
    Ns   = sorted(df["N"].unique())
    v0s  = sorted(df["v0"].unique())
    Ts   = sorted(df["T"].unique())
    epss = sorted(df["eps"].unique())

    # fracción de sticking media por (T, eps, N, v0) — una operación
    grp = (df.groupby(["T", "eps", "N", "v0"])["sticking"]
             .mean()
             .reset_index())

    fig, axes = plt.subplots(len(Ts), len(epss),
                              figsize=(5 * len(epss), 4 * len(Ts)),
                              squeeze=False)
    fig.suptitle("Fracción de sticking  (COR < 0.02)\npor (N, v0, T*, ε_ps)",
                 fontsize=12)

    for i, T in enumerate(Ts):
        for j, eps in enumerate(epss):
            ax = axes[i][j]
            sub = grp.query("T == @T and eps == @eps")
            # pivot → matriz imagen con filas=N, columnas=v0
            mat = (sub.pivot(index="N", columns="v0", values="sticking")
                      .reindex(index=Ns, columns=v0s)
                      .fillna(0)
                      .to_numpy())

            im = ax.imshow(mat, vmin=0, vmax=1, cmap="RdYlGn_r",
                           aspect="auto", origin="upper")
            ax.set_xticks(range(len(v0s)))
            ax.set_xticklabels([f"{v:.2g}" for v in v0s], fontsize=8)
            ax.set_yticks(range(len(Ns)))
            ax.set_yticklabels([f"N={N}" for N in Ns])
            ax.set_title(f"T*={T}, ε_ps={eps}")
            ax.set_xlabel("v₀  [σ/τ]")
            for ni in range(len(Ns)):
                for vi in range(len(v0s)):
                    val = mat[ni, vi]
                    ax.text(vi, ni, f"{val:.1f}", ha="center", va="center",
                            fontsize=8, color="black" if 0.2 < val < 0.8 else "white")
            fig.colorbar(im, ax=ax, label="f_sticking")

    fig.tight_layout()
    out = OUT_DIR / "fig1_sticking_map.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out.name}")


def fig2_balance_energia(df: pd.DataFrame) -> None:
    """Histograma del residuo |COR² + α_int + α_surf - 1|."""
    balances = df["balance"].to_numpy()

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(balances, bins=60, color="steelblue", edgecolor="white", linewidth=0.4)
    ax.axvline(0, color="red", lw=1.5, ls="--", label="balance perfecto")
    ax.set_xlabel("|COR² + α_int + α_surf − 1|")
    ax.set_ylabel("Número de runs")
    ax.set_title("Verificación de balance energético\n(por construcción debe ser ≡ 0)")
    ax.set_yscale("log")
    ax.legend()
    fig.suptitle(
        f"N = {len(balances)} runs    "
        f"max = {balances.max():.2e}    "
        f"mean = {balances.mean():.2e}",
        fontsize=9,
    )
    fig.tight_layout()
    out = OUT_DIR / "fig2_balance_energia.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out.name}")


def fig3_COR_vs_v0(df: pd.DataFrame) -> None:
    """COR(v0) para cada N, separado por (T*, ε_ps), ángulo=0°."""
    sub0 = df[df["ang"] == 0.0]
    Ns   = sorted(sub0["N"].unique())
    Ts   = sorted(sub0["T"].unique())
    epss = sorted(sub0["eps"].unique())

    # media y std de COR sobre réplicas — groupby con dos estadísticos
    grp = (sub0.groupby(["T", "eps", "N", "v0"])["COR"]
               .agg(mean="mean", std="std")
               .reset_index())

    fig, axes = plt.subplots(len(Ts), len(epss),
                              figsize=(6 * len(epss), 4 * len(Ts)),
                              squeeze=False, sharey=True, sharex=True)
    fig.suptitle("COR vs v₀  (ángulo = 0°)", fontsize=12)

    for i, T in enumerate(Ts):
        for j, eps in enumerate(epss):
            ax = axes[i][j]
            sub = grp.query("T == @T and eps == @eps")
            for N in Ns:
                g = sub[sub["N"] == N].sort_values("v0")
                if g.empty:
                    continue
                ax.errorbar(g["v0"], g["mean"], yerr=g["std"].fillna(0),
                            fmt=MARKERS[N] + "-", color=COLORS_N[N],
                            capsize=3, label=f"N={N}", markersize=5)
            ax.set_title(f"T*={T}, ε_ps={eps}")
            ax.set_xlabel("v₀  [σ/τ]")
            ax.set_ylabel("COR")
            ax.set_ylim(-0.05, 1.05)
            ax.legend(fontsize=8, ncol=2)

    fig.tight_layout()
    out = OUT_DIR / "fig3_COR_vs_v0.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out.name}")


def fig4_alpha_int_vs_N(df: pd.DataFrame) -> None:
    """α_E_int(N) para cada v0 y ε_ps, fijando T*=0.1 y ang=0°."""
    T_fixed, ang_fixed = 0.1, 0.0
    sub0 = df.query("T == @T_fixed and ang == @ang_fixed")
    epss = sorted(sub0["eps"].unique())
    Ns   = sorted(sub0["N"].unique())

    grp = (sub0.groupby(["eps", "v0", "N"])["alpha_int"]
               .agg(mean="mean", std="std")
               .reset_index())

    fig, axes = plt.subplots(1, len(epss),
                              figsize=(6 * len(epss), 5),
                              squeeze=False, sharey=True)
    fig.suptitle(
        f"Energía interna del clúster  α_E_int(N)\nT*={T_fixed}, ángulo={ang_fixed}°",
        fontsize=12,
    )

    for j, eps in enumerate(epss):
        ax = axes[0][j]
        sub = grp[grp["eps"] == eps]
        for v0 in sorted(sub["v0"].unique()):
            g = sub[sub["v0"] == v0].sort_values("N")
            if g.empty:
                continue
            ax.errorbar(g["N"], g["mean"], yerr=g["std"].fillna(0),
                        fmt="o-", capsize=3, label=f"v₀={v0:.1f}", markersize=5)
        ax.axhline(0, color="gray", lw=0.8, ls=":")
        ax.axhline(1, color="gray", lw=0.8, ls="--", label="α=1")
        ax.set_xscale("log")
        ax.set_xlabel("N [átomos]")
        ax.set_ylabel(r"$\alpha_{E,\mathrm{int}}$")
        ax.set_title(f"ε_ps = {eps}")
        ax.set_xticks(Ns)
        ax.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
        ax.legend(fontsize=8)

    fig.tight_layout()
    out = OUT_DIR / "fig4_alpha_int_vs_N.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out.name}")


def fig5_eps_comparacion(df: pd.DataFrame) -> None:
    """Comparación ε_ps=0.30 vs 1.00 para N=1 y N=13. T*=0.1, ang=0°."""
    T_fixed, ang_fixed = 0.1, 0.0
    Ns_plot = [N for N in [1, 13] if N in df["N"].values]

    # agg múltiples columnas y estadísticos en un solo groupby
    grp = (df.query("T == @T_fixed and ang == @ang_fixed")
             .groupby(["N", "eps", "v0"])[["COR", "alpha_int"]]
             .agg(["mean", "std"])
             .reset_index())
    grp.columns = ["N", "eps", "v0",
                   "COR_mean", "COR_std",
                   "alpha_int_mean", "alpha_int_std"]

    fig, axes = plt.subplots(len(Ns_plot), 2,
                              figsize=(10, 4 * len(Ns_plot)),
                              sharex=True)
    if len(Ns_plot) == 1:
        axes = axes.reshape(1, -1)
    fig.suptitle(
        f"Efecto de ε_ps: homoatómico (1.00) vs débil (0.30)\nT*={T_fixed}, ángulo={ang_fixed}°",
        fontsize=12,
    )

    for i, N in enumerate(Ns_plot):
        for j, (mean_col, std_col, label) in enumerate([
            ("COR_mean",       "COR_std",       "COR"),
            ("alpha_int_mean", "alpha_int_std", r"$\alpha_{E,\mathrm{int}}$"),
        ]):
            ax = axes[i][j]
            for eps in ["0.30", "1.00"]:
                g = grp[(grp["N"] == N) & (grp["eps"] == eps)].sort_values("v0")
                if g.empty:
                    continue
                ax.errorbar(g["v0"], g[mean_col], yerr=g[std_col].fillna(0),
                            fmt="o-", capsize=3, color=COLORS_EPS[eps],
                            label=f"ε_ps={eps}", markersize=5)
            ax.set_title(f"N={N}, {label}")
            ax.set_xlabel("v₀  [σ/τ]")
            ax.set_ylabel(label)
            if mean_col == "COR_mean":
                ax.set_ylim(-0.05, 1.05)
            ax.legend(fontsize=9)

    fig.tight_layout()
    out = OUT_DIR / "fig5_eps_comparacion.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out.name}")


def fig6_T_comparacion(df: pd.DataFrame) -> None:
    """Comparación T* para N=1 y N=13. ε_ps=1.00, ang=0°."""
    eps_fixed, ang_fixed = "1.00", 0.0
    Ns_plot = [N for N in [1, 13] if N in df["N"].values]
    Ts = sorted(df["T"].unique())

    grp = (df.query("eps == @eps_fixed and ang == @ang_fixed")
             .groupby(["N", "T", "v0"])[["COR", "alpha_int"]]
             .agg(["mean", "std"])
             .reset_index())
    grp.columns = ["N", "T", "v0",
                   "COR_mean", "COR_std",
                   "alpha_int_mean", "alpha_int_std"]

    fig, axes = plt.subplots(len(Ns_plot), 2,
                              figsize=(10, 4 * len(Ns_plot)),
                              sharex=True)
    if len(Ns_plot) == 1:
        axes = axes.reshape(1, -1)
    fig.suptitle(
        f"Efecto de T*: comparación de temperaturas de superficie\nε_ps={eps_fixed}, ángulo={ang_fixed}°",
        fontsize=12,
    )

    for i, N in enumerate(Ns_plot):
        for j, (mean_col, std_col, label) in enumerate([
            ("COR_mean",       "COR_std",       "COR"),
            ("alpha_int_mean", "alpha_int_std", r"$\alpha_{E,\mathrm{int}}$"),
        ]):
            ax = axes[i][j]
            for T in Ts:
                g = grp[(grp["N"] == N) & (grp["T"] == T)].sort_values("v0")
                if g.empty:
                    continue
                ax.errorbar(g["v0"], g[mean_col], yerr=g[std_col].fillna(0),
                            fmt="o-", capsize=3, color=COLORS_T[T],
                            label=f"T*={T}", markersize=5)
            ax.set_title(f"N={N}, {label}")
            ax.set_xlabel("v₀  [σ/τ]")
            ax.set_ylabel(label)
            if mean_col == "COR_mean":
                ax.set_ylim(-0.05, 1.05)
            ax.legend(fontsize=9)

    fig.tight_layout()
    out = OUT_DIR / "fig6_T_comparacion.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out.name}")


def fig7_power_law_COR(df: pd.DataFrame) -> None:
    """Ley de potencia COR(v₀) ~ A·v₀^{-β} en escala log-log (ang=0°).

    El ajuste se realiza por regresión log-log:
        ln(COR) = ln(A) − β · ln(v₀)

    Solo se incluyen puntos con COR > 0.02 (excluir casos de sticking
    donde la ley de potencia no aplica).  Se genera un subgráfico por
    cada combinación (T*, ε_ps); dentro de cada panel, una curva por N.
    """
    sub0 = df.query("ang == 0.0 and not sticking and COR > 0.0")
    if sub0.empty:
        print("  fig7: sin datos no-sticking para ley de potencia — omitida")
        return

    Ns   = sorted(sub0["N"].unique())
    Ts   = sorted(sub0["T"].unique())
    epss = sorted(sub0["eps"].unique())

    grp = (sub0.groupby(["T", "eps", "N", "v0"])["COR"]
               .agg(mean="mean", std="std")
               .reset_index())

    fig, axes = plt.subplots(len(Ts), len(epss),
                              figsize=(6 * len(epss), 4 * len(Ts)),
                              squeeze=False)
    fig.suptitle(r"Ley de potencia  COR $\sim A\,v_0^{-\beta}$  (ángulo = 0°, log-log)",
                 fontsize=12)

    for i, T in enumerate(Ts):
        for j, eps in enumerate(epss):
            ax = axes[i][j]
            for N in Ns:
                g = (grp.query("T == @T and eps == @eps and N == @N")
                        .sort_values("v0"))
                if len(g) < 3:
                    continue

                ax.errorbar(g["v0"], g["mean"], yerr=g["std"].fillna(0),
                            fmt=MARKERS[N] + "", color=COLORS_N[N],
                            capsize=3, markersize=5, linestyle="none",
                            label=f"N={N}")

                # log-log linear fit: ln(COR) = ln(A) - β·ln(v0)
                lnv = np.log(g["v0"].to_numpy())
                lnc = np.log(np.maximum(g["mean"].to_numpy(), 1e-6))
                beta, lnA = np.polyfit(lnv, lnc, 1)
                A = np.exp(lnA)
                v_fit = np.linspace(g["v0"].min(), g["v0"].max(), 50)
                ax.plot(v_fit, A * v_fit ** beta, "--", color=COLORS_N[N],
                        linewidth=1.2, alpha=0.7,
                        label=f"  β={abs(beta):.2f}")

            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_title(f"T*={T}, ε_ps={eps}")
            ax.set_xlabel(r"$v_0$  [σ/τ]")
            ax.set_ylabel("COR")
            ax.legend(fontsize=7, ncol=2)

    fig.tight_layout()
    out = OUT_DIR / "fig7_power_law_COR.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out.name}")


# ═════════════════════════════════════════════════════════════════════════════
# 3. Tabla resumen, estadísticas y exportación
# ═════════════════════════════════════════════════════════════════════════════

def print_resumen_tabla(df: pd.DataFrame) -> None:
    """Imprime tabla de medias por (N, T*, v0, eps) para ang=0°."""
    summary = (
        df[df["ang"] == 0.0]
        .groupby(["N", "T", "v0", "eps"], sort=True)
        .agg(
            COR       =("COR",        "mean"),
            alpha_int =("alpha_int",  "mean"),
            alpha_surf=("alpha_surf", "mean"),
            stk_pct   =("sticking",   lambda x: 100.0 * x.mean()),
            n         =("COR",        "count"),
        )
        .reset_index()
    )

    header = (f"{'N':>4}  {'T*':>5}  {'v0':>5}  {'eps':>5}  "
              f"{'COR':>6}  {'α_int':>7}  {'α_surf':>7}  {'stk%':>5}  {'n':>3}")
    sep = "=" * len(header)
    print(f"\n{sep}")
    print("TABLA RESUMEN  (ángulo=0°, promedio sobre réplicas)")
    print(sep)
    print(header)
    print("-" * len(header))
    for _, row in summary.iterrows():
        print(f"{int(row['N']):>4}  {row['T']:>5.1f}  {row['v0']:>5.2f}  {row['eps']:>5}  "
              f"{row['COR']:>6.3f}  {row['alpha_int']:>7.3f}  {row['alpha_surf']:>7.3f}  "
              f"{row['stk_pct']:>5.0f}%  {int(row['n']):>3}")
    print(sep)


def print_describe(df: pd.DataFrame) -> None:
    """Estadísticas descriptivas de las columnas físicas con describe()."""
    cols = ["COR", "alpha_int", "alpha_surf", "E_inc", "We", "balance"]
    print("\n── Estadísticas descriptivas (describe) ─────────────────────────")
    print(df[[c for c in cols if c in df.columns]]
          .describe()
          .round(4)
          .to_string())
    print()


def exportar_csv(df: pd.DataFrame) -> None:
    """Exporta el DataFrame completo y la tabla resumen a CSV."""
    raw_path = OUT_DIR / "resultados_raw.csv"
    df.to_csv(raw_path, index=False)
    print(f"  → {raw_path.name}  ({len(df)} filas, {len(df.columns)} columnas)")

    resumen_path = OUT_DIR / "resultados_resumen.csv"
    (df[df["ang"] == 0.0]
     .groupby(["N", "T", "v0", "eps"])
     .agg(
         COR_mean       =("COR",        "mean"),
         COR_std        =("COR",        "std"),
         alpha_int_mean =("alpha_int",  "mean"),
         alpha_int_std  =("alpha_int",  "std"),
         alpha_surf_mean=("alpha_surf", "mean"),
         alpha_surf_std =("alpha_surf", "std"),
         We_mean        =("We",         "mean"),
         stk_pct        =("sticking",   lambda x: 100.0 * x.mean()),
         n              =("COR",        "count"),
     )
     .reset_index()
     .to_csv(resumen_path, index=False))
    print(f"  → {resumen_path.name}")


# ═════════════════════════════════════════════════════════════════════════════
# Main
# ═════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("Cargando resultados...")
    df = load_all_results()

    if df.empty:
        print("No se encontraron resultados. Verificar DATA_DIR.")
        sys.exit(1)

    print_describe(df)

    print(f"\nGenerando figuras en {OUT_DIR}/")
    fig1_sticking_map(df)
    fig2_balance_energia(df)
    fig3_COR_vs_v0(df)
    fig4_alpha_int_vs_N(df)
    fig5_eps_comparacion(df)
    fig6_T_comparacion(df)
    fig7_power_law_COR(df)

    print_resumen_tabla(df)

    print(f"\nExportando datos a {OUT_DIR}/")
    exportar_csv(df)

    print(f"\nListo. Archivos en {OUT_DIR}/")
    print("  fig1_sticking_map.pdf      — ¿cuándo hay sticking?")
    print("  fig2_balance_energia.pdf   — ¿conserva energía el análisis?")
    print("  fig3_COR_vs_v0.pdf         — COR vs v0 para cada N")
    print("  fig4_alpha_int_vs_N.pdf    — energía interna del clúster vs N")
    print("  fig5_eps_comparacion.pdf   — efecto de ε_ps")
    print("  fig6_T_comparacion.pdf     — efecto de T* de la superficie")
    print("  fig7_power_law_COR.pdf     — ley de potencia COR ~ A·v₀^{-β}")
    print("  resultados_raw.csv         — todos los runs (una fila por run, incluye We)")
    print("  resultados_resumen.csv     — medias por condición (incluye We_mean)")
