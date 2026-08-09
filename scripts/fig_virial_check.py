"""
fig_virial_check.py — Contraste directo del teorema del virial (T_B ~ Delta V_B)
contra los datos del sweep con orientaciones (data/outputs/terminado/cluster_orient).

Para la nota docs/manuscrito/equiparticion_y_virial.tex:
  - Panel (a): alpha_TB y alpha_VB vs v_imp para N=79 (cluster canonico),
    media +/- sigma sobre 5 orientaciones. El virial armonico predice
    que ambas curvas coincidan.
  - Panel (b): razon alpha_VB/alpha_TB vs v_imp para N=13,55,79,87,
    con la prediccion virial (razon=1) como referencia. La razon por-run
    evita dividir por el promedio (mismo criterio que fig_equiparticion_v2).

Genera: data/outputs/terminado/fig_virial_check.pdf
"""
import numpy as np
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from estilo_figuras import (aplicar_estilo, panel_label, cargar_run,
                            ANCHO, COLOR_N, MARKER_N)

BASE = Path("data/outputs/terminado/cluster_orient")
H0S = [11.03, 18.64, 28.15, 41.47, 58.6, 79.53]
A_TB_MIN = 0.02   # descarta razones VB/TB cuando TB es consistente con ruido


def cargar_orient(N, h0):
    out = []
    d = BASE / f"N{N}" / f"h0_{h0}"
    if not d.exists():
        return out
    for od in sorted(d.glob("o*")):
        r = cargar_run(od / "bounce.dat", M=N)
        if r is not None:
            out.append(r)
    return out


def stats_N(N):
    """v_imp, a_TB (mean,std), a_VB (mean,std), razon (mean,std) por h0."""
    rows = []
    for h0 in H0S:
        runs = cargar_orient(N, h0)
        if not runs:
            continue
        v = np.mean([r["v_imp"] for r in runs])
        aTB = np.array([r["a_TB"] for r in runs])
        aVB = np.array([r["a_VB"] for r in runs])
        good = aTB > A_TB_MIN
        razon = aVB[good] / aTB[good] if good.any() else np.array([])
        rows.append(dict(
            v=v,
            aTB=(aTB.mean(), aTB.std()),
            aVB=(aVB.mean(), aVB.std()),
            razon=(razon.mean(), razon.std()) if len(razon) else (np.nan, np.nan),
            n_razon=int(good.sum()),
        ))
    return rows


def main():
    aplicar_estilo()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(ANCHO, 2.9))

    # --- Panel (a): alpha_TB, alpha_VB vs v para N=79 -------------------------
    N79 = stats_N(79)
    v = np.array([r["v"] for r in N79])
    aTB_m = np.array([r["aTB"][0] for r in N79])
    aTB_s = np.array([r["aTB"][1] for r in N79])
    aVB_m = np.array([r["aVB"][0] for r in N79])
    aVB_s = np.array([r["aVB"][1] for r in N79])

    ax1.errorbar(v, aTB_m, yerr=aTB_s, marker="o", color="#d62728", ms=4.5,
                 capsize=2, lw=1.2, label=r"$\alpha_{T_B}$ (cinética)")
    ax1.errorbar(v, aVB_m, yerr=aVB_s, marker="s", color="#1f77b4", ms=4.5,
                 capsize=2, lw=1.2, label=r"$\alpha_{\Delta V_B}$ (potencial)")
    ax1.set_xlabel(r"$v_{\rm imp}\ (\sigma/\tau)$")
    ax1.set_ylabel(r"fracción de $E_0$")
    ax1.legend(loc="upper left", fontsize=7.5)
    ax1.set_ylim(bottom=0)
    panel_label(ax1, "(a)")

    # --- Panel (b): razon alpha_VB/alpha_TB vs v para varios N ----------------
    for N in [13, 55, 79, 87]:
        rows = stats_N(N)
        vv = np.array([r["v"] for r in rows])
        rm = np.array([r["razon"][0] for r in rows])
        rs = np.array([r["razon"][1] for r in rows])
        m = ~np.isnan(rm)
        ax2.errorbar(vv[m], rm[m], yerr=rs[m], marker=MARKER_N[N],
                     color=COLOR_N[N], ms=5, capsize=2, lw=1.0,
                     label=f"$N={N}$")
    ax2.axhline(1.0, color="0.3", ls="--", lw=1.0)
    ax2.text(1.85, 1.05, "virial\n" r"($\alpha_{\Delta V_B}=\alpha_{T_B}$)",
              fontsize=7, color="0.3", ha="right")
    ax2.set_xlabel(r"$v_{\rm imp}\ (\sigma/\tau)$")
    ax2.set_ylabel(r"$\alpha_{\Delta V_B}/\alpha_{T_B}$")
    ax2.set_ylim(0, 2.4)
    ax2.legend(loc="upper left", fontsize=7, ncol=2)
    panel_label(ax2, "(b)")

    fig.tight_layout()
    out = Path("data/outputs/terminado/fig_virial_check.pdf")
    fig.savefig(out)
    print(f"Figura guardada: {out}")

    # --- Resumen numerico para la nota -----------------------------------------
    print("\nN=79, v_imp / a_TB / a_VB / razon:")
    for r in N79:
        print(f"  v={r['v']:.3f}  aTB={r['aTB'][0]:.3f}±{r['aTB'][1]:.3f}  "
              f"aVB={r['aVB'][0]:.3f}±{r['aVB'][1]:.3f}  "
              f"razon={r['razon'][0]:.2f}±{r['razon'][1]:.2f}")


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    main()
