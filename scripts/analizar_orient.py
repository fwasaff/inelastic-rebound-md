"""
analizar_orient.py — Análisis del sweep con promedio sobre orientaciones.

Responde a los puntos M1 y M3 del arbitraje:

  M3 (barras de error): cada (N, h0) se corrió para 5 orientaciones iniciales
      del cúmulo. Aquí se reporta media ± σ de cada fracción α.

  M1 (reformular equipartición): en vez de aplicar el conteo de modos a E0
      (inconsistente: ignora el canal CM y sobreestima α_surf), se aplica a la
      energía DISIPADA. La fracción de energía no recuperada que queda en la
      bola,
          f_ball = α_ball / (α_ball + α_surf),
      es la cantidad que la equipartición entre 3(N-1) modos internos y 3·N_c
      osciladores de superficie debería predecir:
          f_ball = (3N-3) / (3N-3 + 3 N_c^eff).
      Se invierte para obtener N_c^eff(N) medido a partir de la partición.

Genera:
  data/outputs/cluster_orient/tabla_orient.txt   (resumen con σ)
  data/outputs/fig_equiparticion_v2.pdf
"""
import numpy as np
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from estilo_figuras import (aplicar_estilo, panel_label, cargar_run,
                            ANCHO, COLOR_N, MARKER_N)

BASE = Path("data/outputs/cluster_orient")
NS = [1, 2, 7, 13, 55, 79, 87]
H0S = [11.03, 18.64, 28.15, 41.47, 58.6, 79.53]
V_PLATEAU = 1.1


def cargar_orient(N, h0):
    """Lista de dicts de partición, uno por orientación válida."""
    out = []
    d = BASE / f"N{N}" / f"h0_{h0}"
    if not d.exists():
        return out
    for od in sorted(d.glob("o*")):
        r = cargar_run(od / "bounce.dat", M=N)
        if r is not None:
            out.append(r)
    return out


def main():
    aplicar_estilo()

    # --- Recolectar todo ------------------------------------------------------
    data = {N: {h0: cargar_orient(N, h0) for h0 in H0S} for N in NS}

    # --- Tabla por (N,h0) con media±σ sobre orientaciones ---------------------
    lines = []
    hdr = (f"{'N':>4} {'h0':>7} {'nO':>3} {'v_imp':>6} "
           f"{'e1':>13} {'a_CM':>13} {'a_TB':>13} {'a_VB':>13} "
           f"{'a_surf':>13} {'a_ball':>13}")
    lines.append(hdr); lines.append("-" * len(hdr))
    for N in NS:
        for h0 in H0S:
            runs = data[N][h0]
            if not runs:
                continue
            def ms(key):
                a = np.array([r[key] for r in runs])
                return a.mean(), a.std()
            v = np.mean([r["v_imp"] for r in runs])
            row = f"{N:>4} {h0:>7} {len(runs):>3} {v:>6.3f}"
            for k in ["e1", "a_CM", "a_TB", "a_VB", "a_surf", "a_ball"]:
                m, s = ms(k)
                row += f" {m:>6.3f}±{s:<5.3f}"
            lines.append(row)
        lines.append("")

    # --- Plateau por N (pool de todas las orientaciones y v>V_PLATEAU) --------
    lines.append("=" * 70)
    lines.append("PLATEAU (v_imp > %.1f), pool de orientaciones y alturas:" % V_PLATEAU)
    lines.append(f"{'N':>4} {'n':>4} {'a_TB':>13} {'a_ball':>13} {'a_surf':>13} "
                 f"{'f_ball':>13} {'Nc_eff':>10}")
    lines.append("-" * 74)
    plateau = {}
    for N in NS:
        runs = [r for h0 in H0S for r in data[N][h0] if r["v_imp"] > V_PLATEAU]
        if not runs:
            continue
        aTB = np.array([r["a_TB"] for r in runs])
        aball = np.array([r["a_ball"] for r in runs])
        asurf = np.array([r["a_surf"] for r in runs])
        # f_ball por run (evita dividir medias); filtrar denominadores ~0
        denom = aball + asurf
        good = denom > 0.02
        fball = aball[good] / denom[good]
        # N_c^eff por inversión de f_ball (solo N>1)
        if N > 1:
            fb_m = fball.mean()
            Nc_eff = (3 * N - 3) * (1 - fb_m) / (3 * fb_m) if fb_m > 0 else np.nan
            # propagación de error simple
            Nc_arr = (3 * N - 3) * (1 - fball) / (3 * fball)
            Nc_s = np.nanstd(Nc_arr)
        else:
            fb_m, Nc_eff, Nc_s = 0.0, np.nan, np.nan
        plateau[N] = dict(
            n=len(runs),
            aTB=(aTB.mean(), aTB.std()),
            aball=(aball.mean(), aball.std()),
            asurf=(asurf.mean(), asurf.std()),
            fball=(fb_m, fball.std() if len(fball) else 0.0),
            Nc=(Nc_eff, Nc_s),
        )
        p = plateau[N]
        lines.append(
            f"{N:>4} {p['n']:>4} "
            f"{p['aTB'][0]:>6.3f}±{p['aTB'][1]:<5.3f} "
            f"{p['aball'][0]:>6.3f}±{p['aball'][1]:<5.3f} "
            f"{p['asurf'][0]:>6.3f}±{p['asurf'][1]:<5.3f} "
            f"{p['fball'][0]:>6.3f}±{p['fball'][1]:<5.3f} "
            f"{p['Nc'][0]:>6.1f}" + ("" if np.isnan(p['Nc'][0]) else f"±{p['Nc'][1]:.0f}"))

    txt = "\n".join(lines)
    out_txt = BASE / "tabla_orient.txt"
    out_txt.write_text(txt)
    print(txt)
    print(f"\nTabla guardada: {out_txt}")

    # --- Figura: f_ball(N) con barras de error y ajuste equipartición ---------
    fig_equiparticion(plateau)
    return plateau


def fig_equiparticion(plateau):
    import matplotlib.pyplot as plt
    Ns = np.array([N for N in NS if N in plateau and N > 1])
    fb = np.array([plateau[N]["fball"][0] for N in Ns])
    fb_s = np.array([plateau[N]["fball"][1] for N in Ns])

    # Ajuste de N_c^eff constante por mínimos cuadrados sobre f_ball(N)
    def model(N, Nc):
        return (3 * N - 3) / (3 * N - 3 + 3 * Nc)
    from scipy.optimize import curve_fit
    m = Ns >= 7                       # excluir dímeros (sin huella real)
    popt, pcov = curve_fit(model, Ns[m], fb[m], p0=[20],
                           sigma=np.maximum(fb_s[m], 1e-3), absolute_sigma=True)
    Nc_fit = popt[0]; Nc_err = np.sqrt(pcov[0, 0])

    fig, ax = plt.subplots(figsize=(ANCHO * 0.62, 2.9))
    Ngrid = np.logspace(np.log10(1.5), np.log10(200), 300)
    ax.plot(Ngrid, model(Ngrid, Nc_fit), "k--", lw=1.3,
            label=rf"$\dfrac{{3N-3}}{{3N-3+3N_c^{{\rm ef}}}}$, $N_c^{{\rm ef}}={Nc_fit:.0f}$")
    ax.axhline(1.0, color="0.6", ls=":", lw=0.8)
    ax.text(150, 0.94, r"$f_{\rm ball}\to 1$", fontsize=7.5, color="0.4", ha="right")
    for N in Ns:
        ax.errorbar(N, plateau[N]["fball"][0], yerr=plateau[N]["fball"][1],
                    marker=MARKER_N[N], color=COLOR_N[N], ms=5, capsize=2, lw=1.0,
                    mfc="none" if N == 13 else None)
    ax.set_xscale("log")
    ax.set_xlabel(r"$N$ (átomos en la bola)")
    ax.set_ylabel(r"$f_{\rm ball}=\alpha_{\rm ball}/(\alpha_{\rm ball}+\alpha_{\rm surf})$")
    ax.set_xlim(1.5, 220)
    ax.set_ylim(0, 1.05)
    ax.set_xticks([2, 7, 13, 55, 79, 87])
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}"))
    ax.legend(loc="lower right", fontsize=7.5)
    out = Path("data/outputs/fig_equiparticion_v2.pdf")
    fig.savefig(out)
    print(f"Figura guardada: {out}")
    print(f"  N_c^eff ajustado (N>=7) = {Nc_fit:.1f} ± {Nc_err:.1f}")


if __name__ == "__main__":
    main()
