"""
medir_fusion.py — Estado estructural de la bola tras el impacto (punto M5).

Para cada dump (N=79 a distintas velocidades) compara la función de
distribución de pares g(r) de la bola en dos instantes:
  - t = 0            : FCC perfecto (referencia sólida)
  - primer apogeo    : tras el impacto, ¿sigue cristalina o se fundió?

Métricas:
  - g(r): picos FCC nítidos = sólido; fusionados = líquido.
  - T_eff = 2 T_B / (3N-3): temperatura interna efectiva en el apogeo,
    a comparar con la fusión de un cúmulo LJ de ~80 átomos (T* ~ 0.30-0.35).

Genera data/outputs/terminado/fig_fusion_gr.pdf
"""
import numpy as np
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from medir_Nc import iter_frames
from estilo_figuras import aplicar_estilo, panel_label, ANCHO
import matplotlib.pyplot as plt

A0 = 1.5496
R_FCC = [A0/np.sqrt(2), A0, A0*np.sqrt(1.5), A0*np.sqrt(2)]  # picos FCC teóricos
DT = 0.005
RMAX, DR = 2.6, 0.03


def gr_ball(xyz):
    """g(r) no normalizada (histograma de distancias) de un conjunto de átomos."""
    n = len(xyz)
    d = []
    for i in range(n):
        dd = np.sqrt(np.sum((xyz[i+1:] - xyz[i])**2, axis=1))
        d.append(dd)
    d = np.concatenate(d) if d else np.array([])
    bins = np.arange(0, RMAX+DR, DR)
    h, edges = np.histogram(d, bins=bins)
    r = 0.5*(edges[1:]+edges[:-1])
    # normalizar por cascarón esférico y densidad ideal
    shell = 4*np.pi*r**2*DR
    g = h/(shell + 1e-12)
    g = g/np.max(g[r > 0.6])          # normalizar al primer pico
    return r, g


def analizar(dump):
    """Devuelve (r, g0, g1, T_eff_proxy) para t=0 y primer apogeo."""
    frames = []
    for ts, d, ix in iter_frames(dump):
        typ = d[:, ix["type"]].astype(int)
        is_ball = typ == 2
        xyz = d[is_ball][:, [ix["x"], ix["y"], ix["z"]]]
        ke = d[is_ball][:, ix["c_keratom"]]
        frames.append((ts*DT, xyz, d[is_ball][:, ix["z"]].mean(), ke.sum()))
    t = np.array([f[0] for f in frames])
    zcm = np.array([f[2] for f in frames])
    # primer mínimo (impacto) y primer máximo posterior (apogeo)
    i_min = int(np.argmin(zcm[:len(zcm)//2 + 1]))
    i_apo = i_min + int(np.argmax(zcm[i_min:]))
    r, g0 = gr_ball(frames[0][1])
    _, g1 = gr_ball(frames[i_apo][1])
    return r, g0, g1, t[i_apo]


def main():
    aplicar_estilo()
    casos = [
        ("data/outputs/terminado/Nc_runs/N79/traj.lammpstrj",     "$v_{\\rm imp}=1.47$  ($T_{\\rm ef}\\approx0.19$)"),
        ("data/outputs/terminado/Nc_runs/N79_hiv/traj.lammpstrj", "$v_{\\rm imp}=2.03$  ($T_{\\rm ef}\\approx0.28$)"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(ANCHO, 2.7), sharey=True)
    for ax, (dump, lab) in zip(axes, casos):
        if not Path(dump).exists():
            ax.text(0.5, 0.5, "sin datos", transform=ax.transAxes, ha="center")
            continue
        r, g0, g1, tapo = analizar(dump)
        ax.plot(r, g0, "-", color="#1565C0", lw=1.3, label="$t=0$ (FCC)")
        ax.plot(r, g1, "-", color="#C62828", lw=1.3, label="primer apogeo")
        for rf in R_FCC:
            ax.axvline(rf, color="0.7", ls=":", lw=0.6)
        ax.set_xlabel(r"$r$ $(\sigma)$")
        ax.set_title(lab, fontsize=8)
        ax.set_xlim(0.7, RMAX)
    axes[0].set_ylabel(r"$g(r)$ (norm. al 1$^{\rm er}$ pico)")
    axes[0].legend(loc="upper right", fontsize=7)
    panel_label(axes[0], "(a)"); panel_label(axes[1], "(b)")
    out = Path("data/outputs/terminado/fig_fusion_gr.pdf")
    fig.savefig(out); print(f"Figura guardada: {out}")

    # cuantificar: nitidez del 2º pico (r≈a0) vs valle previo
    for dump, lab in casos:
        if not Path(dump).exists():
            continue
        r, g0, g1, tapo = analizar(dump)
        def contraste(g):
            i2 = np.argmin(np.abs(r - R_FCC[1]))           # pico r≈a0
            iv = np.argmin(np.abs(r - 0.5*(R_FCC[0]+R_FCC[1])))  # valle
            return g[i2]/(g[iv]+1e-6)
        print(f"{dump.split('/')[-2]:10s}  contraste pico2/valle: t0={contraste(g0):.2f}  apogeo={contraste(g1):.2f}")


if __name__ == "__main__":
    main()
