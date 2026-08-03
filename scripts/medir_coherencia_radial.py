"""
medir_coherencia_radial.py — ¿Onda cerca del frente, calor cerca del impacto?

Idea de Gonzalo: en un instante fijo, la energía que acaba de llegar al
borde del disturbio todavía es un pulso elástico coherente (ecuación de
onda); la energía cerca del punto de impacto tuvo más tiempo (~duración del
contacto, ~10 tau) para dispersarse por anarmonicidad y termalizar (ecuación
de calor). No es un exponente que cambia con el tiempo -- es una estructura
espacial en un solo snapshot: anillos discretos (coherente) afuera, perfil
liso y monótono (termalizado) adentro.

Mide, para cada radio, cuánta "textura" oscilatoria tiene el perfil KE(r)
en una ventana deslizante: resta una versión suavizada (pasa-bajos) del
perfil crudo y toma la energía del residuo -- alta textura = todavía
anillos/coherente: baja textura = liso/termalizado.

Uso:
    python3 scripts/medir_coherencia_radial.py \\
        --traj data/outputs/ondas_N79_libre/traj.lammpstrj \\
        --bounce data/outputs/ondas_N79_libre/bounce.dat \\
        --out data/outputs/fig_coherencia_radial.png \\
        --times 3 5 7
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import leer_dump as ld

p = argparse.ArgumentParser()
p.add_argument("--traj", required=True)
p.add_argument("--bounce", required=True)
p.add_argument("--out", required=True)
p.add_argument("--times", type=float, nargs="+", default=[3.0, 5.0, 7.0],
               help="instantes tras el impacto a examinar [tau]")
p.add_argument("--dt", type=float, default=0.005)
p.add_argument("--nbins", type=int, default=120, help="bins radiales finos")
p.add_argument("--smooth-win", type=int, default=9,
               help="ancho (en bins) del suavizado pasa-bajos")
args = p.parse_args()

cx, cy, info = ld.centro_bola(args.traj)
Lx, Ly = info.Lx, info.Ly

t_imp = ld.t_impacto(args.bounce)
times = info.steps * args.dt

r_max = min(Lx, Ly) / 2
# bins de AREA igual (no ancho igual): r_k = r_max*sqrt(k/nbins). Con bins de
# ancho igual, el área (~ #átomos) de cada anillo crece como r, así que cerca
# de r=0 hay muchos menos átomos por bin y el ruido de muestreo es mayor --
# eso fabricaría "textura" espuria cerca del centro sin ninguna física real.
# Con área igual, el ruido de muestreo es comparable en todo r.
r_edges = r_max * np.sqrt(np.linspace(0, 1, args.nbins + 1))
r_mid = 0.5 * (r_edges[:-1] + r_edges[1:])


def smooth(y, w):
    kernel = np.ones(w) / w
    pad = w // 2
    yp = np.pad(y, pad, mode="edge")
    return np.convolve(yp, kernel, mode="valid")[:len(y)]


# Pasada 2: de cada instante pedido solo sobrevive su perfil radial (nbins
# floats), no las posiciones de los 1.5e5 átomos.
print("Pasada 2/2: perfiles radiales en los instantes pedidos ...")
j_pedidos = [int(np.argmin(np.abs(times - (t_imp + t_ev)))) for t_ev in args.times]
perfiles = {}
for j, fr in enumerate(ld.iter_frames(args.traj, ("type", "x", "y", "ke"))):
    if j not in j_pedidos:
        continue
    c = fr.col
    surf = fr.data[fr.data[:, c["type"]] == 1]
    r = np.hypot(surf[:, c["x"]] - cx, surf[:, c["y"]] - cy)
    perfiles[j] = ld.media_por_bin(r, surf[:, c["ke"]], r_edges)

fig, axes = plt.subplots(2, len(args.times), figsize=(5.2 * len(args.times), 7.0),
                          sharex="col")

for col, (t_ev, j) in enumerate(zip(args.times, j_pedidos)):
    prof = perfiles[j]
    valid = np.isfinite(prof)
    prof_i = np.interp(r_mid, r_mid[valid], prof[valid])

    trend = smooth(prof_i, args.smooth_win)
    residual = prof_i - trend
    texture = smooth(residual ** 2, args.smooth_win) ** 0.5   # RMS local del residuo

    r_front_now = 8.28 * t_ev  # c_sonido de equiparticion_y_virial.tex

    ax_top, ax_bot = axes[0, col], axes[1, col]
    ax_top.plot(r_mid, prof_i, color="#1565C0", lw=1.0)
    ax_top.plot(r_mid, trend, color="#B71C1C", lw=1.3, ls="--", label="tendencia (pasa-bajos)")
    ax_top.axvline(r_front_now, color="0.5", ls=":", lw=1.2,
                    label=r"$r=c(t{-}t_{\rm imp})$")
    ax_top.set_title(rf"$t-t_{{\rm imp}}={t_ev:.0f}\,\tau$")
    ax_top.set_ylabel("KE / átomo  [ε]" if col == 0 else "")
    ax_top.legend(fontsize=7)

    ax_bot.fill_between(r_mid, texture, color="#2E7D32", alpha=0.5)
    ax_bot.axvline(r_front_now, color="0.5", ls=":", lw=1.2)
    ax_bot.set_xlabel(r"$r$  [$\sigma$]")
    ax_bot.set_ylabel("textura local (RMS residuo)" if col == 0 else "")

axes[0, 0].annotate("cerca del impacto\n(¿liso = termalizado?)", xy=(0.05, 0.9),
                     xycoords="axes fraction", fontsize=8, color="0.3", ha="left")
fig.suptitle("¿Onda en el frente, calor cerca del impacto? Perfil radial de KE y su "
             "'textura' oscilatoria local (línea punteada gris = posición balística)")
fig.tight_layout()
fig.savefig(args.out, dpi=200, facecolor="white")
print(f"\nFigura guardada: {args.out}")
