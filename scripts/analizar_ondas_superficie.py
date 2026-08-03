"""
analizar_ondas_superficie.py — Diagnóstico: ¿la energía en la superficie se
propaga como onda (balístico) o como calor (difusivo)? Y hasta qué punto la
losa 8x8x6 con bordes periódicos (x,y) contamina la respuesta por wraparound.

No usa OVITO (no hace falta 3D renderizado): parsea el dump de texto
directamente con numpy, calcula la distancia mínima-imagen (x,y) de cada
átomo de superficie al punto de impacto, y arma:

  (1) un diagrama espacio-tiempo <KE>(r,t) por capas radiales — el "cono de
      luz" balístico, con la velocidad del sonido teórica superpuesta
      (de V''(r_m) en equiparticion_y_virial.tex) y una marca en el tiempo
      en que el frente debería alcanzar el borde de la celda periódica;
  (2) tres instantáneas 2D (x,y) de KE por átomo, para ver directamente el
      patrón de interferencia una vez que la onda se topa con su propia
      imagen periódica.

Uso:
    python3 scripts/analizar_ondas_superficie.py \\
        --traj data/outputs/viz_N79_h28/traj.lammpstrj \\
        --bounce data/outputs/viz_N79_h28/bounce.dat \\
        --out data/outputs/fig_ondas_superficie.png
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
p.add_argument("--dt", type=float, default=0.005)
p.add_argument("--nbins", type=int, default=40)
p.add_argument("--vmax", type=float, default=0.01, help="tope color diagrama r-t")
p.add_argument("--vmax-snap", type=float, default=0.02, help="tope color instantáneas 2D")
p.add_argument("--tmax", type=float, default=20.0, help="ventana tras impacto a graficar [tau]")
args = p.parse_args()

# --- Lectura del dump: streaming en dos pasadas (ver scripts/leer_dump.py) ---
# Pasada 1: centroide de la bola. Pasada 2: KE(r,t) columna a columna y solo
# las 3 instantáneas que se grafican. Nunca hay más de un frame en memoria.
cx, cy, info = ld.centro_bola(args.traj)
steps, Lx, Ly = info.steps, info.Lx, info.Ly

# --- Evento de impacto, desde bounce.dat (para alinear t=0 del diagrama) ----
t_imp = ld.t_impacto(args.bounce)
print(f"  t_imp = {t_imp:.1f} tau")

# velocidad del sonido teórica (orden de magnitud): c = a_nn * sqrt(V''/m)
Vpp, m_atom = 57.15, 1.0
a_nn = 1.5496 / np.sqrt(2)     # vecino más cercano fcc, a0=1.5496 sigma
c_sound = a_nn * np.sqrt(Vpp / m_atom)
t_cross = (Lx / 2) / c_sound
t_wrap = t_imp + t_cross
print(f"  c_sonido (estimado) = {c_sound:.2f} sigma/tau")
print(f"  la onda alcanza el borde periódico en t ~ {t_wrap:.1f} tau "
      f"({t_cross:.2f} tau después del impacto)")

# --- Diagrama espacio-tiempo: <KE>(r,t) por capas radiales, distancia mínima-
# imagen en (x,y) [la caja es periódica en x,y, libre en z] --------------------
r_max = np.hypot(Lx / 2, Ly / 2)
r_edges = np.linspace(0, r_max, args.nbins + 1)
r_mid = 0.5 * (r_edges[:-1] + r_edges[1:])
times = steps * args.dt
KE_rt = np.full((args.nbins, len(steps)), np.nan)

# instantáneas 2D: se eligen los frames ahora para guardarlos al vuelo
snap_times = [t_imp + 0.5 * t_cross, t_imp + 1.0 * t_cross, t_imp + 2.5 * t_cross]
snap_labels = [f"medio cruce ({0.5*t_cross:.2f}τ tras impacto)",
               f"~1 cruce de la celda ({t_cross:.2f}τ)",
               f"~2.5 cruces (wraparound)"]
snap_j = [int(np.argmin(np.abs(times - t_ev))) for t_ev in snap_times]
snaps = {}

print("Pasada 2/2: KE(r,t) e instantáneas ...")
for j, fr in enumerate(ld.iter_frames(args.traj, ("type", "x", "y", "ke"))):
    c = fr.col
    surf = fr.data[fr.data[:, c["type"]] == 1]
    x, y, ke = surf[:, c["x"]], surf[:, c["y"]], surf[:, c["ke"]]
    dx = x - cx; dx -= Lx * np.round(dx / Lx)
    dy = y - cy; dy -= Ly * np.round(dy / Ly)
    KE_rt[:, j] = ld.media_por_bin(np.hypot(dx, dy), ke, r_edges)
    if j in snap_j:
        snaps[j] = (x.copy(), y.copy(), ke.copy())

# --- Figura -------------------------------------------------------------------
fig = plt.figure(figsize=(11, 8.5))
gs = fig.add_gridspec(2, 3, height_ratios=[2.0, 1.0], hspace=0.38, wspace=0.25)

ax0 = fig.add_subplot(gs[0, :])
t_win = (times > t_imp - 5) & (times < t_imp + args.tmax)
im = ax0.pcolormesh(times[t_win], r_mid, KE_rt[:, t_win],
                     cmap="inferno", shading="auto", vmin=0, vmax=args.vmax)
ax0.axvline(t_imp, color="cyan", ls="--", lw=1.1, label=r"$t_{\rm imp}$")
t_line = np.linspace(t_imp, t_imp + r_max / c_sound, 80)
ax0.plot(t_line, c_sound * (t_line - t_imp), color="lime", ls=":", lw=1.8,
          label=rf"$r=c\,(t-t_{{\rm imp}})$, $c\approx{c_sound:.1f}\,\sigma/\tau$")
ax0.axhline(Lx / 2, color="white", ls="--", lw=1.0, alpha=0.8,
             label=r"borde de la celda periódica ($L_x/2$)")
ax0.set_xlabel(r"$t$  [$\tau$]")
ax0.set_ylabel(r"$r$ (distancia mínima-imagen al impacto)  [$\sigma$]")
ax0.set_ylim(0, r_max)
ax0.set_title("Diagrama espacio-tiempo: energía cinética por átomo de superficie, "
              r"$\langle {\rm KE}\rangle(r,t)$ — cono de luz balístico y wraparound")
ax0.legend(loc="upper right", fontsize=8, framealpha=0.9)
cb = fig.colorbar(im, ax=ax0, pad=0.01)
cb.set_label("KE / átomo  [ε]  (promedio por capa radial)")

# --- Instantáneas 2D (x,y) --------------------------------------------------
for col, (j, lab) in enumerate(zip(snap_j, snap_labels)):
    x, y, ke = snaps[j]
    ax = fig.add_subplot(gs[1, col])
    sc = ax.scatter(x, y, c=ke, s=11, cmap="inferno", vmin=0, vmax=args.vmax_snap)
    ax.plot(cx, cy, "+", color="cyan", ms=10, mew=1.5)
    ax.set_aspect("equal")
    ax.set_title(f"{lab}\n$t={times[j]:.2f}\\,\\tau$", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])

fig.suptitle("Superficie 8×8×6, bordes periódicos en (x,y): la onda del impacto "
             "se topa con su propia imagen antes de disiparse", fontsize=10.5, y=0.995)
fig.savefig(args.out, dpi=200, facecolor="white", bbox_inches="tight")
print(f"\nFigura guardada: {args.out}")
