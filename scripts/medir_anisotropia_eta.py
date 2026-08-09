"""
medir_anisotropia_eta.py — Mide eta(r,t), la anisotropía de velocidad local
(Sec. 6 de docs/manuscrito/onda_vs_calor_fonones.tex):

    eta(r,t) = (<v_r^2> - <v_perp^2>) / (<v_r^2> + <v_perp^2>)

con v_r la componente de la velocidad de cada átomo a lo largo del radial
3D desde el punto de impacto, y v_perp^2 el promedio POR COMPONENTE de las
dos direcciones transversales restantes (v_perp^2 = (|v|^2 - v_r^2)/2 --
dividir por 2, no dejarlo como suma, es lo que hace que un sistema
isotrópico dé eta=0 exactamente, no eta=-1/3).

eta=1: movimiento puramente radial (coherente, exergía pura -- "sonido").
eta=0: velocidades isotrópicas (equilibrio termodinámico local -- "calor").

Produce dos figuras:
  (1) mapa espacio-tiempo eta(r,t), mismo estilo que el de KE.
  (2) comparación directa de la hipótesis de Gonzalo: eta justo en el
      frente balístico (recién llegado) vs. eta bien atrás, cerca del
      punto de impacto (con más tiempo para dispersarse).

Uso:
    python3 scripts/medir_anisotropia_eta.py \\
        --traj data/outputs/activo/ondas_N79_libre/traj.lammpstrj \\
        --bounce data/outputs/activo/ondas_N79_libre/bounce.dat \\
        --out data/outputs/activo/fig_anisotropia_eta.png \\
        --tmax 20
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import leer_dump as ld

p = argparse.ArgumentParser()
p.add_argument("--traj", required=True)
p.add_argument("--bounce", required=True)
p.add_argument("--out", required=True)
p.add_argument("--tmax", type=float, default=20.0, help="ventana tras impacto [tau]")
p.add_argument("--dt", type=float, default=0.005)
p.add_argument("--nbins", type=int, default=50, help="bins radiales (área igual)")
p.add_argument("--c-sound", type=float, default=8.28, help="velocidad del sonido [sigma/tau]")
args = p.parse_args()

# --- Centro de impacto (x,y) = centroide de la bola; z = tope de la
# superficie sin perturbar (frame 0) -----------------------------------------
cx, cy, info = ld.centro_bola(args.traj)
Lx, Ly, cz = info.Lx, info.Ly, info.z_surf0
print(f"  punto de impacto (x,y,z) = ({cx:.2f}, {cy:.2f}, {cz:.2f}) sigma")

t_imp = ld.t_impacto(args.bounce)
times = info.steps * args.dt
sel = [j for j, tau in enumerate(times) if t_imp <= tau <= t_imp + args.tmax]

r_max = min(Lx, Ly) / 2
r_edges = r_max * np.sqrt(np.linspace(0, 1, args.nbins + 1))   # área igual
r_mid = 0.5 * (r_edges[:-1] + r_edges[1:])


def descomponer(fr):
    """(r_2d, v_r^2, v_perp^2) por átomo de superficie, respecto del impacto."""
    c = fr.col
    surf = fr.data[fr.data[:, c["type"]] == 1]
    x, y, z = (surf[:, c[k]] for k in ("x", "y", "z"))
    vx, vy, vz = (surf[:, c[k]] for k in ("vx", "vy", "vz"))
    dx, dy, dz = x - cx, y - cy, z - cz
    r3 = np.sqrt(dx**2 + dy**2 + dz**2)
    r3 = np.where(r3 < 1e-6, 1e-6, r3)
    v_r = (vx * dx + vy * dy + vz * dz) / r3
    v2_perp = (vx**2 + vy**2 + vz**2 - v_r**2) / 2.0
    # bineo por radio horizontal (consistente con figs previas)
    return np.hypot(dx, dy), v_r**2, v2_perp


# Pasada 2: el mapa eta(r,t) y la comparación frente/estela salen del MISMO
# recorrido -- antes eran dos bucles sobre la lista de frames en RAM, lo que
# obligaba a tenerlos todos cargados y encima calculaba eta_profile dos veces.
print("Pasada 2/2: eta(r,t) y comparación frente/estela ...")
ETA_rt = np.full((args.nbins, len(sel)), np.nan)
t_arr = times[sel]
eta_front, eta_wake = [], []
i2 = 0
for j, fr in enumerate(ld.iter_frames(args.traj, ("type", "x", "y", "z", "vx", "vy", "vz"))):
    if j not in sel:
        continue
    r_2d, v_r2, v2_perp = descomponer(fr)

    prof_r2 = ld.media_por_bin(r_2d, v_r2, r_edges)
    prof_perp2 = ld.media_por_bin(r_2d, v2_perp, r_edges)
    with np.errstate(invalid="ignore", divide="ignore"):
        ETA_rt[:, i2] = (prof_r2 - prof_perp2) / (prof_r2 + prof_perp2)

    r_front = args.c_sound * max(times[j] - t_imp, 0.05)
    m_front = (r_2d >= 0.75 * r_front) & (r_2d < r_front)
    m_wake = r_2d < 0.30 * r_front
    for m, acc in ((m_front, eta_front), (m_wake, eta_wake)):
        if m.sum() > 5:
            a, b = np.mean(v_r2[m]), np.mean(v2_perp[m])
            acc.append((a - b) / (a + b))
        else:
            acc.append(np.nan)
    i2 += 1

eta_front = np.array(eta_front)
eta_wake = np.array(eta_wake)

# --- Figura -------------------------------------------------------------------
fig = plt.figure(figsize=(13, 5))
gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1.0], wspace=0.3)

ax0 = fig.add_subplot(gs[0])
im = ax0.pcolormesh(t_arr - t_imp, r_mid, ETA_rt, cmap="RdBu_r",
                     vmin=-1, vmax=1, shading="auto")
t_line = np.linspace(0.3, args.tmax, 80)
ax0.plot(t_line, args.c_sound * t_line, color="lime", ls=":", lw=1.6,
          label=r"$r=c(t-t_{\rm imp})$")
ax0.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]")
ax0.set_ylabel(r"$r$ (horizontal, desde el impacto)  [$\sigma$]")
ax0.set_title(r"$\eta(r,t)$: rojo=coherente (sonido), azul=isotrópico (calor)")
ax0.legend(fontsize=8, loc="upper left")
cb = fig.colorbar(im, ax=ax0, pad=0.01)
cb.set_label(r"$\eta$")

ax1 = fig.add_subplot(gs[1])
ax1.plot(t_arr - t_imp, eta_front, "o-", ms=4, color="#B71C1C",
          label=r"$\eta$ en el frente (0.75--1.0 $\times\,r_{\rm frente}$)")
ax1.plot(t_arr - t_imp, eta_wake, "s-", ms=4, color="#1565C0",
          label=r"$\eta$ en la estela ($r<0.30\times r_{\rm frente}$, cerca del impacto)")
ax1.axhline(1.0, color="0.5", ls="--", lw=1.0)
ax1.axhline(0.0, color="0.5", ls="--", lw=1.0)
ax1.set_ylim(-0.5, 1.1)
ax1.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]")
ax1.set_ylabel(r"$\eta$")
ax1.set_title("Hipótesis de Gonzalo: ¿frente coherente, estela térmica?")
ax1.legend(fontsize=7.5, loc="lower left")

fig.suptitle(r"Anisotropía de velocidad $\eta(\mathbf{r},t)$ en la superficie "
             f"({os.path.basename(os.path.dirname(args.traj))}, "
             f"ventana {args.tmax:.0f}" r"$\tau$)")
fig.tight_layout()
fig.savefig(args.out, dpi=200, facecolor="white")
print(f"\nFigura guardada: {args.out}")

print("\nResumen eta_frente vs eta_estela (primeros y últimos puntos):")
for i2 in [0, len(sel) // 4, len(sel) // 2, 3 * len(sel) // 4, len(sel) - 1]:
    print(f"  t-t_imp={t_arr[i2]-t_imp:5.2f}  eta_frente={eta_front[i2]:.3f}  "
          f"eta_estela={eta_wake[i2]:.3f}")
