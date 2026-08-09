"""
medir_eco_fondo.py — Mapa <KE>(z,t) por plano atomico bajo el punto de impacto,
para medir CUANDO llega el eco del fondo rigido y fijar la ventana limpia.

Motivo (objecion 5 de la bitacora): la losa movil de in.bounce_cluster tiene
solo ~4.65 sigma de espesor (7 planos FCC entre z=3.87 y z=8.54; las dos capas
de abajo llevan `fix setforce 0 0 0`, o sea pared perfectamente rigida). El
pulso que baja se refleja SIN perdida y vuelve arriba en

    t_eco = 2 * espesor / c = 2*4.65/8.28 = 1.12 tau

que esta DENTRO de todas las ventanas usadas hasta hoy (3.3 tau, 7.5 tau). Todo
lo medido despues de ese instante es superficie + su propia imagen especular.

Este script no supone c: mide el tiempo de llegada por el maximo de <KE> en el
plano de fondo y el retorno en el plano de arriba.

Uso:
    python3 scripts/medir_eco_fondo.py \\
        --traj data/outputs/ondas_N79_libre/traj.lammpstrj \\
        --bounce data/outputs/ondas_N79_libre/bounce.dat \\
        --out data/outputs/fig_eco_fondo.png --tmax 7.5
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
p.add_argument("--tmax", type=float, default=7.5)
p.add_argument("--dt", type=float, default=0.005)
p.add_argument("--rcut", type=float, default=12.0,
               help="radio horizontal [sigma] en torno al impacto; fuera de el "
                    "el pulso ya es lateral y el rebote vertical se diluye")
p.add_argument("--z-fixed-top", type=float, default=3.0992,
               help="tope de la region congelada (2 celdas * a0/2 * 2) [sigma]")
args = p.parse_args()

cx, cy, info = ld.centro_bola(args.traj)
z_top = info.z_surf0
t_imp = ld.t_impacto(args.bounce)
times = info.steps * args.dt
sel = {j for j, t in enumerate(times) if t_imp <= t <= t_imp + args.tmax}
print(f"  t_imp = {t_imp:.2f} tau; capa superior z = {z_top:.3f} sigma; "
      f"fondo congelado z <= {args.z_fixed_top:.3f} sigma")

# Planos FCC: espaciado a0/2 = 0.7748. Se binean con ese paso para que cada bin
# sea UN plano atomico y no una mezcla.
dz = 0.7748
z_edges = np.arange(-dz / 2, z_top + dz, dz)
z_mid = 0.5 * (z_edges[:-1] + z_edges[1:])

print("Pasada 2/2: <KE>(z,t) bajo el impacto ...")
KE_zt, t_rel = [], []
for j, fr in enumerate(ld.iter_frames(args.traj, ("type", "x", "y", "z", "ke"))):
    if j not in sel:
        continue
    c = fr.col
    d = fr.data
    m = d[:, c["type"]] == 1
    d = d[m]
    r2d = np.hypot(d[:, c["x"]] - cx, d[:, c["y"]] - cy)
    cerca = r2d < args.rcut
    KE_zt.append(ld.media_por_bin(d[cerca, c["z"]], d[cerca, c["ke"]], z_edges))
    t_rel.append(times[j] - t_imp)

KE_zt = np.array(KE_zt).T                     # (nz, nt)
t_rel = np.array(t_rel)

# --- Tiempos de llegada medidos ------------------------------------------------
# Plano movil mas bajo y plano de arriba (se ignoran los bins vacios).
vivos = [k for k in range(len(z_mid)) if np.isfinite(KE_zt[k]).any()]
moviles = [k for k in vivos if z_mid[k] > args.z_fixed_top]
k_bot, k_top = moviles[0], moviles[-1]
espesor = z_mid[k_top] - z_mid[k_bot]

perfil_bot = np.nan_to_num(KE_zt[k_bot])
perfil_top = np.nan_to_num(KE_zt[k_top])
t_bot = t_rel[int(np.argmax(perfil_bot))]

# El retorno arriba = primer maximo local del plano superior DESPUES del pico
# de abajo (el primer maximo de todos es el impacto mismo).
post = t_rel > t_bot
t_top2 = t_rel[post][int(np.argmax(perfil_top[post]))] if post.any() else np.nan

print("\n--- RESULTADO -------------------------------------------------")
print(f"  planos moviles: {len(moviles)}  (z = {z_mid[k_bot]:.3f} a {z_mid[k_top]:.3f} sigma)")
print(f"  espesor movil = {espesor:.3f} sigma")
print(f"  pico de <KE> en el plano de FONDO   : t-t_imp = {t_bot:.2f} tau")
print(f"  pico de <KE> en el plano de ARRIBA tras el fondo: t-t_imp = {t_top2:.2f} tau")
if t_bot > 0:
    print(f"  c implicada por la bajada  = {espesor/t_bot:.2f} sigma/tau")
print(f"  ida y vuelta con c=8.28    = {2*espesor/8.28:.2f} tau")
print(f"  => VENTANA LIMPIA del eco vertical: t - t_imp < {2*espesor/8.28:.2f} tau")
print(f"     Con dump cada {times[1]-times[0]:.2f} tau eso son "
      f"{int((2*espesor/8.28)/(times[1]-times[0]))} frames.")
print("---------------------------------------------------------------")

npz = os.path.splitext(args.out)[0] + ".npz"
np.savez(npz, KE_zt=KE_zt, t_rel=t_rel, z_mid=z_mid, t_imp=t_imp,
         espesor=espesor, k_bot=k_bot, k_top=k_top)
print(f"Datos: {npz}")

# --- Figura --------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(12, 4.4), gridspec_kw=dict(width_ratios=[1.25, 1]))

ax = axs[0]
im = ax.pcolormesh(t_rel, z_mid, np.log10(np.maximum(KE_zt, 1e-6)),
                   cmap="inferno", shading="auto")
ax.axhline(args.z_fixed_top, color="cyan", ls="--", lw=1.4, label="tope congelado")
ax.axvline(2 * espesor / 8.28, color="lime", ls=":", lw=1.6,
           label=rf"ida y vuelta $=2d/c={2*espesor/8.28:.2f}\,\tau$")
ax.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]")
ax.set_ylabel(r"$z$  [$\sigma$]")
ax.set_title(rf"$\log_{{10}}\langle KE\rangle(z,t)$, $r<{args.rcut:.0f}\,\sigma$")
ax.legend(fontsize=7.5, loc="upper right")
fig.colorbar(im, ax=ax, pad=0.01)

ax = axs[1]
ax.plot(t_rel, perfil_top, "o-", ms=4, color="#B71C1C",
        label=rf"plano de arriba ($z={z_mid[k_top]:.2f}\,\sigma$)")
ax.plot(t_rel, perfil_bot, "s-", ms=4, color="#1565C0",
        label=rf"plano de fondo movil ($z={z_mid[k_bot]:.2f}\,\sigma$)")
ax.axvline(2 * espesor / 8.28, color="0.3", ls=":", lw=1.6)
ax.axvspan(2 * espesor / 8.28, t_rel.max(), color="0.6", alpha=0.22)
ax.text(2 * espesor / 8.28 + 0.1, ax.get_ylim()[1] * 0.8,
        "contaminado por\nel eco del fondo", fontsize=8, color="0.25")
ax.set_yscale("log")
ax.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]")
ax.set_ylabel(r"$\langle KE\rangle$  [$\epsilon$]")
ax.set_title("El pulso baja, rebota en la pared rigida y vuelve")
ax.legend(fontsize=7.5, loc="lower right")

fig.suptitle("Objecion 5: la losa movil tiene "
             rf"{espesor:.2f}$\,\sigma$ y el fondo es una pared perfectamente reflectante  "
             f"({os.path.basename(os.path.dirname(args.traj))})")
fig.tight_layout()
fig.savefig(args.out, dpi=200, facecolor="white")
print(f"Figura guardada: {args.out}")
