"""
medir_polarizacion_frente.py — Mide la razon de polarizacion

    rho_a = rms(v_z) / rms(v_rho)

del frente de onda en la capa superior de la superficie (v_rho = componente
radial HORIZONTAL desde el impacto, v_z = vertical, v_phi = tangencial).

Motivo (objecion 4 de la bitacora): el criterio eta de la Sec. 6 de
onda_vs_calor_fonones.tex supone "coherente => radial". Eso vale para una onda
longitudinal en un medio infinito, pero NO en una superficie libre, donde el
modo coherente es eliptico en el plano sagital. Con la definicion de
medir_anisotropia_eta.py (v_perp^2 = (|v|^2 - v_r^2)/2) y v_phi = 0:

    eta = (2 - rho_a^2) / (2 + rho_a^2)

de modo que una onda de superficie PERFECTAMENTE coherente con el valor
Rayleigh tipico rho_a ~ 1.5 da eta ~ -0.06, que el criterio lee como "calor".
Este script mide rho_a para refutar (o salvar) esa calibracion CON UN DATO,
no con un argumento.

Diferencia clave con medir_anisotropia_eta.py: el anillo de muestreo NO se
coloca en c*t nominal, sino sobre el frente MEDIDO por energia (mismo umbral
del 15% del pico que usa medir_frente_onda.py). Un anillo colocado en el
c=8.28 nominal puede estar muestreando la zona quieta.

Uso:
    python3 scripts/medir_polarizacion_frente.py \\
        --traj data/outputs/activo/ondas_N79_libre/traj.lammpstrj \\
        --bounce data/outputs/activo/ondas_N79_libre/bounce.dat \\
        --out data/outputs/activo/fig_polarizacion_frente.png \\
        --tmax 7.5
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
p.add_argument("--tmax", type=float, default=7.5, help="ventana tras impacto [tau]")
p.add_argument("--dt", type=float, default=0.005)
p.add_argument("--nbins", type=int, default=60, help="bins radiales (area igual)")
p.add_argument("--umbral", type=float, default=0.15,
               help="fraccion del pico de KE que define el frente")
p.add_argument("--dz-capa", type=float, default=0.5,
               help="espesor [sigma] de la capa superior que se analiza")
args = p.parse_args()

# --- Pasada 1: centro de impacto y nivel de la superficie ---------------------
cx, cy, info = ld.centro_bola(args.traj)
Lx, Ly, z_top = info.Lx, info.Ly, info.z_surf0
print(f"  impacto (x,y) = ({cx:.2f}, {cy:.2f}) sigma; tope de superficie z = {z_top:.3f} sigma")

t_imp = ld.t_impacto(args.bounce)
times = info.steps * args.dt
sel = {j for j, t in enumerate(times) if t_imp <= t <= t_imp + args.tmax}
print(f"  t_imp = {t_imp:.2f} tau; {len(sel)} frames en la ventana de {args.tmax:.1f} tau")

r_max = min(Lx, Ly) / 2
r_edges = r_max * np.sqrt(np.linspace(0, 1, args.nbins + 1))   # bins de area igual
r_mid = 0.5 * (r_edges[:-1] + r_edges[1:])

# --- Pasada 2 ------------------------------------------------------------------
# Todo se calcula sobre la CAPA SUPERIOR: es donde vive el modo de superficie y
# es donde el marco (rho, phi, z) tiene sentido geometrico. Mezclar capas
# profundas (donde el radial 3D se inclina hacia z) es la objecion 1.
print("Pasada 2/2: polarizacion del frente ...")
COLS = ("type", "x", "y", "z", "vx", "vy", "vz", "ke")

t_rel, rho_a, frac_SH, r_frente, eta_medido, eta_pred, n_frente = ([] for _ in range(7))

for j, fr in enumerate(ld.iter_frames(args.traj, COLS)):
    if j not in sel:
        continue
    c = fr.col
    d = fr.data
    capa = (d[:, c["type"]] == 1) & (d[:, c["z"]] > z_top - args.dz_capa)
    d = d[capa]

    dx, dy = d[:, c["x"]] - cx, d[:, c["y"]] - cy
    r2d = np.hypot(dx, dy)
    r2d_safe = np.where(r2d < 1e-6, 1e-6, r2d)
    ux, uy = dx / r2d_safe, dy / r2d_safe          # versor radial horizontal

    vx, vy, vz = (d[:, c[k]] for k in ("vx", "vy", "vz"))
    v_rho = vx * ux + vy * uy                       # radial horizontal
    v_phi = -vx * uy + vy * ux                      # tangencial (SH)
    ke = d[:, c["ke"]]

    # --- El frente se localiza por ENERGIA, no por c*t -------------------------
    perfil_ke = ld.media_por_bin(r2d, ke, r_edges)
    if np.all(np.isnan(perfil_ke)):
        continue
    pico = np.nanmax(perfil_ke)
    sobre = np.where(perfil_ke >= args.umbral * pico)[0]
    if sobre.size == 0:
        continue
    k_frente = sobre[-1]                            # bin mas externo aun "encendido"
    rf = r_mid[k_frente]

    # Anillo sobre el frente medido: los bins encendidos mas externos, con un
    # ancho de ~1 longitud de onda tipica (8 sigma, ver P4) para que el promedio
    # cuadratico sobre la fase sea el de la elipse completa.
    m = (r2d >= max(rf - 8.0, 0.0)) & (r2d <= rf) & (ke >= args.umbral * pico * 0.0)
    if m.sum() < 200:
        continue

    a_rho = np.sqrt(np.mean(v_rho[m] ** 2))
    a_z = np.sqrt(np.mean(vz[m] ** 2))
    a_phi = np.sqrt(np.mean(v_phi[m] ** 2))

    # eta tal como lo define medir_anisotropia_eta.py, pero con el radial
    # horizontal (en la capa superior y a r >> dz el radial 3D coincide):
    #   v_perp^2 = (|v|^2 - v_r^2)/2 = (v_z^2 + v_phi^2)/2
    e_med = (a_rho**2 - 0.5 * (a_z**2 + a_phi**2)) / (a_rho**2 + 0.5 * (a_z**2 + a_phi**2))
    ra = a_z / a_rho

    t_rel.append(times[j] - t_imp)
    r_frente.append(rf)
    rho_a.append(ra)
    frac_SH.append(a_phi**2 / (a_rho**2 + a_z**2 + a_phi**2))
    eta_medido.append(e_med)
    eta_pred.append((2 - ra**2) / (2 + ra**2))      # prediccion onda coherente, v_phi=0
    n_frente.append(int(m.sum()))

t_rel = np.array(t_rel); rho_a = np.array(rho_a); frac_SH = np.array(frac_SH)
r_frente = np.array(r_frente); eta_medido = np.array(eta_medido)
eta_pred = np.array(eta_pred); n_frente = np.array(n_frente)

if t_rel.size == 0:
    raise SystemExit("No hubo frames validos: revisa --tmax / --umbral")

# --- Numero con barra ----------------------------------------------------------
# La ventana [0.3, 1.1] tau es la unica limpia del eco del fondo rigido
# (objecion 5: losa movil de 4.65 sigma, ida y vuelta = 1.12 tau con c=8.28).
# Se excluye t < 0.3 tau porque el frente aun esta bajo la bola.
lim = (t_rel >= 0.3) & (t_rel <= 1.1)
if lim.sum() >= 3:
    ra_m, ra_s = float(np.mean(rho_a[lim])), float(np.std(rho_a[lim]))
else:
    lim = t_rel <= 2.0
    ra_m, ra_s = float(np.mean(rho_a[lim])), float(np.std(rho_a[lim]))

eta_de_ra = (2 - ra_m**2) / (2 + ra_m**2)

# Los per-frame se guardan para no releer 4 GB en cada re-analisis.
npz = os.path.splitext(args.out)[0] + ".npz"
np.savez(npz, t_rel=t_rel, rho_a=rho_a, frac_SH=frac_SH, r_frente=r_frente,
         eta_medido=eta_medido, eta_pred=eta_pred, n_frente=n_frente,
         t_imp=t_imp, cx=cx, cy=cy, z_top=z_top)
print(f"\nDatos por frame: {npz}")

print("\n  t-t_imp   r_frente   rho_a   frac_SH   eta_med   eta_pred   n")
for k in range(len(t_rel)):
    print(f"  {t_rel[k]:6.2f}   {r_frente[k]:7.2f}   {rho_a[k]:5.2f}   "
          f"{frac_SH[k]:6.3f}   {eta_medido[k]:+7.3f}   {eta_pred[k]:+7.3f}   {n_frente[k]:5d}")

print("\n--- RESULTADO -------------------------------------------------")
print(f"  rho_a = rms(v_z)/rms(v_rho) en el frente = {ra_m:.3f} +/- {ra_s:.3f}")
print(f"      (media sobre {lim.sum()} frames, t-t_imp en "
      f"[{t_rel[lim].min():.2f}, {t_rel[lim].max():.2f}] tau)")
print(f"  eta que ESA polarizacion produce si la onda es 100% coherente: {eta_de_ra:+.3f}")
print(f"  fraccion tangencial (SH) media: {np.mean(frac_SH[lim]):.3f}")
print(f"  atomos por anillo: {n_frente[lim].min()}-{n_frente[lim].max()}")
print("\n  Referencia: Rayleigh tipico (nu~0.3) rho_a ~ 1.5 -> eta = -0.06")
print("              longitudinal puro            rho_a = 0   -> eta = +1.00")
print("              isotropo (calor)             rho_a = 1   -> eta =  0.00 (con SH)")
print("---------------------------------------------------------------")

# --- Figura --------------------------------------------------------------------
fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))

ax = axs[0]
ax.plot(t_rel, rho_a, "o-", ms=4, color="#B71C1C")
ax.axhspan(1.3, 1.8, color="#4CAF50", alpha=0.18,
           label=r"Rayleigh tipico ($\nu\approx0.25$–$0.35$)")
ax.axhline(ra_m, color="#B71C1C", ls="--", lw=1.0,
           label=rf"medido $={ra_m:.2f}\pm{ra_s:.2f}$")
ax.axvspan(1.12, t_rel.max(), color="0.6", alpha=0.25)
ax.text(1.18, ax.get_ylim()[1] * 0.95, "eco del fondo\nrigido (1.12$\\tau$)",
        fontsize=7, va="top", color="0.25")
ax.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]")
ax.set_ylabel(r"$\rho_a=\mathrm{rms}(v_z)/\mathrm{rms}(v_\rho)$")
ax.set_title("Polarizacion del frente")
ax.legend(fontsize=7.5, loc="lower right")

ax = axs[1]
ax.plot(t_rel, eta_medido, "o-", ms=4, color="#1565C0", label=r"$\eta$ medido en el anillo")
ax.plot(t_rel, eta_pred, "s--", ms=3.5, color="#F9A825",
        label=r"$\eta=(2-\rho_a^2)/(2+\rho_a^2)$  (onda 100% coherente)")
ax.axhline(0.0, color="0.5", ls=":", lw=1.0)
ax.axvspan(1.12, t_rel.max(), color="0.6", alpha=0.25)
ax.set_ylim(-1.05, 1.05)
ax.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]")
ax.set_ylabel(r"$\eta$")
ax.set_title(r"$\eta$ no distingue: la curva coherente pasa por el dato")
ax.legend(fontsize=7.5, loc="upper left")

ax = axs[2]
ax.plot(t_rel, r_frente, "o-", ms=4, color="#2E7D32", label="frente medido (umbral KE)")
ax.plot(t_rel, 8.28 * t_rel, ":", lw=1.6, color="0.4", label=r"$c_{\rm nom}t$ ($8.28\,\sigma/\tau$)")
ax.axvspan(1.12, t_rel.max(), color="0.6", alpha=0.25)
ax.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]")
ax.set_ylabel(r"$r_{\rm frente}$  [$\sigma$]")
ax.set_title("Donde esta el frente de verdad")
ax.legend(fontsize=7.5, loc="upper left")

fig.suptitle("Objecion 4: el modo coherente de superficie es eliptico, no radial  "
             f"({os.path.basename(os.path.dirname(args.traj))}, capa superior, "
             f"umbral {args.umbral:.0%})")
fig.tight_layout()
fig.savefig(args.out, dpi=200, facecolor="white")
print(f"\nFigura guardada: {args.out}")
