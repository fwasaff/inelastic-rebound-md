"""
medir_chi_coherencia.py — Separa sonido de calor con la descomposicion estandar
de la MD de choque: se resta la velocidad de flujo local ANTES de llamar "calor"
a lo que queda.

    v_i = u(celda) + dv_i         u = (1/n) sum_i v_i  (masas iguales)

    E_coh  = 1/2 n |u|^2          movimiento colectivo  -> sonido
    E_term = 1/2 sum_i |dv_i|^2   fluctuacion en torno al flujo -> calor

    chi = E_coh / (E_coh + E_term)        1 = coherente, 0 = termalizado

Por que esto y no eta: eta = (<v_r^2>-<v_perp^2>)/(<v_r^2>+<v_perp^2>) mide
POLARIZACION, no coherencia. Medido en este mismo repositorio
(scripts/medir_polarizacion_frente.py): eta y la prediccion de una onda 100%
coherente con la razon de elipse rho_a medida correlacionan a 0.9989, con un
residuo del 6% de la variacion total. eta no contiene informacion mas alla de
rho_a, y una onda de Rayleigh limpia (rho_a~1.5) da eta~-0.06, que el criterio
viejo leia como "calor". chi es agnostico a la polarizacion: una Rayleigh
coherente da chi~1.

TRES PERILLAS, y las tres se reportan en vez de elegirse:

1. Tamano de celda l. Demasiado grande borra el movimiento colectivo (chi
   subestimado); demasiado pequeno mete ruido termico en u (chi sobreestimado).
   Se corren tres. Restriccion fisica: l << lambda ~ 9.7 sigma (medido en
   medir_pulso_contacto.py).

2. Piso de muestreo. Una celda de n atomos en equilibrio TERMICO PURO no da
   chi=0 sino chi = 1/n, porque u fluctua: Var(u_x)=sigma_v^2/n, luego
   E[|u|^2]=3 sigma_v^2/n y E_coh = (3/2) sigma_v^2 mientras E_tot = (3/2) n
   sigma_v^2. El piso se dibuja en la figura; sin el, cualquier chi pequeno
   parece "algo de coherencia".

3. Umbral de energia. Las celdas por debajo del umbral (donde el frente no ha
   llegado) tienen chi indefinido y ruidoso. Se reporta para varios umbrales.

Nota sobre la fase: para una onda coherente E_coh y E_tot oscilan JUNTAS (la KE
total pasa por cero cada medio periodo, pero tambien lo hace la parte
colectiva), asi que el COCIENTE es robusto a la fase. Lo unico que hace falta es
no dividir 0/0, y de eso se encarga el umbral de energia.

Celdas CARTESIANAS, no anillos: promediar el vector (vx,vy,vz) sobre un anillo
completo cancela vx y vy por simetria y daria u_horizontal=0 para una onda
radial perfectamente coherente. El bineo por r se hace DESPUES, sobre las
energias ya calculadas celda a celda.

Uso:
    python3 scripts/medir_chi_coherencia.py \\
        --traj data/outputs/activo/ondas_N79_libre/traj.lammpstrj \\
        --bounce data/outputs/activo/ondas_N79_libre/bounce.dat \\
        --out data/outputs/activo/fig_chi_coherencia.png --tmax 7.5
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
p.add_argument("--celdas", type=float, nargs="+", default=[3.10, 4.65, 6.20],
               help="lados de celda [sigma] (2a0, 3a0, 4a0)")
p.add_argument("--umbrales", type=float, nargs="+", default=[0.05, 0.15, 0.30],
               help="fraccion del pico de KE por celda bajo la cual se descarta")
p.add_argument("--nbins", type=int, default=40)
p.add_argument("--z-fixed-top", type=float, default=3.0992,
               help="las capas congeladas (v=0 siempre) hay que EXCLUIRLAS: "
                    "meterlas baja u artificialmente y falsea chi")
p.add_argument("--t-eco", type=float, default=1.8,
               help="llegada medida del eco del fondo [tau], para sombrear")
args = p.parse_args()

cx, cy, info = ld.centro_bola(args.traj)
Lx, Ly = info.Lx, info.Ly
t_imp = ld.t_impacto(args.bounce)
times = info.steps * args.dt
sel = {j for j, t in enumerate(times) if t_imp <= t <= t_imp + args.tmax}
print(f"  t_imp = {t_imp:.2f} tau; {len(sel)} frames; celdas {args.celdas} sigma")

r_max = min(Lx, Ly) / 2
r_edges = r_max * np.sqrt(np.linspace(0, 1, args.nbins + 1))
r_mid = 0.5 * (r_edges[:-1] + r_edges[1:])
NC, NU = len(args.celdas), len(args.umbrales)


def chi_por_celda(d, c, lado):
    """(r_celda, E_coh, E_tot, n) por celda cubica de lado `lado`."""
    x, y, z = (d[:, c[k]] for k in ("x", "y", "z"))
    v = np.stack([d[:, c[k]] for k in ("vx", "vy", "vz")], axis=1)

    ix = np.floor((x - x.min()) / lado).astype(np.int64)
    iy = np.floor((y - y.min()) / lado).astype(np.int64)
    # En z la losa movil mide solo ~5.4 sigma. Con una rejilla de paso `lado`
    # anclada en z.min() la ultima fila queda como una lonja parcial (con
    # lado=4.65 son 0.79 sigma), y esas celdas casi vacias distorsionan el piso
    # 1/n. Se reparte el espesor en un numero ENTERO de capas iguales.
    espesor = z.max() - z.min()
    nz = max(1, int(round(espesor / lado)))
    dz = espesor / nz
    iz = np.clip(((z - z.min()) / dz).astype(np.int64), 0, nz - 1)
    nx, ny = ix.max() + 1, iy.max() + 1
    cid = (iz * ny + iy) * nx + ix
    cid, inv = np.unique(cid, return_inverse=True)
    ncel = cid.size

    n = np.bincount(inv, minlength=ncel).astype(np.float64)
    su = np.stack([np.bincount(inv, weights=v[:, k], minlength=ncel) for k in range(3)], axis=1)
    s2 = np.bincount(inv, weights=(v ** 2).sum(axis=1), minlength=ncel)

    E_tot = 0.5 * s2                                   # KE total de la celda
    E_coh = 0.5 * (su ** 2).sum(axis=1) / n            # = 1/2 n |u|^2
    # r de la celda = centroide de sus atomos (mas fiel que el centro geometrico
    # en las celdas del borde, que estan parcialmente llenas)
    cxm = np.bincount(inv, weights=x, minlength=ncel) / n
    cym = np.bincount(inv, weights=y, minlength=ncel) / n
    return np.hypot(cxm - cx, cym - cy), E_coh, E_tot, n


print("Pasada 2/2: chi(r,t) ...")
CHI = np.full((NC, NU, args.nbins, len(sel)), np.nan)
CHI_frente = np.full((NC, NU, len(sel)), np.nan)
CHI_estela = np.full((NC, NU, len(sel)), np.nan)
PISO = np.zeros((NC, NU))
t_rel = np.full(len(sel), np.nan)
i2 = 0

for j, fr in enumerate(ld.iter_frames(args.traj, ("type", "x", "y", "z", "vx", "vy", "vz"))):
    if j not in sel:
        continue
    c = fr.col
    d = fr.data
    movil = (d[:, c["type"]] == 1) & (d[:, c["z"]] > args.z_fixed_top)
    d = d[movil]
    t_rel[i2] = times[j] - t_imp

    for a, lado in enumerate(args.celdas):
        rc, E_coh, E_tot, n = chi_por_celda(d, c, lado)
        ke_cel = E_tot / n                              # KE media por atomo
        pico = ke_cel.max()
        for b, umb in enumerate(args.umbrales):
            m = ke_cel >= umb * pico
            if m.sum() < 4:
                continue
            PISO[a, b] += np.mean(1.0 / n[m])           # se promedia y se
                                                        # normaliza al final
            # chi por bin radial: se agregan ENERGIAS, no se promedian chi
            # (promediar cocientes le da el mismo peso a una celda muerta que
            # al frente entero)
            num = ld.media_por_bin(rc[m], E_coh[m], r_edges)
            den = ld.media_por_bin(rc[m], E_tot[m], r_edges)
            with np.errstate(invalid="ignore", divide="ignore"):
                CHI[a, b, :, i2] = num / den

            # Frente = celdas encendidas mas externas (1 lambda_R de ancho);
            # estela = tercio interior del radio alcanzado.
            rf = rc[m].max()
            f = m & (rc >= rf - 9.7)
            w = m & (rc < 0.33 * rf)
            for mask, out in ((f, CHI_frente), (w, CHI_estela)):
                if mask.sum() >= 4:
                    out[a, b, i2] = E_coh[mask].sum() / E_tot[mask].sum()
    i2 += 1

PISO /= len(sel)

npz = os.path.splitext(args.out)[0] + ".npz"
np.savez(npz, CHI=CHI, CHI_frente=CHI_frente, CHI_estela=CHI_estela, PISO=PISO,
         t_rel=t_rel, r_mid=r_mid, celdas=np.array(args.celdas),
         umbrales=np.array(args.umbrales), t_imp=t_imp)
print(f"Datos: {npz}")

# --- Resultado -----------------------------------------------------------------
lim = t_rel <= args.t_eco          # ventana limpia del eco del fondo
print(f"\n--- RESULTADO (ventana limpia t-t_imp <= {args.t_eco} tau, "
      f"{lim.sum()} frames) ---")
print(f"{'celda':>7} {'n_arm':>7} {'piso':>7} {'umbral':>7} "
      f"{'chi_frente':>11} {'chi_corr':>9} {'chi_estela':>11}")
# chi_corr descuenta la fuga termica al canal colectivo. Si una fraccion
# chi_true de la energia es coherente, la parte termica restante aporta
# ademas E_term/n al numerador, luego
#     chi_med = chi_true + (1-chi_true)/n   =>   chi_true = (chi_med - 1/n)/(1 - 1/n)
CHI_CORR = np.zeros((NC, NU))
for a, lado in enumerate(args.celdas):
    for b, umb in enumerate(args.umbrales):
        cf = np.nanmean(CHI_frente[a, b, lim])
        ce = np.nanmean(CHI_estela[a, b, lim])
        pi = PISO[a, b]
        CHI_CORR[a, b] = (cf - pi) / (1 - pi)
        print(f"{lado:7.2f} {1/pi:7.1f} {pi:7.3f} {umb:7.2f} "
              f"{cf:11.3f} {CHI_CORR[a,b]:9.3f} {ce:11.3f}")
print("-" * 70)
if NC >= 3:
    b0 = min(1, NU - 1)
    d = np.diff(CHI_CORR[:, b0])
    print(f"  chi_corr vs tamano de celda: "
          f"{'MONOTONA DECRECIENTE (sin plateau)' if np.all(d < 0) else 'no monotona'}")
    print(f"  variacion de chi_corr en el rango de celdas: "
          f"{CHI_CORR[:,b0].min():.3f} a {CHI_CORR[:,b0].max():.3f} "
          f"(factor {CHI_CORR[:,b0].max()/max(CHI_CORR[:,b0].min(),1e-9):.1f})")
    print("  Sin plateau, chi NO es un numero: es una funcion de la perilla.")
print("-" * 70)

# --- Figura --------------------------------------------------------------------
fig, axs = plt.subplots(2, 2, figsize=(12.5, 8.4))
b_ref = min(1, NU - 1)                                  # umbral central
COL = [plt.cm.viridis(v) for v in np.linspace(0.05, 0.85, max(NC, NU))]

ax = axs[0, 0]
im = ax.pcolormesh(t_rel, r_mid, CHI[0, b_ref], cmap="magma", vmin=0, vmax=1,
                   shading="auto")
ax.plot(t_rel, 4.6 * t_rel, ":", lw=1.6, color="cyan", label=r"$c_R t$ (medido)")
ax.axvline(args.t_eco, color="lime", ls="--", lw=1.5, label="eco del fondo")
ax.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]")
ax.set_ylabel(r"$r$  [$\sigma$]")
ax.set_title(rf"$\chi(r,t)$  (celda ${args.celdas[0]:.2f}\,\sigma$, "
             rf"umbral {args.umbrales[b_ref]:.0%})")
ax.legend(fontsize=7.5, loc="upper left")
fig.colorbar(im, ax=ax, pad=0.01).set_label(r"$\chi$")

ax = axs[0, 1]
for a, lado in enumerate(args.celdas):
    ax.plot(t_rel, CHI_frente[a, b_ref], "o-", ms=3.5, color=COL[a],
            label=rf"$\ell={lado:.2f}\,\sigma$ ($n\approx{1/PISO[a,b_ref]:.0f}$)")
    ax.axhline(PISO[a, b_ref], color=COL[a], ls=":", lw=1.0)
ax.axvspan(args.t_eco, t_rel.max(), color="0.6", alpha=0.25)
ax.text(args.t_eco + 0.1, 0.05, "eco del fondo", fontsize=8, color="0.25")
ax.set_ylim(0, 1.05)
ax.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]")
ax.set_ylabel(r"$\chi$")
ax.set_title(r"Frente. Punteado = piso de muestreo $1/n$")
ax.legend(fontsize=7.5, loc="lower left")

ax = axs[1, 0]
for a, lado in enumerate(args.celdas):
    ax.plot(t_rel, CHI_estela[a, b_ref], "s-", ms=3.5, color=COL[a],
            label=rf"$\ell={lado:.2f}\,\sigma$")
    ax.axhline(PISO[a, b_ref], color=COL[a], ls=":", lw=1.0)
ax.axvspan(args.t_eco, t_rel.max(), color="0.6", alpha=0.25)
ax.set_ylim(0, 1.05)
ax.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]")
ax.set_ylabel(r"$\chi$")
ax.set_title(r"Estela ($r<0.33\,r_{\rm max}$). H1 pide que baje al piso")
ax.legend(fontsize=7.5, loc="lower left")

ax = axs[1, 1]
lados = np.array(args.celdas)
cf = np.array([np.nanmean(CHI_frente[a, b_ref, lim]) for a in range(NC)])
ax.plot(lados, cf, "o-", ms=5, color="#B71C1C", label=r"$\chi$ medido (frente)")
ax.plot(lados, CHI_CORR[:, b_ref], "s--", ms=5, color="#1565C0",
        label=r"$\chi$ corregido por el piso")
ax.plot(lados, PISO[:, b_ref], "^:", ms=4, color="0.45", label=r"piso $1/n$")
ax.axvspan(0, 1.1, color="#EF9A9A", alpha=0.35)
ax.text(0.15, 0.92, "$n\\lesssim4$\n(piso domina)", fontsize=7, color="#B71C1C")
ax.axvspan(5.0, lados.max() * 1.05, color="#90CAF9", alpha=0.35)
ax.text(5.2, 0.92, "$\\ell\\gtrsim a_{\\rm fuente}$\n(promedia la fase)",
        fontsize=7, color="#1565C0")
ax.set_xlim(0, lados.max() * 1.05)
ax.set_ylim(0, 1.05)
ax.set_xlabel(r"lado de celda $\ell$  [$\sigma$]")
ax.set_ylabel(r"$\chi$")
ax.set_title(r"Convergencia: sin plateau, $\chi$ no es un numero")
ax.legend(fontsize=7.5, loc="center right")

fig.suptitle(r"$\chi=E_{\rm coh}/(E_{\rm coh}+E_{\rm term})$: sonido vs calor "
             "sin suponer polarizacion  "
             f"({os.path.basename(os.path.dirname(args.traj))})")
fig.tight_layout()
fig.savefig(args.out, dpi=200, facecolor="white")
print(f"Figura guardada: {args.out}")
