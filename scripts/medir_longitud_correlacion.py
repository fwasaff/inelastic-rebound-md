"""
medir_longitud_correlacion.py — Correlacion espacial del campo de velocidades en
la capa superior de la superficie:

    C(d) = < v_i . v_j >_{|r_ij| = d}  /  < |v|^2 >

Por que este observable y no chi: chi = E_coh/(E_coh+E_term) exige un tamano de
celda de grano grueso l con a_red << l << lambda. Medido en este repositorio
(medir_chi_coherencia.py), chi cae monotonamente de 0.73 a 0.17 al variar l
entre 1.55 y 6.20 sigma, SIN plateau: no existe ese l, asi que chi no es un
numero sino una funcion de la perilla.

C(d) no tiene perillas. Y distingue las dos hipotesis por la FORMA, no por un
valor:

  - Onda coherente de longitud lambda: v(x) ~ A cos(kx - wt), luego
    C(d) ~ cos(k d): OSCILA con periodo lambda y no decae en pocas distancias
    interatomicas.
  - Campo termico: C(d) cae a ~0 en una o dos distancias interatomicas y no
    vuelve a subir.

PERO la normalizacion importa, y mucho. Se calculan TRES estimadores sobre los
mismos pares porque el resultado cambia de raiz segun cual se use (esto ES la
objecion 8 de la bitacora, cerrada el 2026-08-09):

  1. direccion, sum(v_i^.v_j^)/N -- acotada pero CIEGA A LA AMPLITUD.
  2. amplitud,  sum(v_i.v_j)/sum(|v_i||v_j|) -- LA BUENA: pesa por amplitud y
     sigue acotada, porque Cauchy-Schwarz vale termino a termino.
  3. ingenua,   sum(v_i.v_j)/(N<|v|^2>) -- la dominan unos pocos atomos
     rapidos del crater; llego a dar |C| ~ 5, imposible.

Con el estimador 1 la coherencia lejana (8-15 sigma) parecia decaer de 0.50 a 0
con tau ~ 0.25 tau. Con el estimador 2 ese decaimiento NO EXISTE: es ruido en
torno a cero desde el primer frame. El decaimiento eran los atomos QUIETOS del
campo lejano, cuyas velocidades minusculas estan debilmente alineadas por el
campo elastico cuasiestatico y que el estimador 1 pesa igual que los del frente.

Resultado con el estimador bueno: xi = 2.12 sigma (~2 vecinos), C ~ 0 mas alla
de 4 sigma, y plano en el tiempo. La energia que entra en la superficie no
tiene estructura coherente mas alla de un par de distancias interatomicas --
dentro de una ventana que, por lo demas, NO es valida para la pregunta (el
forzamiento sigue activo y el eco del fondo llega a 1.8 tau).

Uso:
    python3 scripts/medir_longitud_correlacion.py \\
        --traj data/outputs/activo/ondas_N79_libre/traj.lammpstrj \\
        --bounce data/outputs/activo/ondas_N79_libre/bounce.dat \\
        --out data/outputs/activo/fig_correlacion_velocidad.png --tmax 3.0
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
p.add_argument("--tmax", type=float, default=3.0)
p.add_argument("--dt", type=float, default=0.005)
p.add_argument("--rmax-region", type=float, default=30.0,
               help="radio [sigma] de la region analizada en torno al impacto")
p.add_argument("--rmin-region", type=float, default=5.0,
               help="se excluye el crater: unos pocos atomos bajo la bola llevan "
                    "|v| ~100x la tipica y, al correlacionar velocidades sin "
                    "normalizar, un solo par de ellos domina su bin y produce "
                    "|C|>1 (imposible por Cauchy-Schwarz). Ver bitacora 2026-08-09")
p.add_argument("--dmax", type=float, default=25.0, help="separacion maxima [sigma]")
p.add_argument("--nd", type=int, default=100, help="bins de separacion")
p.add_argument("--dz-capa", type=float, default=0.5)
p.add_argument("--t-eco", type=float, default=1.8)
p.add_argument("--max-atomos", type=int, default=4000,
               help="submuestreo: N^2 pares crece rapido y 4000 ya da 8e6 pares")
args = p.parse_args()

cx, cy, info = ld.centro_bola(args.traj)
z_top = info.z_surf0
t_imp = ld.t_impacto(args.bounce)
times = info.steps * args.dt
sel = {j for j, t in enumerate(times) if t_imp <= t <= t_imp + args.tmax}
print(f"  t_imp = {t_imp:.2f} tau; {len(sel)} frames")

d_edges = np.linspace(0, args.dmax, args.nd + 1)
d_mid = 0.5 * (d_edges[:-1] + d_edges[1:])
rng = np.random.default_rng(0)

print("Pasada 2/2: C(d) en la capa superior ...")
CD, CD_amp, CD_cru, t_rel = [], [], [], []

for j, fr in enumerate(ld.iter_frames(args.traj, ("type", "x", "y", "z", "vx", "vy", "vz"))):
    if j not in sel:
        continue
    c = fr.col
    d = fr.data
    m = ((d[:, c["type"]] == 1) & (d[:, c["z"]] > z_top - args.dz_capa))
    d = d[m]
    x, y = d[:, c["x"]], d[:, c["y"]]
    rr0 = np.hypot(x - cx, y - cy)
    dentro = (rr0 < args.rmax_region) & (rr0 >= args.rmin_region)
    d = d[dentro]
    if d.shape[0] < 100:
        for L in (CD, CD_amp, CD_cru):
            L.append(np.full(args.nd, np.nan))
        t_rel.append(times[j] - t_imp)
        continue
    if d.shape[0] > args.max_atomos:
        d = d[rng.choice(d.shape[0], args.max_atomos, replace=False)]

    pos = np.stack([d[:, c["x"]], d[:, c["y"]]], axis=1).astype(np.float64)
    v = np.stack([d[:, c[k]] for k in ("vx", "vy", "vz")], axis=1).astype(np.float64)

    # Marco CILINDRICO local (rho, phi, z) centrado en el impacto.
    # Correlacionar los vectores crudos mide la geometria, no la coherencia: en
    # una onda radial dos atomos en lados opuestos del anillo se mueven en
    # sentidos opuestos, asi que v_i.v_j < 0 por pura rotacion del versor, aunque
    # la onda sea perfectamente coherente. En la base cilindrica ese giro se
    # descuenta y lo que queda es la correlacion de fase.
    dxr, dyr = pos[:, 0] - cx, pos[:, 1] - cy
    rr = np.maximum(np.hypot(dxr, dyr), 1e-9)
    ur, us = dxr / rr, dyr / rr
    v = np.stack([v[:, 0] * ur + v[:, 1] * us,       # v_rho
                  -v[:, 0] * us + v[:, 1] * ur,      # v_phi
                  v[:, 2]], axis=1)                  # v_z

    nv = np.maximum(np.linalg.norm(v, axis=1, keepdims=True), 1e-12)
    vhat = v / nv                                    # versores
    vmod = nv[:, 0]                                  # |v| por atomo

    # Pares i<j: la matriz completa de 4000x4000 son 16e6 float64 = 128 MB,
    # asi que se recorre por bloques de filas.
    #
    # Se acumulan TRES estimadores sobre los mismos pares, porque la objecion 8
    # es precisamente cual de ellos vale:
    #   num_dir  = sum(v_i^ . v_j^)          -> C_dir = num_dir/cnt
    #        acotada, pero CIEGA A LA AMPLITUD: un atomo quieto pesa igual que
    #        uno del frente. Es el mismo defecto que la objecion 3 le achacaba
    #        a eta.
    #   num_amp / den_amp = sum(v_i.v_j) / sum(|v_i||v_j|)
    #        LA BUENA: pesa por amplitud (los atomos energeticos mandan) y sigue
    #        acotada en [-1,1], porque Cauchy-Schwarz da |v_i.v_j| <= |v_i||v_j|
    #        TERMINO A TERMINO, asi que el cociente de las sumas tambien lo esta.
    #   num_cru  = sum(v_i.v_j) / (cnt * <|v|^2>)
    #        la ingenua: es la que daba |C| ~ 5 dominada por el crater.
    num_dir = np.zeros(args.nd)
    num_amp = np.zeros(args.nd)
    den_amp = np.zeros(args.nd)
    num_cru = np.zeros(args.nd)
    cnt = np.zeros(args.nd)
    B = 500
    for i0 in range(0, pos.shape[0], B):
        i1 = min(i0 + B, pos.shape[0])
        dd = np.linalg.norm(pos[i0:i1, None, :] - pos[None, :, :], axis=2)
        vv_dir = vhat[i0:i1] @ vhat.T
        vv_raw = v[i0:i1] @ v.T
        mm = vmod[i0:i1, None] * vmod[None, :]
        # solo j > i, para no contar el par dos veces ni el autopar (d=0)
        jj = np.arange(pos.shape[0])[None, :] > np.arange(i0, i1)[:, None]
        k = np.digitize(dd[jj], d_edges) - 1
        ok = (k >= 0) & (k < args.nd)
        kk = k[ok]
        num_dir += np.bincount(kk, weights=vv_dir[jj][ok], minlength=args.nd)
        num_amp += np.bincount(kk, weights=vv_raw[jj][ok], minlength=args.nd)
        den_amp += np.bincount(kk, weights=mm[jj][ok], minlength=args.nd)
        num_cru += np.bincount(kk, weights=vv_raw[jj][ok], minlength=args.nd)
        cnt += np.bincount(kk, minlength=args.nd)

    # Se exigen >=200 pares por bin para que un bin flaco no invente estructura.
    v2 = float(np.mean(vmod ** 2))
    with np.errstate(invalid="ignore", divide="ignore"):
        vale = cnt >= 200
        CD.append(np.where(vale, num_dir / np.maximum(cnt, 1), np.nan))
        CD_amp.append(np.where(vale, num_amp / np.maximum(den_amp, 1e-30), np.nan))
        CD_cru.append(np.where(vale, num_cru / np.maximum(cnt * v2, 1e-30), np.nan))
    t_rel.append(times[j] - t_imp)

CD = np.array(CD)
CD_amp = np.array(CD_amp)
CD_cru = np.array(CD_cru)
t_rel = np.array(t_rel)

# --- Numeros -------------------------------------------------------------------
lim = t_rel <= args.t_eco
C_lim = np.nanmean(CD[lim], axis=0)

# Longitud de correlacion: primer cruce por 1/e
bajo = np.where(C_lim < 1 / np.e)[0]
xi = d_mid[bajo[0]] if bajo.size else np.nan

# Deteccion de onda: SOLO cuenta si C(d) baja de -0.05 en al menos 3 bins
# consecutivos. Un plateau positivo con ruido produce cruces por cero
# espurios; con el criterio laxo (primer cambio de signo) este mismo dataset
# daba "lambda = 16 sigma" y luego "24 sigma" segun el preprocesado, o sea que
# el numero lo fijaba el ruido, no la fisica.
neg = np.nan_to_num(C_lim) < -0.05
racha = np.convolve(neg.astype(int), np.ones(3, int), "same")
hay_onda = bool((racha >= 3).any())
if hay_onda:
    i_min = int(np.nanargmin(C_lim))
    lam = 2 * d_mid[i_min]
else:
    i_min, lam = -1, np.nan

# Meseta lejana: valor de C mas alla de varias xi. Si el campo estuviera
# termalizado, C -> 0 ahi. Un plateau positivo delata movimiento colectivo de
# toda la region (deriva comun), que es coherencia de escala grande.
lejos = (d_mid >= 8.0) & (d_mid <= 15.0)
C_far = np.nanmean(C_lim[lejos])
C_far_t = np.array([np.nanmean(CD[k][lejos]) for k in range(len(t_rel))])
C_far_amp = np.array([np.nanmean(CD_amp[k][lejos]) for k in range(len(t_rel))])
C_far_cru = np.array([np.nanmean(CD_cru[k][lejos]) for k in range(len(t_rel))])
print(f"\n  Rango de cada estimador (control de acotamiento, |C|<=1):")
for nom, A in (("direccion", CD), ("amplitud", CD_amp), ("ingenuo", CD_cru)):
    print(f"     {nom:10s}: min {np.nanmin(A):+8.3f}  max {np.nanmax(A):+8.3f}"
          f"{'   <-- SE SALE DE [-1,1]' if np.nanmax(np.abs(A)) > 1.01 else ''}")

a_red = 1.5496 / np.sqrt(2)          # vecino mas cercano FCC = a0/sqrt(2)

print("\n--- RESULTADO (ventana limpia t-t_imp <= "
      f"{args.t_eco} tau, {lim.sum()} frames) ---")
print(f"  distancia al primer vecino (FCC)      a = {a_red:.2f} sigma")
print(f"  longitud de correlacion (C = 1/e)    xi = {xi:.2f} sigma  = {xi/a_red:.1f} a")
print(f"  meseta lejana <C>_(8-15 sigma)          = {C_far:+.3f}")
print()
if hay_onda:
    print(f"  ONDA RESUELTA: minimo en {d_mid[i_min]:.2f} sigma -> lambda ~ {lam:.1f} sigma")
else:
    print("  NO hay onda periodica resuelta: C(d) no baja de -0.05 en ningun")
    print("  tramo. No hay longitud de onda medible en este observable, y por")
    print("  tanto tampoco un lambda con el que justificar el tamano de celda")
    print("  de grano grueso de chi.")
print()
print("  Coherencia de gran escala <C>(8-15 sigma) frente al tiempo:")
print(f"     {'t-t_imp':>8} {'direccion':>10} {'AMPLITUD':>10} {'ingenuo':>9}")
for k in range(len(t_rel)):
    marca = "" if t_rel[k] <= args.t_eco else "  <- tras el eco"
    print(f"     {t_rel[k]:8.2f} {C_far_t[k]:+10.3f} {C_far_amp[k]:+10.3f} "
          f"{C_far_cru[k]:+9.3f}{marca}")

# Tiempo de decaimiento del estimador BUENO (amplitud), por ajuste exponencial
# sobre la ventana limpia y mientras la senal siga por encima del ruido.
mfit = (t_rel <= args.t_eco) & (C_far_amp > 0.02)
if mfit.sum() >= 3:
    A = np.polyfit(t_rel[mfit], np.log(C_far_amp[mfit]), 1)
    print(f"\n  AJUSTE (estimador de amplitud, {mfit.sum()} puntos, t<={args.t_eco}):")
    print(f"     tau_decaimiento = {-1/A[0]:.3f} tau")
else:
    print("\n  AJUSTE: no hay suficientes puntos sobre el ruido en la ventana limpia")
print("-" * 66)

npz = os.path.splitext(args.out)[0] + ".npz"
np.savez(npz, CD=CD, t_rel=t_rel, d_mid=d_mid, C_lim=C_lim, xi=xi,
         CD_amp=CD_amp, CD_cru=CD_cru, C_far=C_far, C_far_t=C_far_t,
         C_far_amp=C_far_amp, C_far_cru=C_far_cru, hay_onda=hay_onda,
         a_red=a_red, t_imp=t_imp)
print(f"Datos: {npz}")

# --- Figura --------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))

ax = axs[0]
for k in range(len(t_rel)):
    col = plt.cm.viridis(k / max(len(t_rel) - 1, 1))
    ax.plot(d_mid, CD[k], "-", lw=1.0, color=col, alpha=0.75)
ax.plot(d_mid, C_lim, "-", lw=2.4, color="#B71C1C",
        label=rf"media, $t-t_{{\rm imp}}\leq{args.t_eco}\,\tau$")
ax.axhline(0, color="0.5", lw=0.8)
ax.axhline(1 / np.e, color="0.5", ls=":", lw=0.9)
ax.axvline(a_red, color="#1565C0", ls="--", lw=1.2,
           label=rf"1er vecino $a={a_red:.2f}\,\sigma$")
ax.axvline(xi, color="#B71C1C", ls=":", lw=1.4, label=rf"$\xi={xi:.2f}\,\sigma$")
ax.set_xlim(0, args.dmax)
ax.set_xlabel(r"separacion $d$  [$\sigma$]")
ax.set_ylabel(r"$C(d)=\langle \mathbf{v}_i\cdot\mathbf{v}_j\rangle/\langle|v|^2\rangle$")
ax.set_title("Correlacion de velocidad (color = tiempo)")
ax.legend(fontsize=7.5)

ax = axs[1]
ax.plot(t_rel, C_far_amp, "o-", ms=5, color="#B71C1C",
        label="pesada por amplitud (la buena)")
ax.plot(t_rel, C_far_t, "s--", ms=4, color="#1565C0",
        label="solo direccion (ciega a la amplitud)")
ax.legend(fontsize=7.5)
ax.axhline(0, color="0.5", lw=0.8)
ax.axvspan(args.t_eco, t_rel.max(), color="0.6", alpha=0.25)
ax.text(args.t_eco + 0.05, 0.9 * np.nanmax(C_far_t), "eco del\nfondo",
        fontsize=8, color="0.25")
ax.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]")
ax.set_ylabel(r"$\langle C\rangle_{8\text{--}15\,\sigma}$")
ax.set_title("Coherencia de gran escala frente al tiempo\n"
             "(sin celda, sin piso, sin umbral)")

fig.suptitle("Sin perillas: la forma de $C(d)$ decide onda vs calor  "
             f"({os.path.basename(os.path.dirname(args.traj))}, capa superior, "
             rf"$r<{args.rmax_region:.0f}\,\sigma$)")
fig.tight_layout()
fig.savefig(args.out, dpi=200, facecolor="white")
print(f"Figura guardada: {args.out}")
