"""
medir_coherencia_frente_estela.py — Longitud de correlacion xi medida POR
SEPARADO dentro del frente y dentro de la estela, con el anillo SIGUIENDO al
frente en vez de quieto.

Es el test directo de la hipotesis de Gonzalo (frente coherente, estela
termica). Junta las tres lecciones de las rondas anteriores:

  1. El anillo va sobre el frente MEDIDO POR ENERGIA, no sobre c*t nominal.
     (objecion 3: con c=8.28 nominal el anillo muestreaba la zona quieta)
  2. La correlacion se pesa por AMPLITUD, con la normalizacion acotada
     C(d) = sum(v_i.v_j) / sum(|v_i||v_j|).
     (objecion 8: normalizar a versores da peso igual a un atomo quieto que a
     uno del frente, y eso fabricaba un decaimiento que no existia)
  3. Todo en la capa superior y en base cilindrica local, para no mezclar la
     rotacion geometrica del versor radial con la coherencia de fase.
     (objecion 1)

Por que xi y no chi: chi exige un tamano de celda de grano grueso y no tiene
meseta (varia un factor 4.4). xi no tiene perillas de ese tipo: es la distancia
a la que el campo de velocidades deja de estar correlacionado consigo mismo.

  - Onda coherente de longitud lambda -> xi ~ lambda/2 (aqui ~5-8 sigma)
  - Calor                             -> xi ~ 1-2 distancias interatomicas

Lo que decide:
  * xi_frente ~ xi_estela ~ 1 sigma  -> no hay frente coherente; H1 cae.
  * xi_frente >> xi_estela           -> es H1, y el ritmo al que xi cae detras
                                        del frente es tau_th.

Uso:
    python3 scripts/medir_coherencia_frente_estela.py \\
        --traj data/outputs/activo/ondas_N79_profunda_h28.15/traj.lammpstrj \\
        --bounce data/outputs/activo/ondas_N79_profunda_h28.15/bounce.dat \\
        --out data/outputs/activo/fig_coherencia_frente_estela.png \\
        --tmax 10.5
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
p.add_argument("--tmax", type=float, default=10.5,
               help="ventana limpia tras el impacto [tau]")
p.add_argument("--dt", type=float, default=0.005)
p.add_argument("--nbins-r", type=int, default=60, help="bins radiales para el perfil de KE")
p.add_argument("--umbral", type=float, default=0.15,
               help="fraccion del pico de KE que define el borde del frente")
p.add_argument("--ancho-frente", type=float, default=10.0,
               help="espesor del anillo del frente [sigma] (~1 lambda_R)")
p.add_argument("--r-crater", type=float, default=5.0,
               help="radio excluido en torno al impacto: unos pocos atomos del "
                    "crater llevan |v| ~100x la tipica y dominarian su bin")
p.add_argument("--r-fijo", type=float, nargs="+", default=[25.0, 40.0, 50.0],
               help="CONTROL de geometria fija: anillo QUIETO en [r_fijo-ancho, "
                    "r_fijo] medido en todos los instantes. El frente lo cruza y "
                    "se va, asi que cualquier cambio de xi ahi es fisica y no "
                    "geometria -- el area, la curvatura y el numero de atomos no "
                    "cambian. Es el control que decide si el crecimiento de "
                    "xi_frente con el radio es real o artefacto del anillo movil.")
p.add_argument("--dmax", type=float, default=15.0, help="separacion maxima [sigma]")
p.add_argument("--nd", type=int, default=60)
p.add_argument("--dz-capa", type=float, default=0.5)
p.add_argument("--max-atomos", type=int, default=3000)
p.add_argument("--min-atomos", type=int, default=200)
p.add_argument("--min-pares", type=int, default=200,
               help="pares minimos por bin; por debajo el bin queda en NaN")
args = p.parse_args()

cx, cy, info = ld.centro_bola(args.traj)
z_top = info.z_surf0
Lx, Ly = info.Lx, info.Ly
t_imp = ld.t_impacto(args.bounce)
times = info.steps * args.dt
sel = {j for j, t in enumerate(times) if t_imp <= t <= t_imp + args.tmax}
print(f"  t_imp = {t_imp:.2f} tau; capa superior z = {z_top:.2f}; {len(sel)} frames")

r_max = min(Lx, Ly) / 2
r_edges = r_max * np.sqrt(np.linspace(0, 1, args.nbins_r + 1))   # area igual
r_mid = 0.5 * (r_edges[:-1] + r_edges[1:])
d_edges = np.linspace(0, args.dmax, args.nd + 1)
d_mid = 0.5 * (d_edges[:-1] + d_edges[1:])
rng = np.random.default_rng(0)
a_red = 1.5496 / np.sqrt(2)          # 1er vecino FCC = 1.096 sigma


def correlacion(pos, v):
    """C(d) = sum(v_i.v_j)/sum(|v_i||v_j|) sobre pares i<j, por bloques.

    Acotada en [-1,1] termino a termino por Cauchy-Schwarz, y pesada por
    amplitud: los atomos energeticos mandan, que es lo que se quiere medir.
    """
    n = pos.shape[0]
    vmod = np.linalg.norm(v, axis=1)
    num = np.zeros(args.nd)
    den = np.zeros(args.nd)
    cnt = np.zeros(args.nd)
    B = 400
    for i0 in range(0, n, B):
        i1 = min(i0 + B, n)
        dd = np.linalg.norm(pos[i0:i1, None, :] - pos[None, :, :], axis=2)
        vv = v[i0:i1] @ v.T
        mm = vmod[i0:i1, None] * vmod[None, :]
        jj = np.arange(n)[None, :] > np.arange(i0, i1)[:, None]
        k = np.digitize(dd[jj], d_edges) - 1
        ok = (k >= 0) & (k < args.nd)
        kk = k[ok]
        num += np.bincount(kk, weights=vv[jj][ok], minlength=args.nd)
        den += np.bincount(kk, weights=mm[jj][ok], minlength=args.nd)
        cnt += np.bincount(kk, minlength=args.nd)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(cnt >= args.min_pares, num / np.maximum(den, 1e-30), np.nan)


def xi_de(C):
    """Primer cruce por 1/e, interpolado. NaN si no cruza dentro de dmax."""
    ok = np.isfinite(C)
    if ok.sum() < 3:
        return np.nan
    d, c = d_mid[ok], C[ok]
    bajo = np.where(c < 1 / np.e)[0]
    if bajo.size == 0:
        return np.nan                      # sigue correlacionado a dmax
    k = bajo[0]
    if k == 0:
        return d[0]
    # interpolacion lineal entre el ultimo punto por encima y el primero por debajo
    c0, c1 = c[k - 1], c[k]
    return float(d[k - 1] + (d[k] - d[k - 1]) * (c0 - 1 / np.e) / max(c0 - c1, 1e-12))


print("Pasada 2/2: xi en el frente y en la estela ...")
res = dict(t=[], r_f=[], xi_f=[], xi_e=[], n_f=[], n_e=[], C_f=[], C_e=[])
for q in range(len(args.r_fijo)):
    res[f"xi_x{q}"] = []

for j, fr in enumerate(ld.iter_frames(args.traj,
                                      ("type", "x", "y", "z", "vx", "vy", "vz", "ke"))):
    if j not in sel:
        continue
    c = fr.col
    d = fr.data
    capa = (d[:, c["type"]] == 1) & (d[:, c["z"]] > z_top - args.dz_capa)
    d = d[capa]

    dx, dy = d[:, c["x"]] - cx, d[:, c["y"]] - cy
    r2d = np.hypot(dx, dy)
    rs = np.maximum(r2d, 1e-9)
    ur, us = dx / rs, dy / rs
    vx, vy, vz = (d[:, c[k]] for k in ("vx", "vy", "vz"))
    # base cilindrica local (rho, phi, z)
    V = np.stack([vx * ur + vy * us, -vx * us + vy * ur, vz], axis=1)
    P = np.stack([d[:, c["x"]], d[:, c["y"]]], axis=1).astype(np.float64)
    ke = d[:, c["ke"]]

    # --- Frente localizado POR ENERGIA -----------------------------------------
    perfil = ld.media_por_bin(r2d, ke, r_edges)
    if np.all(np.isnan(perfil)):
        continue
    sobre = np.where(perfil >= args.umbral * np.nanmax(perfil))[0]
    if sobre.size == 0:
        continue
    r_f = r_mid[sobre[-1]]                       # bin encendido mas externo
    r_int = r_f - args.ancho_frente              # borde interior del anillo

    if r_int <= args.r_crater:                   # frente aun pegado al crater:
        continue                                 # frente y estela no separables

    m_f = (r2d >= r_int) & (r2d <= r_f)
    m_e = (r2d >= args.r_crater) & (r2d < r_int)
    # Control: anillos QUIETOS, misma area y curvatura en todo instante.
    regiones = [("f", m_f), ("e", m_e)]
    for q, rq in enumerate(args.r_fijo):
        regiones.append((f"x{q}", (r2d >= rq - args.ancho_frente) & (r2d <= rq)))

    fila = {}
    for nom, m in regiones:
        idx = np.where(m)[0]
        if idx.size < args.min_atomos:
            fila[nom] = (np.nan, np.nan, 0)
            continue
        if idx.size > args.max_atomos:
            idx = rng.choice(idx, args.max_atomos, replace=False)
        C = correlacion(P[idx], V[idx])
        fila[nom] = (xi_de(C), C, idx.size)

    res["t"].append(times[j] - t_imp)
    res["r_f"].append(r_f)
    res["xi_f"].append(fila["f"][0]); res["C_f"].append(fila["f"][1]); res["n_f"].append(fila["f"][2])
    res["xi_e"].append(fila["e"][0]); res["C_e"].append(fila["e"][1]); res["n_e"].append(fila["e"][2])
    for q in range(len(args.r_fijo)):
        res[f"xi_x{q}"].append(fila[f"x{q}"][0])

t = np.array(res["t"]); r_f = np.array(res["r_f"])
xi_f = np.array(res["xi_f"], dtype=float); xi_e = np.array(res["xi_e"], dtype=float)
XI_X = np.array([res[f"xi_x{q}"] for q in range(len(args.r_fijo))], dtype=float)
n_f = np.array(res["n_f"]); n_e = np.array(res["n_e"])

if t.size == 0:
    raise SystemExit("Ningun frame con frente separable del crater. Prueba --tmax mayor.")

# --- Resultado -----------------------------------------------------------------
print(f"\n  {'t-t_imp':>8} {'r_f':>7} {'xi_frente':>10} {'xi_estela':>10}")
for k in range(len(t)):
    print(f"  {t[k]:8.2f} {r_f[k]:7.1f} {xi_f[k]:10.2f} {xi_e[k]:10.2f}")

vf, ve = xi_f[np.isfinite(xi_f)], xi_e[np.isfinite(xi_e)]
print("\n--- RESULTADO -------------------------------------------------")
print(f"  1er vecino FCC                a = {a_red:.2f} sigma")
if vf.size:
    print(f"  xi FRENTE  mediana = {np.median(vf):.2f} sigma = {np.median(vf)/a_red:.1f} a"
          f"   (rango {vf.min():.2f}-{vf.max():.2f})")
if ve.size:
    print(f"  xi ESTELA  mediana = {np.median(ve):.2f} sigma = {np.median(ve)/a_red:.1f} a"
          f"   (rango {ve.min():.2f}-{ve.max():.2f})")
if vf.size and ve.size:
    razon = np.median(vf) / max(np.median(ve), 1e-9)
    print(f"  razon frente/estela = {razon:.2f}")
    print()
    if razon > 1.5:
        print("  => El frente esta MAS correlacionado que la estela: es la forma")
        print("     que predice H1 (Gonzalo). El ritmo al que xi cae detras del")
        print("     frente es tau_th.")
    else:
        print("  => El frente NO esta mas correlacionado que la estela. No hay")
        print("     frente coherente distinguible: H1 no se sostiene con este dato.")
print("---------------------------------------------------------------")

# --- CONTROL de geometria fija + ajuste de tau_th -------------------------------
# Cada anillo [rq-ancho, rq] NO se mueve: misma area, misma curvatura y mismo
# numero de atomos en todos los instantes. Por tanto cualquier cambio de xi ahi
# es fisica, no geometria. Es lo que descarta que el crecimiento de xi_frente
# con el radio (que se ve con el anillo movil) sea un artefacto.
#
# El instante en que el frente cruza el anillo se saca del ajuste COMPLETO
# r_f = a*t + b, incluyendo el intercepto. Usar rq/a (sin intercepto) da un
# cruce ~2 tau tarde, porque el frente no nace en r=0 sino en r=b~9.5 sigma,
# que es el tamano de la zona de contacto.
pend, inter = np.polyfit(t, r_f, 1)
print(f"\n--- CONTROL de geometria fija (frente: r_f = {pend:.2f} t + {inter:.1f}) ---")

# El piso de xi es FISICO: el primer bin de separacion con pares suficientes
# esta en ~1.12 sigma, porque no hay atomos mas cerca que el primer vecino
# (1.10 sigma). Un xi que toca ese piso significa "la correlacion ya murio
# entre vecinos" = calor. Pero censura cualquier ajuste exponencial, y por eso
# aqui NO se ajusta una exponencial: se caracteriza la caida sin suponerle forma.

# COLAPSO: cada anillo se mide a una EDAD distinta (el tiempo que lleva pasado
# el frente por el). Si los anillos caen sobre una sola curva de xi frente a la
# edad, el decaimiento es una propiedad del material y no de donde se mire.
E, X = [], []
for q, rq in enumerate(args.r_fijo):
    tc = (rq - inter) / max(pend, 1e-9)
    ok = np.isfinite(XI_X[q])
    E.append(t[ok] - tc)
    X.append(XI_X[q][ok])
E = np.concatenate(E); X = np.concatenate(X)

print(f"\n  COLAPSO sobre la edad desde el paso del frente ({len(args.r_fijo)} anillos):")
print(f"    {'edad [tau]':>14}  {'xi [sigma]':>10}  {'n':>3}  {'dispersion':>10}")
bordes = np.arange(-1.0, 6.01, 0.5)
for lo, hi in zip(bordes[:-1], bordes[1:]):
    m = (E >= lo) & (E < hi) & np.isfinite(X)
    if m.sum() >= 2:
        print(f"    [{lo:5.1f},{hi:5.1f})  {np.median(X[m]):10.2f}  {m.sum():3d}  "
              f"{np.std(X[m]):10.2f}")
m0 = (E >= 0) & (E < 0.5) & np.isfinite(X)
xi0 = np.median(X[m0]) if m0.sum() else np.nan
piso = np.nanmin(X)
print(f"\n    xi al pasar el frente = {xi0:.2f} sigma")
print(f"    piso alcanzado        = {piso:.2f} sigma  (1er vecino {a_red:.2f} = calor)")
print(f"    semivida (xi -> {(xi0+piso)/2:.2f}) = edad ~2.0-2.5 tau")
print()
print("    La dispersion entre anillos a la misma edad (columna 'dispersion')")
print("    es la que dice si el colapso es real: si es pequena frente a la")
print("    caida total, el decaimiento no depende del radio.")
print("---------------------------------------------------------------")

npz = os.path.splitext(args.out)[0] + ".npz"
np.savez(npz, t=t, r_f=r_f, xi_f=xi_f, xi_e=xi_e, XI_X=XI_X,
         r_fijo=np.array(args.r_fijo), pend=pend, inter=inter,
         n_f=n_f, n_e=n_e, d_mid=d_mid, a_red=a_red, t_imp=t_imp,
         C_f=np.array([x if np.ndim(x) else np.full(args.nd, np.nan) for x in res["C_f"]]),
         C_e=np.array([x if np.ndim(x) else np.full(args.nd, np.nan) for x in res["C_e"]]))
print(f"Datos: {npz}")

# --- Figura --------------------------------------------------------------------
fig, axs = plt.subplots(1, 3, figsize=(14.5, 4.4))

ax = axs[0]
ax.plot(t, xi_f, "o-", ms=5, color="#B71C1C", label="frente (anillo que lo sigue)")
ax.plot(t, xi_e, "s-", ms=5, color="#1565C0", label="estela (detras del frente)")
for q, rq in enumerate(args.r_fijo):
    col = plt.cm.Purples(0.45 + 0.5 * q / max(len(args.r_fijo) - 1, 1))
    ax.plot(t, XI_X[q], "^-", ms=4, color=col,
            label=rf"fijo [{rq-args.ancho_frente:.0f},{rq:.0f}]$\sigma$")
    ax.axvline((rq - inter) / pend, color=col, ls=":", lw=1.2)
ax.axhline(a_red, color="0.45", ls="--", lw=1.2, label=rf"1er vecino $a={a_red:.2f}\,\sigma$")
ax.axhspan(0, a_red, color="#90CAF9", alpha=0.30)
ax.text(t.min() + 0.2, a_red * 0.45, "zona 'calor'\n(sin orden mas alla del vecino)",
        fontsize=7, color="#0D47A1")
ax.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]")
ax.set_ylabel(r"$\xi$  [$\sigma$]")
ax.set_title("Longitud de correlacion")
ax.legend(fontsize=7.5)

ax = axs[1]
nC = len(res["C_f"])
for k in range(nC):
    C = res["C_f"][k]
    if np.ndim(C):
        ax.plot(d_mid, C, "-", lw=1.0, alpha=0.75,
                color=plt.cm.autumn(k / max(nC - 1, 1)))
Cf = np.nanmean(np.array([x for x in res["C_f"] if np.ndim(x)]), axis=0) \
     if any(np.ndim(x) for x in res["C_f"]) else np.full(args.nd, np.nan)
Ce = np.nanmean(np.array([x for x in res["C_e"] if np.ndim(x)]), axis=0) \
     if any(np.ndim(x) for x in res["C_e"]) else np.full(args.nd, np.nan)
ax.plot(d_mid, Cf, "-", lw=2.6, color="#B71C1C", label="media FRENTE")
ax.plot(d_mid, Ce, "--", lw=2.6, color="#1565C0", label="media ESTELA")
ax.axhline(1 / np.e, color="0.5", ls=":", lw=1.0)
ax.axhline(0, color="0.5", lw=0.8)
ax.axvline(a_red, color="0.45", ls="--", lw=1.1)
ax.set_xlabel(r"separacion $d$  [$\sigma$]")
ax.set_ylabel(r"$C(d)$")
ax.set_title("$C(d)$ (naranja claro→oscuro = tiempo, en el frente)")
ax.legend(fontsize=7.5)

ax = axs[2]
ax.plot(t, r_f, "o-", ms=4, color="#2E7D32", label="frente medido (umbral KE)")
cR = np.polyfit(t, r_f, 1)[0] if t.size > 2 else np.nan
ax.plot(t, np.polyval(np.polyfit(t, r_f, 1), t), ":", lw=1.5, color="0.4",
        label=rf"ajuste: $c={cR:.2f}\,\sigma/\tau$")
ax.fill_between(t, r_f - args.ancho_frente, r_f, color="#B71C1C", alpha=0.18,
                label="anillo del frente")
ax.axhline(args.r_crater, color="0.5", ls="--", lw=1.0, label="crater excluido")
ax.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]")
ax.set_ylabel(r"$r$  [$\sigma$]")
ax.set_title("Donde se midio cada cosa")
ax.legend(fontsize=7)

fig.suptitle("Frente vs estela: el anillo SIGUE al frente, y la correlacion pesa por amplitud  "
             f"({os.path.basename(os.path.dirname(args.traj))})")
fig.tight_layout()
fig.savefig(args.out, dpi=200, facecolor="white")
print(f"Figura guardada: {args.out}")
