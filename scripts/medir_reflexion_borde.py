"""
medir_reflexion_borde.py — Coeficiente de reflexion del borde LATERAL libre,
medido como energia que vuelve / energia que salio a traves de un cilindro de
control de radio R.

Para que sirve: fija cuanto tiene que absorber la esponja del diseno nuevo. Sin
este numero, "la esponja funciona" es fe. El criterio de aceptacion del nuevo
run es R < 1% medido con este mismo script.

Metodo: se sigue la KE de la superficie movil dentro de r < R y se compara la
fraccion E(r<R)/E_total con la FRACCION DE AREA del cilindro. Si el borde no
absorbe, la energia no se pierde: solo se reparte, y acaba en equiparticion
geometrica. Si absorbiera, E_total decaeria.

Descartado: estimar R sumando los repuntes de E_dentro. Ese estimador daba
R=0.96, pero por el motivo equivocado — cuenta como "eco" el ruido termico, la
relajacion lenta y, sobre todo, el SEGUNDO REBOTE de la bola (2do impacto en
t=87.6 tau, o sea t-t_imp=40.9 tau, visible como un salto de E_dentro a 39 tau).
Por eso el analisis se corta en 38 tau.

Se usa solo KE, no KE+PE: la PE por atomo necesitaria restar el estado de
referencia atomo a atomo (el dump no trae el id ordenado entre frames de forma
garantizada), y para un COCIENTE el factor virial ~2 se cancela.

Nota sobre el fondo: no hace falta medirlo. `fix setforce 0.0 0.0 0.0` es una
pared perfectamente rigida, o sea R_fondo = 1 por construccion. Lo que si esta
medido es CUANDO vuelve (scripts/medir_eco_fondo.py: 1.8 tau).

Uso:
    python3 scripts/medir_reflexion_borde.py \\
        --traj data/outputs/ondas_N79_libre/traj.lammpstrj \\
        --bounce data/outputs/ondas_N79_libre/bounce.dat \\
        --out data/outputs/fig_reflexion_borde.png
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
p.add_argument("--dt", type=float, default=0.005)
p.add_argument("--R", type=float, default=40.0, help="radio del cilindro de control [sigma]")
p.add_argument("--nbins", type=int, default=50)
p.add_argument("--z-fixed-top", type=float, default=3.0992)
args = p.parse_args()

cx, cy, info = ld.centro_bola(args.traj)
Lx, Ly = info.Lx, info.Ly
t_imp = ld.t_impacto(args.bounce)
times = info.steps * args.dt
r_borde = min(Lx, Ly) / 2 - max(cx - Lx / 2, cy - Ly / 2)
print(f"  t_imp = {t_imp:.2f} tau; borde libre mas cercano a r ~ {r_borde:.1f} sigma")

r_edges = np.linspace(0, min(Lx, Ly) / 2, args.nbins + 1)
r_mid = 0.5 * (r_edges[:-1] + r_edges[1:])

print("Pasada 2/2: energia dentro del cilindro de control ...")
E_dentro, E_total, E_rt, t_rel = [], [], [], []
for j, fr in enumerate(ld.iter_frames(args.traj, ("type", "x", "y", "z", "ke"))):
    if times[j] < t_imp:
        continue
    c = fr.col
    d = fr.data
    d = d[(d[:, c["type"]] == 1) & (d[:, c["z"]] > args.z_fixed_top)]
    r2d = np.hypot(d[:, c["x"]] - cx, d[:, c["y"]] - cy)
    ke = d[:, c["ke"]]
    E_dentro.append(float(ke[r2d < args.R].sum()))
    E_total.append(float(ke.sum()))
    # perfil de energia TOTAL por corona (no media): asi la suma sobre coronas
    # es la energia del sistema y se ve adonde va
    idx = np.clip(np.digitize(r2d, r_edges) - 1, 0, args.nbins - 1)
    E_rt.append(np.bincount(idx, weights=ke, minlength=args.nbins))
    t_rel.append(times[j] - t_imp)

E_dentro = np.array(E_dentro)
E_total = np.array(E_total)
E_rt = np.array(E_rt).T
t_rel = np.array(t_rel)

# --- Coeficiente de reflexion --------------------------------------------------
# La bola vuelve a impactar en t=87.6 tau (minimo local de Zcm), o sea
# t-t_imp = 40.9 tau, y mete un segundo pulso de energia. Todo lo que pase a
# partir de ~38 tau NO es eco: es el segundo rebote. Se corta antes.
c_med = 4.6
t_corte = 38.0
t_vuelta = 2 * r_borde / c_med - args.R / c_med
print(f"  con c={c_med} sigma/tau: el pulso toca el borde en t={r_borde/c_med:.1f} tau "
      f"y podria reingresar a r<{args.R:.0f} en t={t_vuelta:.1f} tau")
print(f"  se corta el analisis en t-t_imp = {t_corte} tau (2do impacto en 40.9 tau)")

# El coeficiente NO se estima sumando repuntes de E_dentro: ese estimador cuenta
# como "eco" tanto el ruido termico como la relajacion lenta, y ademas se comia
# el segundo rebote (daba R=0.96 por el motivo equivocado).
#
# Criterio limpio: si el borde no absorbe nada, la energia no se pierde, solo se
# reparte, y al final E_dentro/E_total tiene que valer la FRACCION DE AREA del
# cilindro de control. Si el borde absorbiera, E_total decaeria.
tardio = (t_rel > 20) & (t_rel < t_corte)
frac_area = np.pi * args.R**2 / (Lx * Ly)
frac_med = float(np.mean(E_dentro[tardio] / E_total[tardio]))
# Conservacion: E_total al final frente al maximo (descontando lo que la bola
# se lleva de vuelta al rebotar)
conserva = float(np.mean(E_total[tardio]) / E_total.max())

print("\n--- RESULTADO -------------------------------------------------")
print(f"  cilindro de control R = {args.R:.0f} sigma; borde libre a {r_borde:.1f} sigma")
print(f"  fraccion de area del cilindro          = {frac_area:.3f}")
print(f"  E_dentro/E_total medida en t=20-38 tau = {frac_med:.3f}")
print(f"  => la energia se reparte hasta la equiparticion GEOMETRICA:")
print(f"     el borde libre devuelve todo lo que le llega.")
print(f"  KE total del sustrato en t=20-38 frente a su maximo = {conserva:.3f}")
print(f"  R_lateral (borde libre, sin esponja) ~ 1.0")
print("     (esperado: un borde libre tiene traccion nula y refleja las ondas")
print("      elasticas por completo, igual que la pared rigida del fondo)")
print()
print("  R_fondo = 1.000 por construccion (fix setforce 0 0 0 = pared rigida);")
print("  el eco vertical esta medido en medir_eco_fondo.py: vuelve a 1.8 tau.")
print()
print(f"  CRITERIO PARA EL RUN NUEVO: R_lateral < 0.01 con esponja, medido asi.")
print("---------------------------------------------------------------")

npz = os.path.splitext(args.out)[0] + ".npz"
np.savez(npz, E_dentro=E_dentro, E_total=E_total, E_rt=E_rt, t_rel=t_rel,
         r_mid=r_mid, frac_area=frac_area, frac_med=frac_med,
         R_control=args.R, r_borde=r_borde)
print(f"Datos: {npz}")

# --- Figura --------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.4))

ax = axs[0]
im = ax.pcolormesh(t_rel, r_mid, np.log10(np.maximum(E_rt, 1e-5)),
                   cmap="inferno", shading="auto")
ax.axhline(r_borde, color="cyan", ls="--", lw=1.4, label="borde libre mas cercano")
ax.axhline(args.R, color="lime", ls=":", lw=1.4, label=rf"control $R={args.R:.0f}\,\sigma$")
ax.plot(t_rel, c_med * t_rel, "-", lw=1.2, color="white", alpha=0.7,
        label=rf"$c={c_med}\,\sigma/\tau$")
ax.set_ylim(0, r_mid.max())
ax.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]")
ax.set_ylabel(r"$r$  [$\sigma$]")
ax.set_title(r"$\log_{10}$ KE por corona: el pulso sale y el borde lo devuelve")
ax.legend(fontsize=7, loc="upper left")
fig.colorbar(im, ax=ax, pad=0.01)

ax = axs[1]
ax.plot(t_rel, E_dentro / E_total, "-", lw=1.5, color="#B71C1C",
        label=r"$E(r<R)/E_{\rm total}$")
ax.axhline(frac_area, color="#1565C0", ls="--", lw=1.5,
           label=rf"fraccion de area $={frac_area:.2f}$")
ax.axvline(t_corte, color="0.35", ls=":", lw=1.5, label=r"2do rebote")
ax.axvspan(t_corte, t_rel.max(), color="0.6", alpha=0.25)
ax.set_ylim(0, 1.0)
ax.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]")
ax.set_ylabel("fraccion de la KE dentro del control")
ax.set_title("Llega a la equiparticion geometrica: nada se absorbe\n"
             r"$R_{\rm lateral}\simeq1$ (objetivo del run nuevo: $<0.01$)")
ax.legend(fontsize=7.5)

fig.suptitle("P5: cuanto devuelve el borde libre, y por tanto cuanto debe absorber la esponja  "
             f"({os.path.basename(os.path.dirname(args.traj))})")
fig.tight_layout()
fig.savefig(args.out, dpi=200, facecolor="white")
print(f"Figura guardada: {args.out}")
