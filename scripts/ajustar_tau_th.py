"""
ajustar_tau_th.py — Post-procesa los .npz ya generados por
medir_coherencia_frente_estela.py (no vuelve a leer las trayectorias) y hace
un AJUSTE NUMERICO real de tau_th por altura, en vez de leer la tabla a ojo.

Reconstruye, para cada uno de los 3 anillos fijos (r_fijo), la "edad desde
que paso el frente" = t - t_pasada(r_fijo), con t_pasada de la recta
r_f = pend*t + inter ya ajustada dentro del script original. Junta los
puntos de los 3 anillos (edad >= 0) y ajusta

    xi(edad) = piso + (xi0 - piso) * exp(-edad / tau_th)

con piso FIJO en a_red (el "1er vecino FCC" = calor, ya reportado por el
script original) -- no es un parametro libre, es la misma physical floor
que ya se usa para leer la "semivida" a ojo.

v_imp se lee directo de bounce.dat (Vcmz en el frame mas cercano a t_imp),
no de ninguna tabla vieja.

Uso:
    python3 scripts/ajustar_tau_th.py
"""
import glob
import numpy as np
from scipy.optimize import curve_fit

ALTURAS = ["28.15", "56.93", "79.03"]
BASE = "data/outputs/activo/ondas_N79_profunda_h{h}_zfull"


def v_impacto(bounce_path, t_imp):
    # t_imp (leer_dump.t_impacto) es el MINIMO de Zcm -- maxima compresion,
    # con Vcmz ~0 ahi (la bola esta momentaneamente detenida, como una
    # pelota contra un resorte al fondo del rebote). La velocidad de
    # impacto real es el maximo de |Vcmz| durante la caida libre, ANTES
    # de que el contacto empiece a frenarla -- justo antes de t_imp.
    data = np.loadtxt(bounce_path, comments="#")
    t, vcmz = data[:, 0], data[:, 2]
    antes = t <= t_imp
    return float(np.max(np.abs(vcmz[antes])))


def decaimiento(edad, xi0, tau):
    piso = decaimiento.piso
    return piso + (xi0 - piso) * np.exp(-edad / tau)


resultados = []
for h in ALTURAS:
    outdir = BASE.format(h=h)
    npz_path = glob.glob(f"{outdir}/fig_coherencia_frente_estela.npz")[0]
    d = np.load(npz_path)

    t, XI_X, r_fijo = d["t"], d["XI_X"], d["r_fijo"]
    pend, inter, a_red, t_imp = float(d["pend"]), float(d["inter"]), float(d["a_red"]), float(d["t_imp"])

    edades, xis = [], []
    for i, r in enumerate(r_fijo):
        t_pasada = (r - inter) / pend
        edad_i = t - t_pasada
        xi_i = XI_X[i]
        ok = np.isfinite(xi_i) & (edad_i >= 0)
        edades.append(edad_i[ok])
        xis.append(xi_i[ok])
    edad = np.concatenate(edades)
    xi = np.concatenate(xis)

    decaimiento.piso = a_red
    xi0_guess = float(np.max(xi))
    tau_guess = 2.0
    popt, pcov = curve_fit(decaimiento, edad, xi, p0=[xi0_guess, tau_guess],
                            bounds=([a_red, 0.05], [50, 50]))
    xi0, tau_th = popt
    xi0_err, tau_err = np.sqrt(np.diag(pcov))

    v_imp = v_impacto(f"{outdir}/bounce.dat", t_imp)

    # bondad de ajuste
    pred = decaimiento(edad, *popt)
    ss_res = np.sum((xi - pred) ** 2)
    ss_tot = np.sum((xi - np.mean(xi)) ** 2)
    r2 = 1 - ss_res / ss_tot

    resultados.append(dict(h=h, v_imp=v_imp, xi0=xi0, xi0_err=xi0_err,
                            tau_th=tau_th, tau_err=tau_err, n=len(edad), r2=r2,
                            piso=a_red))

print(f"{'h0':>7} {'v_imp':>8} {'xi0 (σ)':>14} {'tau_th (τ)':>16} {'n':>4} {'R²':>6}")
print("-" * 62)
for r in resultados:
    print(f"{r['h']:>7} {r['v_imp']:8.3f} "
          f"{r['xi0']:6.2f} ± {r['xi0_err']:4.2f}  "
          f"{r['tau_th']:6.2f} ± {r['tau_err']:4.2f}   "
          f"{r['n']:4d} {r['r2']:6.3f}")

# Tendencia: pendiente de tau_th vs v_imp, con su propio error (regresion
# ponderada simple 1/sigma^2), para no leerla a ojo tampoco.
v = np.array([r["v_imp"] for r in resultados])
tau = np.array([r["tau_th"] for r in resultados])
tau_e = np.array([r["tau_err"] for r in resultados])
w = 1 / tau_e**2
A = np.vstack([v, np.ones_like(v)]).T
Wm = np.diag(w)
cov_slope = np.linalg.inv(A.T @ Wm @ A)
coef = cov_slope @ A.T @ Wm @ tau
slope, intercept = coef
slope_err = np.sqrt(cov_slope[0, 0])
print()
print(f"tau_th(v) ~ ({slope:+.3f} +/- {slope_err:.3f}) * v + {intercept:.2f}")
print(f"  pendiente / error = {slope/slope_err:+.2f}  "
      f"(|.|<2 => no distinguible de pendiente cero con 3 puntos)")
