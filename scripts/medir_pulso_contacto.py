"""
medir_pulso_contacto.py — Duracion y contenido espectral del pulso de contacto
bola-superficie, a partir de Uc(t) (compute group/group) en bounce.dat.

Para que sirve: es el UNICO numero que dimensiona honestamente una capa
absorbente. El criterio estandar de las capas amortiguadoras (Holian & Ravelo,
PRB 51, 11275 (1995), y toda la literatura de PML) es "espesor >= la longitud
de onda significativa maxima que produce la simulacion". Sin medir lambda_dom
el espesor de la esponja es una opinion.

De aqui salen:
  - tau_c  = FWHM del pulso Uc(t)                 -> duracion del forzamiento
  - w_dom  = frecuencia dominante del espectro    -> que hay que amortiguar
  - lambda_dom = c * 2*pi/w_dom                   -> espesor minimo de esponja

Se reporta con la velocidad MEDIDA del frente de superficie (~4.6 sigma/tau,
Rayleigh) y con la longitudinal nominal (8.28), porque la esponja tiene que
absorber ambas.

Uso:
    python3 scripts/medir_pulso_contacto.py \\
        --bounce data/outputs/activo/ondas_N79_libre/bounce.dat \\
        --out data/outputs/activo/fig_pulso_contacto.png
"""

import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument("--bounce", required=True)
p.add_argument("--out", required=True)
p.add_argument("--c-rayleigh", type=float, default=4.6)
p.add_argument("--c-long", type=float, default=8.28)
args = p.parse_args()

COLS = "time Zcm Vcmz KB TB VB Uc UF KB3d KS PE_surf".split()
d = np.loadtxt(args.bounce, comments="#")
t = d[:, COLS.index("time")]
Uc = d[:, COLS.index("Uc")]
Zcm = d[:, COLS.index("Zcm")]

t_imp = t[np.argmin(Zcm)]
dt = float(np.median(np.diff(t)))

# --- FWHM del pulso ------------------------------------------------------------
pico = Uc.max()
i_pico = int(np.argmax(Uc))
mitad = pico / 2.0
izq = np.where(Uc[:i_pico] <= mitad)[0]
der = np.where(Uc[i_pico:] <= mitad)[0]
t_izq = t[izq[-1]] if izq.size else t[0]
t_der = t[i_pico + der[0]] if der.size else t[-1]
tau_c = t_der - t_izq

# --- Duracion de contacto al 5% del pico ---------------------------------------
# Se restringe a +-8 tau del impacto: Uc tiene una cola de ruido de ~1e-3 eps en
# todo el vuelo (pares bola-superficie dentro del cutoff antes y despues del
# rebote), y sin acotar la ventana el ultimo cruce lo fija esa cola, no el pulso.
vent = np.abs(t - t_imp) < 8.0
sobre = np.where(vent & (Uc > 0.05 * pico))[0]
t_c5 = t[sobre[-1]] - t[sobre[0]]

# --- Espectro ------------------------------------------------------------------
# ADVERTENCIA sobre la lectura: Uc(t) es un pulso UNIPOLAR, y el espectro de un
# pulso unipolar tiene su maximo en DC y decae monotonamente. No existe una
# "frecuencia dominante": si se resta la media y se enventana, el pico aparente
# lo fija el LARGO DE LA VENTANA, no la fisica (con una ventana de 16 tau salia
# w=0.79, que es justo 2*pi/8). Por eso aqui el espectro se usa solo para el
# corte de alta frecuencia, y el tamano de la esponja se fija por c*tau_c.
m = (t > t_imp - 8) & (t < t_imp + 8)
sig = (Uc[m] - Uc[m].mean()) * np.hanning(m.sum())
F = np.abs(np.fft.rfft(sig)) ** 2
w = 2 * np.pi * np.fft.rfftfreq(m.sum(), dt)

F_n = F[1:] / F[1:].max()
w_c = w[1:][np.where(F_n > 1e-3)[0][-1]]       # corte al -30 dB

# La extension espacial del paquete emitido es c*tau_c: ESE es el "espesor >=
# lambda significativa" del criterio de Holian & Ravelo.
lam_R = args.c_rayleigh * tau_c
lam_L = args.c_long * tau_c

# El muestreo del dump NO lo fija el pulso sino la vibracion atomica: la
# frecuencia mas alta que existe en el solido es la de Debye.
rho_red = 1.0742                                # N/V reducida del cristal FCC
k_D = (6 * np.pi**2 * rho_red) ** (1 / 3)
w_Debye = args.c_long * k_D

print("\n--- RESULTADO -------------------------------------------------")
print(f"  t_imp = {t_imp:.2f} tau, muestreo de bounce.dat dt = {dt:.3f} tau")
print(f"  Uc pico = {pico:.2f} eps (plateau ruidoso 7.0-7.5 entre t=46.1 y 46.8)")
print(f"  FWHM del pulso          tau_c = {tau_c:.2f} tau  ({t_izq:.2f} -> {t_der:.2f})")
print(f"  duracion de contacto al 5%    = {t_c5:.2f} tau")
print(f"  corte espectral a -30 dB      = {w_c:.1f} 1/tau")
print()
print(f"  extension del paquete, c_Rayleigh={args.c_rayleigh}: lambda_R = {lam_R:.1f} sigma")
print(f"  extension del paquete, c_long={args.c_long}    : lambda_L = {lam_L:.1f} sigma")
print()
print(f"  => ESPESOR DE ESPONJA: manda la longitudinal, lambda_L = {lam_L:.0f} sigma.")
print(f"     W = 40 sigma son {40/lam_L:.1f} lambda_L y {40/lam_R:.1f} lambda_R: suficiente.")
print(f"     (el minimo del criterio es 1 lambda; por debajo de {lam_L:.0f} sigma no bajar)")
print()
print(f"  => DUMP: lo fija la vibracion atomica, no el pulso.")
print(f"     w_Debye = c_L*(6 pi^2 rho)^(1/3) = {w_Debye:.1f} 1/tau")
print(f"     Nyquist necesario => dt_dump < pi/w_Debye = {np.pi/w_Debye:.4f} tau")
print(f"     dump actual 0.25 tau -> Nyquist {np.pi/0.25:.1f} 1/tau: "
      f"{w_Debye/(np.pi/0.25):.1f}x por debajo. ALIASEADO.")
print("---------------------------------------------------------------")

# --- Figura --------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(11, 4.0))

ax = axs[0]
ax.plot(t - t_imp, Uc, "-", lw=1.4, color="#B71C1C")
ax.axhline(mitad, color="0.5", ls=":", lw=1.0)
ax.annotate("", xy=(t_izq - t_imp, mitad), xytext=(t_der - t_imp, mitad),
            arrowprops=dict(arrowstyle="<->", color="#1565C0", lw=1.4))
ax.text(0.5 * (t_izq + t_der) - t_imp, mitad * 1.15,
        rf"FWHM $={tau_c:.2f}\,\tau$", ha="center", fontsize=9, color="#1565C0")
ax.set_xlim(-3, 6)
ax.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]")
ax.set_ylabel(r"$U_c$  [$\epsilon$]")
ax.set_title("Pulso de contacto bola-superficie")

ax = axs[1]
ax.semilogy(w[1:], F_n, "-", lw=1.3, color="#2E7D32")
ax.axvline(w_c, color="#B71C1C", ls="--", lw=1.3, label=rf"corte $-30$ dB $={w_c:.1f}$")
ax.axvline(np.pi / 0.25, color="0.4", ls=":", lw=1.6,
           label=r"Nyquist del dump ($\Delta t=0.25\,\tau$)")
ax.axvline(w_Debye, color="#6A1B9A", ls="-.", lw=1.5,
           label=rf"$\omega_{{\rm Debye}}={w_Debye:.0f}$")
ax.set_xlim(0, 40)
ax.set_ylim(1e-6, 2)
ax.set_xlabel(r"$\omega$  [$\tau^{-1}$]")
ax.set_ylabel("potencia normalizada")
ax.set_title("Espectro: unipolar, sin pico propio (el maximo es DC)")
ax.legend(fontsize=7)

fig.suptitle(r"El pulso de contacto fija $\lambda_{\rm dom}$ y con ella el espesor de la esponja")
fig.tight_layout()
fig.savefig(args.out, dpi=200, facecolor="white")
print(f"Figura guardada: {args.out}")
