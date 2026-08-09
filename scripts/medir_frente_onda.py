"""
medir_frente_onda.py — ¿Balístico (r ~ t) o difusivo (r ~ sqrt(t))?

Mide la posición del frente de energía r_frente(t) en la ventana limpia de
la corrida sin periodicidad (ondas_N79_libre, Nx=Ny=80), antes de que el eco
del borde libre contamine la señal (~7.5 tau tras el impacto, ver
fig_ondas_superficie_libre.png). Ajusta r_frente = A*(t-t_imp)^p por mínimos
cuadrados sobre log-log; p~1 => ecuación de onda (balístico), p~0.5 =>
ecuación de calor (difusivo).

Uso:
    python3 scripts/medir_frente_onda.py \\
        --traj data/outputs/activo/ondas_N79_libre/traj.lammpstrj \\
        --bounce data/outputs/activo/ondas_N79_libre/bounce.dat \\
        --out data/outputs/activo/fig_frente_onda.png \\
        --t-clean 7.5
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
p.add_argument("--t-clean", type=float, default=7.5,
               help="ventana limpia tras el impacto [tau] (ver diagnóstico previo)")
p.add_argument("--thresh-frac", type=float, default=0.15,
               help="umbral del frente, como fracción del pico KE de cada instante")
p.add_argument("--dt", type=float, default=0.005)
p.add_argument("--nbins", type=int, default=60)
args = p.parse_args()

cx, cy, info = ld.centro_bola(args.traj)
Lx, Ly = info.Lx, info.Ly

t_imp = ld.t_impacto(args.bounce)
t_lo, t_hi = t_imp, t_imp + args.t_clean
print(f"t_imp={t_imp:.2f}  ventana limpia=[{t_lo:.2f},{t_hi:.2f}] tau")

times = info.steps * args.dt

r_max = min(Lx, Ly) / 2
r_edges = np.linspace(0, r_max, args.nbins + 1)
r_mid = 0.5 * (r_edges[:-1] + r_edges[1:])

# Pasada 2: solo se procesan los frames de la ventana limpia; el resto se lee
# y se descarta sin tocarlo.
print("Pasada 2/2: frente de energía en la ventana limpia ...")
t_front, r_front, r_rms = [], [], []
for j, fr in enumerate(ld.iter_frames(args.traj, ("type", "x", "y", "ke"))):
    if not (t_lo <= times[j] <= t_hi):
        continue
    c = fr.col
    surf = fr.data[fr.data[:, c["type"]] == 1]
    x, y, ke = surf[:, c["x"]], surf[:, c["y"]], surf[:, c["ke"]]
    dx = x - cx; dy = y - cy
    r = np.hypot(dx, dy)
    inr = r < r_max
    r_in, ke_in = r[inr], ke[inr]

    prof = np.nan_to_num(ld.media_por_bin(r, ke, r_edges), nan=0.0)
    if prof.max() <= 0:
        continue
    thresh = args.thresh_frac * prof.max()
    above = np.where(prof >= thresh)[0]
    if len(above) == 0:
        continue
    t_front.append(times[j] - t_imp)
    r_front.append(r_mid[above[-1]])

    # r_rms: radio cuadrático medio ponderado por KE excedente (por átomo, no
    # por bin -- evita el sesgo de área creciente: el área de un anillo crece
    # como r, así que un piso frío pero no-cero sesgaría <r^2> hacia afuera si
    # no se resta primero). Línea base = mediana (robusta a la cola caliente).
    baseline = np.median(ke_in)
    excess = np.clip(ke_in - baseline, 0.0, None)
    if excess.sum() > 0:
        r_rms.append(np.sqrt(np.sum(excess * r_in ** 2) / np.sum(excess)))
    else:
        r_rms.append(np.nan)
r_rms = np.array(r_rms)

t_front = np.array(t_front); r_front = np.array(r_front)
m = t_front > 0.3   # descarta el primer instante (r~0, log mal condicionado)
logt, logr = np.log(t_front[m]), np.log(r_front[m])
p_fit, logA = np.polyfit(logt, logr, 1)
A_fit = np.exp(logA)
print(f"\nAjuste r_frente = A * t^p:  p = {p_fit:.3f}   A = {A_fit:.3f}")
print("  (p~1.0 => balístico/ecuación de onda; p~0.5 => difusivo/ecuación de calor)")

# --- Exponente local p(t): ventana deslizante en log-log ----------------------
# Un solo p global solo dice "en promedio, balístico o difusivo". Lo que se
# necesita para ver el CRUCE (misma física, dos regímenes del mismo proceso)
# es si p(t) empieza a derivar de 1 hacia 1/2 con el tiempo -- eso es visible
# aun dentro de una ventana "corta" si el cruce ya empezó, sin necesitar
# ver el régimen difusivo puro todavía.
W = 7  # puntos por ventana deslizante
t_local, p_local = [], []
order = np.argsort(t_front[m])
lt, lr, tf = logt[order], logr[order], t_front[m][order]
for k in range(len(lt) - W):
    seg_t, seg_r = lt[k:k + W], lr[k:k + W]
    pk, _ = np.polyfit(seg_t, seg_r, 1)
    t_local.append(tf[k:k + W].mean())
    p_local.append(pk)
t_local, p_local = np.array(t_local), np.array(p_local)

# comparación: mejor ajuste forzando p=1 y p=0.5
A1 = np.exp(np.mean(logr - 1.0 * logt))
A05 = np.exp(np.mean(logr - 0.5 * logt))
resid1 = np.sum((logr - (logt + np.log(A1))) ** 2)
resid05 = np.sum((logr - (0.5 * logt + np.log(A05))) ** 2)
print(f"  residuo^2 forzando p=1.0 (onda):    {resid1:.4f}  (A={A1:.3f})")
print(f"  residuo^2 forzando p=0.5 (difusión): {resid05:.4f}  (A={A05:.3f})")

# --- Lo mismo con r_rms (ponderado por energía, menos ruidoso) ----------------
mr = (t_front > 0.3) & np.isfinite(r_rms) & (r_rms > 0)
logt_r, logr_r = np.log(t_front[mr]), np.log(r_rms[mr])
p_rms, logA_rms = np.polyfit(logt_r, logr_r, 1)
A_rms = np.exp(logA_rms)
print(f"\nCon r_rms (ponderado por energía):  p = {p_rms:.3f}   A = {A_rms:.3f}")

order_r = np.argsort(t_front[mr])
lt_r, lr_r, tf_r = logt_r[order_r], logr_r[order_r], t_front[mr][order_r]
t_local_rms, p_local_rms = [], []
for k in range(len(lt_r) - W):
    pk, _ = np.polyfit(lt_r[k:k + W], lr_r[k:k + W], 1)
    t_local_rms.append(tf_r[k:k + W].mean())
    p_local_rms.append(pk)
t_local_rms, p_local_rms = np.array(t_local_rms), np.array(p_local_rms)

# --- Figura -------------------------------------------------------------------
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15.5, 4.5))

ax1.plot(t_front, r_front, "o", ms=4, color="#1565C0", label="frente (umbral)")
ax1.plot(t_front[mr], r_rms[mr], "s", ms=4, color="#E65100", alpha=0.7,
          label=r"$r_{\rm rms}$ (ponderado energía)")
tt_fine = np.linspace(t_front.min(), t_front.max(), 100)
ax1.plot(tt_fine, A_fit * tt_fine ** p_fit, "-", color="#B71C1C",
          label=rf"ajuste frente: $p={p_fit:.2f}$")
ax1.plot(tt_fine, A_rms * tt_fine ** p_rms, "-", color="#E65100", lw=1.6,
          label=rf"ajuste $r_{{\rm rms}}$: $p={p_rms:.2f}$")
ax1.plot(tt_fine, A1 * tt_fine, "--", color="0.4", lw=1.2,
          label=r"forzado $p=1$ (onda)")
ax1.plot(tt_fine, A05 * tt_fine ** 0.5, ":", color="0.4", lw=1.6,
          label=r"forzado $p=1/2$ (calor)")
ax1.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]")
ax1.set_ylabel(r"$r$  [$\sigma$]")
ax1.set_title("Escala lineal")
ax1.legend(fontsize=7.5)

ax2.loglog(t_front, r_front, "o", ms=4, color="#1565C0")
ax2.loglog(t_front[mr], r_rms[mr], "s", ms=4, color="#E65100", alpha=0.7)
ax2.loglog(tt_fine, A_fit * tt_fine ** p_fit, "-", color="#B71C1C")
ax2.loglog(tt_fine, A_rms * tt_fine ** p_rms, "-", color="#E65100", lw=1.6)
ax2.loglog(tt_fine, A1 * tt_fine, "--", color="0.4", lw=1.2)
ax2.loglog(tt_fine, A05 * tt_fine ** 0.5, ":", color="0.4", lw=1.6)
ax2.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]  (log)")
ax2.set_ylabel(r"$r$  [$\sigma$]  (log)")
ax2.set_title(f"Escala log-log — frente p={p_fit:.2f}, rms p={p_rms:.2f}")

ax3.plot(t_local, p_local, "o-", ms=4, color="#1565C0", label="frente (umbral)")
ax3.plot(t_local_rms, p_local_rms, "s-", ms=4, color="#E65100",
          label=r"$r_{\rm rms}$ (ponderado energía)")
ax3.axhline(1.0, color="0.4", ls="--", lw=1.2, label=r"$p=1$ (onda)")
ax3.axhline(0.5, color="0.4", ls=":", lw=1.6, label=r"$p=1/2$ (calor)")
ax3.set_xlabel(r"$t-t_{\rm imp}$  [$\tau$]  (centro de ventana deslizante)")
ax3.set_ylabel(r"$p(t)$ local")
ax3.set_title(f"Exponente local (ventana de {W} puntos)")
ax3.legend(fontsize=7.5)
ax3.set_ylim(0, 1.6)

fig.suptitle(rf"Frente de energía en la superficie (umbral={args.thresh_frac:.0%} del pico) "
             rf"vs. $r_{{\rm rms}}$ (ponderado por energía) — ventana limpia {args.t_clean:.1f}τ")
fig.tight_layout()
fig.savefig(args.out, dpi=200, facecolor="white")
print(f"\nFigura guardada: {args.out}")
