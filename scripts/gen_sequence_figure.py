"""
gen_sequence_figure.py — Secuencia de 4 instantáneas del primer rebote (paper).

Sistema: cluster N=79 pre-minimizado, h0=28.15σ, c_ps=0 (in.bounce_cluster
con -var dump_freq 50). Los 4 eventos se detectan automáticamente de
bounce.dat:
  (a) caída libre      (t = t_imp − 8 τ)
  (b) compresión máxima (primer mínimo de Z_cm)
  (c) expansión         (t = t_imp + 6 τ)
  (d) primer apogeo     (primer máximo de Z_cm tras el impacto)

Coloración térmica (fondo claro):
  Bola:       KE interna por átomo (sin CM), escala 0–0.30 ε
  Superficie: KE por átomo, escala 0–0.01 ε (×30 más sensible)

Barra bajo cada panel: partición instantánea de E0 en
  K_CM, T_B, ΔV_B, ΔE_surf, U_c (transitoria) y U_F (gravitatoria).

Uso (desde la raíz del proyecto):
    LD_LIBRARY_PATH=~/ovito_env/lib/python3.12/site-packages/PySide6/Qt/lib \\
    ~/ovito_env/bin/python scripts/gen_sequence_figure.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image

import sys
LANG = "en" if "--en" in sys.argv else "es"

RUN     = "data/outputs/viz_N79_h28"
TRAJ    = f"{RUN}/traj.lammpstrj"
BOUNCE  = f"{RUN}/bounce.dat"
OUT     = ("data/outputs/fig_sequence.png" if LANG == "es"
           else "data/outputs/fig_sequence_en.png")
TMPDIR  = "/tmp/claude-1000/sequence_panels"

DT_DUMP = 0.25   # tau/frame (dump_freq=50, dt=0.005)
DT_PRNT = 0.10   # tau/linea de bounce.dat (print_freq=20)
G       = 0.026
M       = 79.0
W, H    = 900, 640
BG      = (0.93, 0.945, 0.96)
KE_MAX_BALL = 0.30
KE_MAX_SURF = 0.01

os.makedirs(TMPDIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Eventos desde bounce.dat
# ---------------------------------------------------------------------------
COLS = "time Zcm Vcmz KB TB VB Uc UF KB3d KS PE_surf".split()
d = np.loadtxt(BOUNCE, comments="#")
arr = {c: d[:, i] for i, c in enumerate(COLS[:d.shape[1]])}
tt, Zcm = arr["time"], arr["Zcm"]

# primer impacto = mínimo global (el más profundo); apogeo = máximo en los
# 40 τ siguientes (sin scipy: el venv de OVITO solo tiene numpy)
i_imp = int(np.argmin(Zcm))
i_apo = i_imp + int(np.argmax(Zcm[i_imp:i_imp + int(40 / DT_PRNT)]))
t_imp, t_apo = tt[i_imp], tt[i_apo]
z0 = float(Zcm[i_imp])
h0_eff = float(Zcm[0] - z0)
E0 = M * G * h0_eff

T_EVENTS = [t_imp - 5.0, t_imp, t_imp + 6.0, t_apo]
TLABELS  = [rf"$t={t:.0f}\,\tau$" for t in T_EVENTS]
DESCS    = {"es": ["Caída libre", "Compresión máxima", "Expansión", "Primer apogeo"],
            "en": ["Free fall", "Maximum compression", "Expansion", "First apex"]}[LANG]
PANLAB   = ["(a)", "(b)", "(c)", "(d)"]

print(f"t_imp={t_imp:.1f}τ  t_apo={t_apo:.1f}τ  z0={z0:.2f}σ  "
      f"h0_eff={h0_eff:.2f}σ  E0={E0:.1f}ε")

# Partición instantánea (fracciones de E0) en cada evento
def particion(t_ev):
    i = np.argmin(np.abs(tt - t_ev))
    KCM  = arr["KB3d"][i]
    TB   = arr["TB"][i]
    dVB  = arr["VB"][i] - arr["VB"][0]
    Uc   = arr["Uc"][i]
    # PE_surf cruda incluye 0.5*Uc (convención pe/atom)
    dEs  = (arr["KS"][i] + arr["PE_surf"][i] - 0.5 * Uc
            - arr["KS"][0] - arr["PE_surf"][0])
    UF   = max(M * G * (Zcm[i] - z0), 0.0)
    return np.array([KCM, TB, dVB, dEs, Uc, UF]) / E0

FRACS = [particion(t) for t in T_EVENTS]

# ---------------------------------------------------------------------------
# Render OVITO
# ---------------------------------------------------------------------------
from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer

_NODES = np.array([[0.05, 0.05, 0.40], [0.20, 0.55, 0.85],
                   [0.95, 0.45, 0.05], [1.00, 0.95, 0.10]])
_POS = np.array([0.0, 0.33, 0.70, 1.0])

def thermal_color(v):
    v = np.clip(v, 0.0, 1.0)
    out = np.zeros((len(v), 3))
    for k in range(len(_POS) - 1):
        m = (v >= _POS[k]) & (v < _POS[k + 1])
        s = (v[m] - _POS[k]) / (_POS[k + 1] - _POS[k])
        out[m] = (1 - s[:, None]) * _NODES[k] + s[:, None] * _NODES[k + 1]
    out[v >= 1.0] = _NODES[-1]
    return out

def color_thermal(frame, data):
    types = np.array(data.particles["Particle Type"])
    ke = np.array(data.particles["c_keratom"])
    vel = np.array(data.particles["Velocity"])
    colors = np.zeros((len(types), 3))
    surf = types == 1
    colors[surf] = thermal_color(ke[surf] / KE_MAX_SURF)
    ball = types == 2
    if ball.any():
        vb = vel[ball]
        ke_int = 0.5 * ((vb - vb.mean(axis=0)) ** 2).sum(axis=1)
        colors[ball] = thermal_color(ke_int / KE_MAX_BALL)
    data.particles_.create_property("Color", data=colors)

pipeline = import_file(TRAJ, sort_particles=True)
pipeline.source.data.particles_.vis.radius = 0.45
pipeline.source.data.cell_.vis.enabled = False
pipeline.modifiers.append(color_thermal)
pipeline.add_to_scene()

pos0 = np.array(pipeline.compute(0).particles["Position"])
xc, yc = pos0[:, 0].mean(), pos0[:, 1].mean()

renderer = TachyonRenderer()
renderer.ambient_occlusion = True
renderer.shadows = True
renderer.direct_light_intensity = 1.45
renderer.ambient_occlusion_brightness = 1.15

vp = Viewport(type=Viewport.Type.Perspective)
vp.camera_pos = (xc, yc - 70, 11.5)
vp.camera_dir = (0.0, 1.0, 0.0)
vp.camera_up = (0.0, 0.0, 1.0)
vp.fov = 0.33

print("Renderizando paneles...")
pngs = []
for i, t_ev in enumerate(T_EVENTS):
    fr = int(round(t_ev / DT_DUMP))
    out_png = os.path.join(TMPDIR, f"panel_{i}.png")
    if not os.path.exists(out_png) or "--rerender" in sys.argv:
        vp.render_image(size=(W, H), frame=fr, renderer=renderer,
                        background=BG, filename=out_png)
    pngs.append(out_png)
    print(f"  {PANLAB[i]} {DESCS[i]:18s} t={t_ev:6.1f}τ (frame {fr})  "
          + "  ".join(f"{n}={100*f:.0f}%" for n, f in
                      zip(["KCM", "TB", "dVB", "dEs", "Uc", "UF"], FRACS[i])))

# ---------------------------------------------------------------------------
# Composición
# ---------------------------------------------------------------------------
C_CM, C_TB, C_dVB = "#1565C0", "#B71C1C", "#E65100"
C_SURF, C_UC, C_UF = "#2E7D32", "#757575", "#F9A825"
BAR_COLS = [C_CM, C_TB, C_dVB, C_SURF, C_UC, C_UF]
LEG = {"es": [
    r"$K_{\rm CM}$ — KE traslacional del CM",
    r"$T_B$ — KE interna de la bola",
    r"$\Delta V_B$ — PE intra-bola (expansión anarmónica)",
    r"$\Delta E_{\rm surf}$ — energía de la superficie (fonones)",
    r"$U_c$ — PE de contacto (transitoria)",
    r"$U_F$ — PE gravitatoria restante",
], "en": [
    r"$K_{\rm CM}$ — CM translational KE",
    r"$T_B$ — internal KE of the ball",
    r"$\Delta V_B$ — intra-ball PE (anharmonic expansion)",
    r"$\Delta E_{\rm surf}$ — surface energy (phonons)",
    r"$U_c$ — contact PE (transient)",
    r"$U_F$ — remaining gravitational PE",
]}[LANG]

TCOLD = {"es": "0 (fría)", "en": "0 (cold)"}[LANG]
TBALL = {"es": "Bola: KE interna / átomo", "en": "Ball: internal KE / atom"}[LANG]
TSURF = {"es": r"Superficie: KE / átomo ($\times$30 más sensible)",
         "en": r"Surface: KE / atom ($\times$30 more sensitive)"}[LANG]

fig = plt.figure(figsize=(9.0, 9.5))
fig.patch.set_facecolor("white")

POS = [[0.015, 0.620, 0.475, 0.320], [0.510, 0.620, 0.475, 0.320],
       [0.015, 0.220, 0.475, 0.320], [0.510, 0.220, 0.475, 0.320]]

for i, p in enumerate(POS):
    ax = fig.add_axes(p)
    ax.imshow(np.array(Image.open(pngs[i])))
    ax.axis("off")
    ax.text(0.03, 0.97, PANLAB[i], transform=ax.transAxes, fontsize=15,
            fontweight="bold", va="top")
    ax.text(0.97, 0.97, f"{TLABELS[i]} — {DESCS[i]}", transform=ax.transAxes,
            fontsize=10.5, va="top", ha="right",
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=2.5))
    # barra de partición bajo la imagen
    x0, bar_y, bar_h, bar_w = 0.04, -0.075, 0.055, 0.92
    for frac, col in zip(FRACS[i], BAR_COLS):
        w = max(frac, 0.0) * bar_w
        if w < 0.003:
            continue
        ax.add_patch(mpatches.Rectangle(
            (x0, bar_y), w, bar_h, facecolor=col, edgecolor="white",
            lw=0.5, transform=ax.transAxes, clip_on=False, zorder=5))
        if frac >= 0.05:
            ax.text(x0 + w / 2, bar_y + bar_h / 2, f"{100*frac:.0f}%",
                    transform=ax.transAxes, fontsize=9, fontweight="bold",
                    ha="center", va="center", color="white", zorder=6,
                    clip_on=False)
        x0 += w

# leyenda (fila propia, sin solaparse con las colorbars)
patches = [mpatches.Patch(color=c, label=l) for c, l in zip(BAR_COLS, LEG)]
fig.legend(handles=patches, loc="center", ncol=2, fontsize=9,
           bbox_to_anchor=(0.5, 0.125), frameon=False,
           handlelength=1.4, columnspacing=1.2)

# dos colorbars en la franja inferior
cmap = LinearSegmentedColormap.from_list("thermal", _NODES.tolist())
cax_b = fig.add_axes([0.06, 0.038, 0.40, 0.013])
cb = ColorbarBase(cax_b, cmap=cmap, orientation="horizontal")
cb.set_ticks([0, 1])
cb.set_ticklabels([TCOLD, f"{KE_MAX_BALL:.2f} ε"])
cb.ax.tick_params(labelsize=8, pad=2)
cax_b.set_title(TBALL, fontsize=8.5, pad=3)

cax_s = fig.add_axes([0.54, 0.038, 0.40, 0.013])
cb = ColorbarBase(cax_s, cmap=cmap, orientation="horizontal")
cb.set_ticks([0, 1])
cb.set_ticklabels([TCOLD, f"{KE_MAX_SURF:.2f} ε"])
cb.ax.tick_params(labelsize=8, pad=2)
cax_s.set_title(TSURF, fontsize=8.5, pad=3)

fig.savefig(OUT, dpi=250, facecolor="white")
print(f"\nFigura guardada: {OUT} ({os.path.getsize(OUT)//1024} KB)")
