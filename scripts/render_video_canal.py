"""
render_video_canal.py — Renderiza UN video por canal de energía, con la
ecuación del canal y su porcentaje de E0 actualizándose en vivo.

Canales (definiciones de la Sec. 2 del manuscrito):
  cm    alpha_cm(t)   = [K_cm(t) + Mg(z_cm(t)-z0)] / E0   -> h1/h0 al apogeo
  TB    alpha_TB(t)   = T_B(t) / E0                        (bola, KE interna)
  VB    alpha_VB(t)   = [V_B(t)-V_B(0)] / E0                (bola, PE interna)
  surf  alpha_surf(t) = [Delta(K_S+V_surf)] / E0            (superficie)

Coloreado por átomo, específico de cada canal:
  TB    bola: KE interna (relativa a la V_cm instantánea de la bola), Viridis.
              superficie: gris (no pertenece al canal).
  VB    bola: PE relativa a t=0 por átomo, azul-blanco-rojo.
              superficie: gris.
  surf  superficie: energía total por átomo relativa a t=0 (KE + DeltaPE),
              azul-blanco-rojo. bola: gris.
  cm    canal colectivo (del CM, no por átomo): bola en color sólido de
              realce, superficie gris.

Uso:
    LD_LIBRARY_PATH=~/ovito_env/lib/python3.12/site-packages/PySide6/Qt/lib:$LD_LIBRARY_PATH \\
    ~/ovito_env/bin/python scripts/render_video_canal.py \\
        --traj data/outputs/viz_N79_h28/traj.lammpstrj \\
        --bounce data/outputs/viz_N79_h28/bounce.dat \\
        --channel TB --N 79 \\
        --out docs/reunion_gonzalo/videos_reunion/canal_TB.mp4 \\
        --t0 30 --t1 95 --vmax 0.9
"""

import argparse
import os
import shutil
import subprocess
import sys

import numpy as np
from scipy.signal import find_peaks

p = argparse.ArgumentParser()
p.add_argument("--traj", required=True)
p.add_argument("--bounce", required=True, help="bounce.dat de la misma corrida")
p.add_argument("--channel", required=True, choices=["cm", "TB", "VB", "surf"])
p.add_argument("--out", required=True)
p.add_argument("--N", type=int, required=True, help="átomos de la bola")
p.add_argument("--t0", type=float, required=True)
p.add_argument("--t1", type=float, required=True)
p.add_argument("--vmax", type=float, default=0.6,
               help="escala de color del canal (unidades de energía por átomo, eps)")
p.add_argument("--dt", type=float, default=0.005)
p.add_argument("--fps", type=int, default=25)
p.add_argument("--panel", default="900x680")
args = p.parse_args()

G = 0.026
PW, PH = (int(x) for x in args.panel.split("x"))
HEADER_H, FOOTER_H, PLOT_H = 138, 116, 130
FRAMES_DIR = os.path.join(
    os.environ.get("TMPDIR", "/tmp"),
    "canal_frames_" + os.path.splitext(os.path.basename(args.out))[0])
RAW = os.path.join(FRAMES_DIR, "raw")
shutil.rmtree(FRAMES_DIR, ignore_errors=True)
os.makedirs(RAW)
os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

# --- Serie temporal del canal desde bounce.dat ------------------------------
d = np.loadtxt(args.bounce, comments="#")
t_b, Zcm, KB, TB, VB, Uc, UF, KB3d, KS, PEs = (
    d[:, 0], d[:, 1], d[:, 3], d[:, 4], d[:, 5], d[:, 6], d[:, 7], d[:, 8],
    d[:, 9], d[:, 10])
mins, _ = find_peaks(-Zcm, prominence=0.5, distance=50)
z0 = float(Zcm[mins].mean()) if len(mins) else float(Zcm.min())
h0_eff = Zcm[0] - z0
E0 = args.N * G * h0_eff
v_imp = float(np.sqrt(2.0 * G * h0_eff))

if args.channel == "cm":
    pct = (KB3d + args.N * G * (Zcm - z0)) / E0 * 100.0
    eq = (r"$\alpha_{\rm cm}(t)=\dfrac{K_{\rm cm}(t)+Mg\,(z_{\rm cm}(t)-z_0)}"
          r"{E_0}\ \longrightarrow\ \dfrac{h_1}{h_0}$")
    title = r"Canal $\alpha_{\rm cm}$ — energía mecánica del centro de masa"
    detail = ("K_cm(t) = ½ M (Vcmx² + Vcmy² + Vcmz²)  →  variable KB3d en bounce.dat\n"
               "M = N·m (N={0}, m=1);  z₀ = nivel de contacto = mín(z_cm) tras el impacto"
               .format(args.N))
elif args.channel == "TB":
    pct = TB / E0 * 100.0
    eq = r"$\alpha_{T_B}(t)=\dfrac{T_B(t)}{E_0}$"
    title = r"Canal $\alpha_{T_B}$ — energía cinética interna de la bola"
    detail = ("T_B(t) = Σ_(i∈bola) ½ m (v_i − v_cm,bola)²  =  compute ball_ke − KB3d\n"
               "KE interna de la bola: resta la traslación 3D del CM a la KE total de sus átomos")
elif args.channel == "VB":
    pct = (VB - VB[0]) / E0 * 100.0
    eq = r"$\alpha_{\Delta V_B}(t)=\dfrac{V_B(t)-V_B(0)}{E_0}$"
    title = r"Canal $\alpha_{\Delta V_B}$ — energía potencial interna de la bola"
    detail = ("V_B(t) = Σ_(i∈bola) PE_i(t) − ½ U_c(t)  =  compute pe_ball − ½ compute Uc\n"
               "U_c: interacción bola–superficie (group/group); se resta la mitad porque\n"
               "pe/atom la reparte entre los dos átomos de cada par cruzado")
else:
    pct = (KS + PEs - KS[0] - PEs[0]) / E0 * 100.0
    eq = r"$\alpha_{\rm surf}(t)=\dfrac{\Delta(K_S+V_{\rm surf})}{E_0}$"
    title = r"Canal $\alpha_{\rm surf}$ — fonones de la superficie"
    detail = ("K_S(t) = Σ_(i∈surf) ½ m v_i²  (compute surf_ke);  "
               "V_surf(t) = Σ_(i∈surf) PE_i  (compute pe_surf_r)\n"
               "V_surf incluye la mitad de U_c asignada a la superficie; α_surf usa Δ(K_S+V_surf) vs t=0")

DETAIL_LINES = detail.count("\n") + 1
DETAIL_H = 14 + 22 * DETAIL_LINES
TITLE_H, EQ_H = 32, 56
HEAD_MARGIN = 4
HEADER_H = HEAD_MARGIN + TITLE_H + EQ_H + DETAIL_H


def pct_at(tau):
    return float(np.interp(tau, t_b, pct))


# --- OVITO: pipeline + coloreado específico del canal -----------------------
from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer
from ovito.modifiers import ColorCodingModifier

print(f"Cargando {args.traj} ...")
pipe = import_file(args.traj, sort_particles=True)
nframes = pipe.source.num_frames
times = np.array([pipe.compute(f).attributes["Timestep"] * args.dt for f in (0, 1)])
dt_dump = times[1] - times[0]
frames = [f for f in range(nframes)
          if args.t0 <= f * dt_dump + times[0] <= args.t1]
print(f"  {len(frames)} frames ({len(frames)/args.fps:.1f} s)")

d0 = pipe.compute(0)
id0 = np.asarray(d0.particles["Particle Identifier"])
pe0_map = dict(zip(id0.tolist(), np.asarray(d0.particles["c_pe_atom"]).tolist()))

GRAY = np.array([0.80, 0.80, 0.83])
HIGHLIGHT = np.array([0.20, 0.45, 0.78])   # bola en el canal alpha_cm


def compute_channel(frame, data):
    typ = np.asarray(data.particles["Particle Type"])
    ids = np.asarray(data.particles["Particle Identifier"])
    is_ball = typ == 2
    is_surf = typ == 1
    if args.channel == "TB":
        vel = np.asarray(data.particles["Velocity"])
        v_cm = vel[is_ball].mean(axis=0) if is_ball.any() else np.zeros(3)
        v_int = vel.copy(); v_int[is_ball] -= v_cm
        chan = 0.5 * np.sum(v_int**2, axis=1)
        active = is_ball
    elif args.channel == "VB":
        pe = np.asarray(data.particles["c_pe_atom"])
        pe0 = np.array([pe0_map[i] for i in ids])
        chan = pe - pe0
        active = is_ball
    elif args.channel == "surf":
        ke = np.asarray(data.particles["c_keratom"])
        pe = np.asarray(data.particles["c_pe_atom"])
        pe0 = np.array([pe0_map[i] for i in ids])
        chan = ke + (pe - pe0)
        active = is_surf
    else:
        chan = np.zeros(len(typ))
        active = np.zeros(len(typ), dtype=bool)
    data.particles_.create_property("Chan", data=chan)
    data.particles_.create_property("Active", data=active.astype(float))


def flatten_inactive(frame, data):
    active = np.asarray(data.particles["Active"])
    typ = np.asarray(data.particles["Particle Type"])
    if "Color" in data.particles:
        color = np.asarray(data.particles["Color"]).copy()
    else:
        color = np.tile(GRAY, (len(typ), 1))
    if args.channel == "cm":
        color[typ == 2] = HIGHLIGHT
        color[typ == 1] = GRAY
    else:
        color[active < 0.5] = GRAY
    data.particles_.create_property("Color", data=color)


pipe.modifiers.append(compute_channel)
if args.channel != "cm":
    gradient = (ColorCodingModifier.Viridis() if args.channel == "TB"
                else ColorCodingModifier.BlueWhiteRed())
    start = 0.0 if args.channel == "TB" else -args.vmax
    pipe.modifiers.append(ColorCodingModifier(
        property="Chan", start_value=start, end_value=args.vmax,
        gradient=gradient, only_selected=False))
pipe.modifiers.append(flatten_inactive)
pipe.source.data.particles_.vis.radius = 0.45
pipe.source.data.cell_.vis.enabled = False

BG = (0.97, 0.97, 0.985)
data_imp = pipe.compute(frames[len(frames) // 2])
pos = np.array(data_imp.particles["Position"])
xc, yc, _ = pos.mean(axis=0)
zc = 12.0

vp = Viewport()
vp.type = Viewport.Type.Perspective
vp.camera_pos = (xc, yc - 45, zc + 8)
vp.camera_dir = (0.0, 1.0, -0.18)
vp.camera_up = (0.0, 0.0, 1.0)
vp.fov = 0.65

renderer = TachyonRenderer()
renderer.ambient_occlusion = True
renderer.shadows = True
renderer.direct_light_intensity = 1.4
renderer.ambient_occlusion_brightness = 1.0

pipe.add_to_scene()
print("Renderizando panel ...")
for i, fr in enumerate(frames):
    vp.render_image(size=(PW, PH), frame=fr, renderer=renderer,
                     background=BG,
                     filename=os.path.join(RAW, f"f_{i:04d}.png"))
    if i % 50 == 0 or i == len(frames) - 1:
        print(f"  [{100*i//len(frames):3d}%] frame {fr:4d}", flush=True)
pipe.remove_from_scene()

# --- Cabecera (título + ecuación) y barra de color, generadas una vez ------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

DARK = "#1e1e28"


def render_mathtext(text, fontsize, w_px, h_px):
    fig = plt.figure(figsize=(w_px / 100, h_px / 100), dpi=100)
    fig.patch.set_facecolor(BG)
    fig.text(0.5, 0.5, text, ha="center", va="center",
              fontsize=fontsize, color=DARK)
    buf = os.path.join(FRAMES_DIR, "tmp_label.png")
    fig.savefig(buf, facecolor=fig.get_facecolor())
    plt.close(fig)
    return Image.open(buf).convert("RGB").resize((w_px, h_px))


header_title = render_mathtext(title, 15, PW, TITLE_H)
header_eq = render_mathtext(eq, 17, PW, EQ_H)
header_detail = render_mathtext(detail, 10.5, PW, DETAIL_H)

if args.channel != "cm":
    cmap_name = "viridis" if args.channel == "TB" else "bwr"
    vmin = 0.0 if args.channel == "TB" else -args.vmax
    fig, ax = plt.subplots(figsize=(PW * 0.55 / 100, 64 / 100), dpi=100)
    fig.subplots_adjust(left=0.08, right=0.92, top=0.55, bottom=0.42)
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=args.vmax)
    cb = matplotlib.colorbar.ColorbarBase(
        ax, cmap=matplotlib.colormaps[cmap_name], norm=norm,
        orientation="horizontal")
    cb.set_label("energía por átomo [ε]", fontsize=8, color=DARK)
    cb.ax.tick_params(labelsize=7.5, labelcolor=DARK)
    for spine in cb.ax.spines.values():
        spine.set_edgecolor(DARK)
    fig.patch.set_facecolor(BG); ax.set_facecolor(BG)
    cbar_path = os.path.join(FRAMES_DIR, "cbar.png")
    fig.savefig(cbar_path, facecolor=fig.get_facecolor())
    plt.close(fig)
    cbar_img = Image.open(cbar_path).convert("RGB").resize((int(PW * 0.55), 64))
else:
    cbar_img = None

# --- Gráfica embebida: alpha_channel(t) vs t, con marcador que avanza -------
PLOT_COLOR = "#2a5f9e"
plot_fig, plot_ax = plt.subplots(figsize=(PW / 100, PLOT_H / 100), dpi=100)
plot_fig.patch.set_facecolor(BG); plot_ax.set_facecolor(BG)
plot_ax.plot(t_b, pct, color=PLOT_COLOR, lw=1.4, zorder=2)
plot_ax.axhline(0, color="#bbbbbb", lw=0.6, zorder=1)
plot_ax.axvspan(args.t0, args.t1, color=PLOT_COLOR, alpha=0.06, zorder=0)
t_lo, t_hi = float(t_b.min()), float(t_b.max())
y_lo, y_hi = float(pct.min()), float(pct.max())
y_pad = 0.08 * (y_hi - y_lo if y_hi > y_lo else 1.0)
y_lo, y_hi = y_lo - y_pad, y_hi + y_pad
plot_ax.set_xlim(t_lo, t_hi)
plot_ax.set_ylim(y_lo, y_hi)
plot_ax.set_xlabel("t  [τ]", fontsize=8.5, color=DARK)
plot_ax.set_ylabel(f"α ({args.channel})  [% de E₀]", fontsize=8.5, color=DARK)
plot_ax.tick_params(labelsize=7.5, colors=DARK)
for spine in plot_ax.spines.values():
    spine.set_edgecolor(DARK)
plot_fig.subplots_adjust(left=0.11, right=0.98, top=0.94, bottom=0.26)
plot_fig.canvas.draw()
_bbox = plot_ax.get_window_extent()
AX_X0, AX_X1 = _bbox.x0, _bbox.x1
AX_Y0, AX_Y1 = _bbox.y0, _bbox.y1   # coords matplotlib (origen abajo-izq)
plot_path = os.path.join(FRAMES_DIR, "plot_base.png")
plot_fig.savefig(plot_path, facecolor=plot_fig.get_facecolor())
plt.close(plot_fig)
plot_base_img = Image.open(plot_path).convert("RGB")
if plot_base_img.size != (PW, PLOT_H):
    plot_base_img = plot_base_img.resize((PW, PLOT_H))


def plot_marker_xy(tau, val):
    fx = (tau - t_lo) / (t_hi - t_lo)
    fy = (val - y_lo) / (y_hi - y_lo)
    px = AX_X0 + fx * (AX_X1 - AX_X0)
    py = PLOT_H - (AX_Y0 + fy * (AX_Y1 - AX_Y0))
    return px, py


AX_TOP_PX = PLOT_H - AX_Y1
AX_BOT_PX = PLOT_H - AX_Y0

try:
    font_pct = ImageFont.truetype(
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 34)
    font_time = ImageFont.truetype(
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 16)
    font_note = ImageFont.truetype(
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 12)
except OSError:
    font_pct = font_time = font_note = ImageFont.load_default()

W = PW
H = HEADER_H + PH + PLOT_H + FOOTER_H
BG255 = tuple(int(c * 255) for c in BG)

print("Componiendo frames finales ...")
for i, fr in enumerate(frames):
    tau = fr * dt_dump + times[0]
    val = pct_at(tau)
    canvas = Image.new("RGB", (W, H), BG255)
    dr = ImageDraw.Draw(canvas)

    canvas.paste(header_title, (0, HEAD_MARGIN))
    canvas.paste(header_eq, (0, HEAD_MARGIN + TITLE_H))
    canvas.paste(header_detail, (0, HEAD_MARGIN + TITLE_H + EQ_H))
    dr.text((W - 14, 8), f"t = {tau:5.1f} τ", font=font_time,
             fill=(30, 30, 40), anchor="ra")
    dr.text((W - 14, 30), f"h₀ = {h0_eff:.1f} σ   (v_imp = {v_imp:.2f} σ/τ)",
             font=font_note, fill=(70, 70, 80), anchor="ra")

    im = Image.open(os.path.join(RAW, f"f_{i:04d}.png")).convert("RGB")
    canvas.paste(im, (0, HEADER_H))

    plot_y0 = HEADER_H + PH
    canvas.paste(plot_base_img, (0, plot_y0))
    mx, my = plot_marker_xy(tau, val)
    dr.line([(mx, plot_y0 + AX_TOP_PX), (mx, plot_y0 + AX_BOT_PX)],
            fill=(200, 60, 60), width=1)
    dr.ellipse([mx - 4, plot_y0 + my - 4, mx + 4, plot_y0 + my + 4],
               fill=(200, 30, 30))

    foot_y = HEADER_H + PH + PLOT_H
    dr.text((16, foot_y + FOOTER_H / 2), f"{val:5.1f} %  de  E₀",
             font=font_pct, fill=(20, 60, 120), anchor="lm")
    if cbar_img is not None:
        canvas.paste(cbar_img, (W - cbar_img.width - 14, foot_y + 20))
    else:
        dr.text((W - 14, foot_y + FOOTER_H / 2), "canal colectivo del CM\n"
                 "(no se reparte entre átomos)", font=font_note,
                 fill=(70, 70, 80), anchor="rm", align="right")

    canvas.save(os.path.join(FRAMES_DIR, f"frame_{i:04d}.png"))
    if i % 50 == 0 or i == len(frames) - 1:
        print(f"  [{100*i//len(frames):3d}%] frame {fr:4d}  {args.channel}={val:.1f}%",
              flush=True)

r = subprocess.run([
    "ffmpeg", "-y", "-framerate", str(args.fps),
    "-i", os.path.join(FRAMES_DIR, "frame_%04d.png"),
    "-c:v", "libx264", "-preset", "slow", "-crf", "17",
    "-pix_fmt", "yuv420p", "-movflags", "+faststart", args.out],
    capture_output=True, text=True)
if r.returncode != 0:
    print("ERROR ffmpeg:", r.stderr[-600:])
    sys.exit(1)
print(f"\nVideo listo ({os.path.getsize(args.out)/1e6:.1f} MB): {args.out}")
