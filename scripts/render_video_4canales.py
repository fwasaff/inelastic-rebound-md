"""
render_video_4canales.py — Combina los 4 canales de energía en un solo video,
grilla 2x2 sincronizada: alpha_cm | alpha_TB / alpha_VB | alpha_surf. Cada
panel lleva su ecuación y su % de E0 en vivo (misma lógica de coloreado que
render_video_canal.py, unificada en un solo render).

Uso:
    LD_LIBRARY_PATH=~/ovito_env/lib/python3.12/site-packages/PySide6/Qt/lib:$LD_LIBRARY_PATH \\
    ~/ovito_env/bin/python scripts/render_video_4canales.py \\
        --traj data/outputs/viz_N79_h28/traj.lammpstrj \\
        --bounce data/outputs/viz_N79_h28/bounce.dat \\
        --out docs/reunion_gonzalo/videos_reunion/video_4canales_solido.mp4 \\
        --N 79 --t0 30 --t1 95
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
p.add_argument("--bounce", required=True)
p.add_argument("--out", required=True)
p.add_argument("--N", type=int, required=True)
p.add_argument("--t0", type=float, required=True)
p.add_argument("--t1", type=float, required=True)
p.add_argument("--vmax-tb", type=float, default=0.9)
p.add_argument("--vmax-vb", type=float, default=0.7)
p.add_argument("--vmax-surf", type=float, default=0.04)
p.add_argument("--dt", type=float, default=0.005)
p.add_argument("--fps", type=int, default=25)
p.add_argument("--cell", default="620x420", help="tamaño del panel 3D en cada celda")
args = p.parse_args()

G = 0.026
CW, CH = (int(x) for x in args.cell.split("x"))
HEAD_H, FOOT_H, DIV = 70, 56, 6
CHANNELS = ["cm", "TB", "VB", "surf"]   # orden de lectura: sup-izq, sup-der, inf-izq, inf-der

FRAMES_DIR = os.path.join(
    os.environ.get("TMPDIR", "/tmp"),
    "canal4_frames_" + os.path.splitext(os.path.basename(args.out))[0])
RAW = {c: os.path.join(FRAMES_DIR, c) for c in CHANNELS}
shutil.rmtree(FRAMES_DIR, ignore_errors=True)
for d_ in RAW.values():
    os.makedirs(d_)
os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

# --- Series temporales de cada canal (bounce.dat) ---------------------------
d = np.loadtxt(args.bounce, comments="#")
t_b, Zcm, KB, TB, VB, Uc, UF, KB3d, KS, PEs = (
    d[:, 0], d[:, 1], d[:, 3], d[:, 4], d[:, 5], d[:, 6], d[:, 7], d[:, 8],
    d[:, 9], d[:, 10])
mins, _ = find_peaks(-Zcm, prominence=0.5, distance=50)
z0 = float(Zcm[mins].mean()) if len(mins) else float(Zcm.min())
h0_eff = Zcm[0] - z0
E0 = args.N * G * h0_eff

PCT = {
    "cm":   (KB3d + args.N * G * (Zcm - z0)) / E0 * 100.0,
    "TB":   TB / E0 * 100.0,
    "VB":   (VB - VB[0]) / E0 * 100.0,
    "surf": (KS + PEs - KS[0] - PEs[0]) / E0 * 100.0,
}
EQ = {
    "cm":   r"$\alpha_{\rm cm}=\dfrac{h_1}{h_0}$  (en el apogeo)",
    "TB":   r"$\alpha_{T_B}(t)=\dfrac{T_B(t)}{E_0}$",
    "VB":   r"$\alpha_{\Delta V_B}(t)=\dfrac{V_B(t)-V_B(0)}{E_0}$",
    "surf": r"$\alpha_{\rm surf}(t)=\dfrac{\Delta(K_S+V_{\rm surf})}{E_0}$",
}
TITLE = {
    "cm":   "CM: energía mecánica",
    "TB":   "Bola: calor cinético",
    "VB":   "Bola: calor potencial",
    "surf": "Superficie: fonones",
}


def pct_at(channel, tau):
    return float(np.interp(tau, t_b, PCT[channel]))


# --- OVITO: una pipeline por canal, coloreado específico --------------------
from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer
from ovito.modifiers import ColorCodingModifier

print(f"Cargando {args.traj} ...")
pipe_ref = import_file(args.traj, sort_particles=True)
nframes = pipe_ref.source.num_frames
times = np.array([pipe_ref.compute(f).attributes["Timestep"] * args.dt for f in (0, 1)])
dt_dump = times[1] - times[0]
frames = [f for f in range(nframes)
          if args.t0 <= f * dt_dump + times[0] <= args.t1]
print(f"  {len(frames)} frames ({len(frames)/args.fps:.1f} s)")

d0 = pipe_ref.compute(0)
id0 = np.asarray(d0.particles["Particle Identifier"])
pe0_map = dict(zip(id0.tolist(), np.asarray(d0.particles["c_pe_atom"]).tolist()))

GRAY = np.array([0.80, 0.80, 0.83])
HIGHLIGHT = np.array([0.20, 0.45, 0.78])
VMAX = {"TB": args.vmax_tb, "VB": args.vmax_vb, "surf": args.vmax_surf}


def make_channel_fn(channel):
    def compute_channel(frame, data):
        typ = np.asarray(data.particles["Particle Type"])
        ids = np.asarray(data.particles["Particle Identifier"])
        is_ball = typ == 2
        is_surf = typ == 1
        if channel == "TB":
            vel = np.asarray(data.particles["Velocity"])
            v_cm = vel[is_ball].mean(axis=0) if is_ball.any() else np.zeros(3)
            v_int = vel.copy(); v_int[is_ball] -= v_cm
            chan = 0.5 * np.sum(v_int**2, axis=1)
            active = is_ball
        elif channel == "VB":
            pe = np.asarray(data.particles["c_pe_atom"])
            pe0 = np.array([pe0_map[i] for i in ids])
            chan = pe - pe0
            active = is_ball
        elif channel == "surf":
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
        typ = np.asarray(data.particles["Particle Type"])
        if "Color" in data.particles:
            color = np.asarray(data.particles["Color"]).copy()
        else:
            color = np.tile(GRAY, (len(typ), 1))
        if channel == "cm":
            color[typ == 2] = HIGHLIGHT
            color[typ == 1] = GRAY
        else:
            active = np.asarray(data.particles["Active"])
            color[active < 0.5] = GRAY
        data.particles_.create_property("Color", data=color)

    return compute_channel, flatten_inactive


pipes = {}
for ch in CHANNELS:
    pp = import_file(args.traj, sort_particles=True)
    comp_fn, flat_fn = make_channel_fn(ch)
    pp.modifiers.append(comp_fn)
    if ch != "cm":
        gradient = ColorCodingModifier.Viridis() if ch == "TB" else ColorCodingModifier.BlueWhiteRed()
        start = 0.0 if ch == "TB" else -VMAX[ch]
        pp.modifiers.append(ColorCodingModifier(
            property="Chan", start_value=start, end_value=VMAX[ch], gradient=gradient))
    pp.modifiers.append(flat_fn)
    pp.source.data.particles_.vis.radius = 0.45
    pp.source.data.cell_.vis.enabled = False
    pipes[ch] = pp

BG = (0.97, 0.97, 0.985)
data_imp = pipes["cm"].compute(frames[len(frames) // 2])
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

for ch in CHANNELS:
    print(f"Renderizando panel {ch} ...")
    pipes[ch].add_to_scene()
    for i, fr in enumerate(frames):
        vp.render_image(size=(CW, CH), frame=fr, renderer=renderer,
                         background=BG, filename=os.path.join(RAW[ch], f"f_{i:04d}.png"))
        if i % 80 == 0 or i == len(frames) - 1:
            print(f"  [{100*i//len(frames):3d}%] frame {fr:4d}", flush=True)
    pipes[ch].remove_from_scene()

# --- Cabeceras (título+ecuación) y colorbars, una vez por canal ------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

DARK = "#1e1e28"


def render_mathtext(text, fontsize, w_px, h_px, tag="lbl"):
    SS = 3  # supersample para texto nítido al reducir a tamaño final
    fig = plt.figure(figsize=(w_px * SS / 100, h_px * SS / 100), dpi=100)
    fig.patch.set_facecolor(BG)
    fig.text(0.5, 0.5, text, ha="center", va="center",
              fontsize=fontsize * SS, color=DARK)
    buf = os.path.join(FRAMES_DIR, f"tmp_{tag}.png")
    fig.savefig(buf, facecolor=fig.get_facecolor())
    plt.close(fig)
    img = Image.open(buf).convert("RGB")
    return img.resize((w_px, h_px), Image.LANCZOS)


headers = {}
for ch in CHANNELS:
    htitle = render_mathtext(TITLE[ch], 13, CW, 24, tag=f"title_{ch}")
    heq = render_mathtext(EQ[ch], 14, CW, 46, tag=f"eq_{ch}")
    headers[ch] = (htitle, heq)

cbars = {}
for ch in ["TB", "VB", "surf"]:
    cmap_name = "viridis" if ch == "TB" else "bwr"
    vmin = 0.0 if ch == "TB" else -VMAX[ch]
    fig, ax = plt.subplots(figsize=(CW * 0.5 / 100, 26 / 100), dpi=100)
    fig.subplots_adjust(left=0.06, right=0.94, top=0.62, bottom=0.38)
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=VMAX[ch])
    cb = matplotlib.colorbar.ColorbarBase(
        ax, cmap=matplotlib.colormaps[cmap_name], norm=norm, orientation="horizontal")
    cb.ax.tick_params(labelsize=6.5, labelcolor=DARK)
    for spine in cb.ax.spines.values():
        spine.set_edgecolor(DARK)
    fig.patch.set_facecolor(BG); ax.set_facecolor(BG)
    cp = os.path.join(FRAMES_DIR, f"cbar_{ch}.png")
    fig.savefig(cp, facecolor=fig.get_facecolor()); plt.close(fig)
    cbars[ch] = Image.open(cp).convert("RGB").resize((int(CW * 0.5), 26))

try:
    font_pct = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 22)
    font_time = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 18)
    font_note = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 10)
except OSError:
    font_pct = font_time = font_note = ImageFont.load_default()

CELL_W = CW
CELL_H = HEAD_H + CH + FOOT_H
TOP = 40
W = CELL_W * 2 + DIV
H = TOP + CELL_H * 2 + DIV
BG255 = tuple(int(c * 255) for c in BG)
POS = {"cm": (0, 0), "TB": (1, 0), "VB": (0, 1), "surf": (1, 1)}

print("Componiendo frames finales ...")
for i, fr in enumerate(frames):
    tau = fr * dt_dump + times[0]
    canvas = Image.new("RGB", (W, H), BG255)
    dr = ImageDraw.Draw(canvas)
    dr.text((W / 2, TOP / 2), f"t = {tau:6.1f} τ", font=font_time,
             fill=(30, 30, 40), anchor="mm")

    for ch, (cx, cy) in POS.items():
        x0 = cx * (CELL_W + DIV)
        y0 = TOP + cy * (CELL_H + DIV)
        htitle, heq = headers[ch]
        canvas.paste(htitle, (x0, y0))
        canvas.paste(heq, (x0, y0 + 24))
        im = Image.open(os.path.join(RAW[ch], f"f_{i:04d}.png")).convert("RGB")
        canvas.paste(im, (x0, y0 + HEAD_H))
        val = pct_at(ch, tau)
        foot_y = y0 + HEAD_H + CH
        dr.text((x0 + 10, foot_y + FOOT_H / 2), f"{val:5.1f} %",
                 font=font_pct, fill=(20, 60, 120), anchor="lm")
        if ch in cbars:
            canvas.paste(cbars[ch], (x0 + CELL_W - cbars[ch].width - 10, foot_y + 14))
        else:
            dr.text((x0 + CELL_W - 10, foot_y + FOOT_H / 2),
                     "canal colectivo (no por átomo)", font=font_note,
                     fill=(70, 70, 80), anchor="rm")

    if DIV > 0:
        dr.rectangle([CELL_W, TOP, CELL_W + DIV, H], fill=(210, 210, 220))
        dr.rectangle([0, TOP + CELL_H, W, TOP + CELL_H + DIV], fill=(210, 210, 220))

    canvas.save(os.path.join(FRAMES_DIR, f"frame_{i:04d}.png"))
    if i % 80 == 0 or i == len(frames) - 1:
        print(f"  [{100*i//len(frames):3d}%] frame {fr:4d}", flush=True)

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
