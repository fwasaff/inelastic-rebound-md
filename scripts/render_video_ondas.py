"""
render_video_ondas.py — Video de cámara cenital (top-down) de la superficie,
coloreada por KE por átomo, para ver directamente el frente balístico del
impacto expandirse y toparse con su propia imagen periódica (wraparound).
Pensado para la corrida de losa grande (Nx=Ny=40) que resuelve el frente
antes de que se contamine (ver fig_ondas_superficie_grande.png).

La bola se elimina de la escena (no aporta al fenómeno que se quiere ver
aquí); solo se renderiza la superficie.

Uso:
    LD_LIBRARY_PATH=~/ovito_env/lib/python3.12/site-packages/PySide6/Qt/lib \\
    ~/ovito_env/bin/python scripts/render_video_ondas.py \\
        --traj data/outputs/ondas_N79_grande/traj.lammpstrj \\
        --bounce data/outputs/ondas_N79_grande/bounce.dat \\
        --out docs/reunion_gonzalo/videos_reunion/ondas_superficie_grande.mp4 \\
        --vmax 0.003
"""

import argparse
import os
import shutil
import subprocess
import sys

import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--traj", required=True)
p.add_argument("--bounce", required=True)
p.add_argument("--out", required=True)
p.add_argument("--t0", type=float, default=None, help="defecto: t_imp-3")
p.add_argument("--t1", type=float, default=None, help="defecto: t_imp+38")
p.add_argument("--vmax", type=float, default=0.003, help="tope color KE/átomo [eps]")
p.add_argument("--dt", type=float, default=0.005)
p.add_argument("--fps", type=int, default=25)
p.add_argument("--panel", default="800x800")
args = p.parse_args()

PW, PH = (int(x) for x in args.panel.split("x"))
HEADER_H, FOOTER_H = 68, 56
FRAMES_DIR = os.path.join(
    os.environ.get("TMPDIR", "/tmp"),
    "ondas_frames_" + os.path.splitext(os.path.basename(args.out))[0])
RAW = os.path.join(FRAMES_DIR, "raw")
shutil.rmtree(FRAMES_DIR, ignore_errors=True)
os.makedirs(RAW)
os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

# --- t_imp desde bounce.dat ---------------------------------------------------
d = np.loadtxt(args.bounce, comments="#")
tt, Zcm = d[:, 0], d[:, 1]
t_imp = float(tt[np.argmin(Zcm)])
t0 = args.t0 if args.t0 is not None else t_imp - 3.0
t1 = args.t1 if args.t1 is not None else t_imp + 38.0
print(f"t_imp={t_imp:.1f} tau  ventana=[{t0:.1f},{t1:.1f}] tau")

# --- OVITO: solo superficie, coloreada por KE/átomo ---------------------------
from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer
from ovito.modifiers import ColorCodingModifier, ExpressionSelectionModifier, DeleteSelectedModifier

print(f"Cargando {args.traj} ...")
pipe = import_file(args.traj, sort_particles=True)
nframes = pipe.source.num_frames
times = np.array([pipe.compute(f).attributes["Timestep"] * args.dt for f in (0, 1)])
dt_dump = times[1] - times[0]
frames = [f for f in range(nframes) if t0 <= f * dt_dump + times[0] <= t1]
print(f"  {len(frames)} frames ({len(frames)/args.fps:.1f} s de video)")

pipe.modifiers.append(ExpressionSelectionModifier(expression="ParticleType==2"))
pipe.modifiers.append(DeleteSelectedModifier())
pipe.modifiers.append(ColorCodingModifier(
    property="c_keratom", start_value=0.0, end_value=args.vmax,
    gradient=ColorCodingModifier.Viridis()))
pipe.source.data.particles_.vis.radius = 0.55
pipe.source.data.cell_.vis.enabled = False

BG = (0.06, 0.06, 0.09)
data0 = pipe.compute(0)
pos0 = np.array(data0.particles["Position"])
xc, yc = pos0[:, 0].mean(), pos0[:, 1].mean()
Lx = data0.cell[0, 0]

vp = Viewport()
vp.type = Viewport.Type.Ortho
vp.camera_pos = (xc, yc, 40.0)
vp.camera_dir = (0.0, 0.0, -1.0)
vp.camera_up = (0.0, 1.0, 0.0)
vp.fov = Lx * 1.03   # cubre toda la celda con margen chico

renderer = TachyonRenderer()
renderer.ambient_occlusion = False   # vista cenital plana: no aporta, y es cara
renderer.shadows = False
renderer.direct_light_intensity = 1.6

pipe.add_to_scene()
print("Renderizando ...")
for i, fr in enumerate(frames):
    vp.render_image(size=(PW, PH), frame=fr, renderer=renderer,
                     background=BG, filename=os.path.join(RAW, f"f_{i:04d}.png"))
    if i % 40 == 0 or i == len(frames) - 1:
        print(f"  [{100*i//len(frames):3d}%] frame {fr:4d}", flush=True)
pipe.remove_from_scene()

# --- Cabecera + colorbar + composición -----------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

DARK = "#e8e8ee"
fig, ax = plt.subplots(figsize=(PW * 0.5 / 100, 34 / 100), dpi=100)
fig.subplots_adjust(left=0.06, right=0.94, top=0.55, bottom=0.42)
norm = matplotlib.colors.Normalize(vmin=0.0, vmax=args.vmax)
cb = matplotlib.colorbar.ColorbarBase(
    ax, cmap=matplotlib.colormaps["viridis"], norm=norm, orientation="horizontal")
cb.set_label("KE / átomo de superficie  [ε]", fontsize=8, color=DARK)
cb.ax.tick_params(labelsize=7.5, labelcolor=DARK)
for spine in cb.ax.spines.values():
    spine.set_edgecolor(DARK)
fig.patch.set_facecolor(BG); ax.set_facecolor(BG)
cbar_path = os.path.join(FRAMES_DIR, "cbar.png")
fig.savefig(cbar_path, facecolor=fig.get_facecolor())
plt.close(fig)
cbar_img = Image.open(cbar_path).convert("RGB").resize((int(PW * 0.5), 34))

try:
    font_time = ImageFont.truetype(
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 18)
    font_title = ImageFont.truetype(
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 15)
except OSError:
    font_time = font_title = ImageFont.load_default()

W = PW
H = HEADER_H + PH + FOOTER_H
BG255 = tuple(int(c * 255) for c in BG)

print("Componiendo frames finales ...")
for i, fr in enumerate(frames):
    tau = fr * dt_dump + times[0]
    canvas = Image.new("RGB", (W, H), BG255)
    dr = ImageDraw.Draw(canvas)
    dr.text((14, 12), "Superficie (vista cenital)", font=font_title,
             fill=(230, 230, 236))
    dr.text((14, 34), "frente balístico del impacto", font=font_time,
             fill=(180, 180, 190))
    dr.text((W - 14, 20), f"t = {tau:5.1f} τ  (t_imp = {t_imp:.1f} τ)",
             font=font_time, fill=(230, 230, 236), anchor="ra")

    im = Image.open(os.path.join(RAW, f"f_{i:04d}.png")).convert("RGB")
    canvas.paste(im, (0, HEADER_H))

    foot_y = HEADER_H + PH
    canvas.paste(cbar_img, (W // 2 - cbar_img.width // 2, foot_y + 14))

    canvas.save(os.path.join(FRAMES_DIR, f"frame_{i:04d}.png"))
    if i % 40 == 0 or i == len(frames) - 1:
        print(f"  [{100*i//len(frames):3d}%] frame {fr:4d}  t={tau:.1f}", flush=True)

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
