"""
render_video_canales.py — Renderiza un video de dos paneles sincronizados
que separan visualmente los canales de energía interna: KE por átomo
(alpha_TB + fonones de superficie) y PE relativa al equilibrio por átomo
(alpha_DeltaVB, compresión/expansión). Fondo claro, paletas perceptualmente
uniformes (Viridis para KE, azul-blanco-rojo divergente para PE).

Reemplaza render_video_reunion.py (fondo negro + paleta Hot, poco legible).

Uso:
    LD_LIBRARY_PATH=~/ovito_env/lib/python3.12/site-packages/PySide6/Qt/lib:$LD_LIBRARY_PATH \\
    ~/ovito_env/bin/python scripts/render_video_canales.py \\
        --traj data/outputs/viz_N79_h28/traj.lammpstrj \\
        --out docs/videos_reunion/video_N79_v116_solido.mp4 \\
        --t0 30 --t1 95 --kemax 0.5 --pemax 0.5
"""

import argparse
import os
import shutil
import subprocess
import sys

import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--traj", required=True)
p.add_argument("--out", required=True)
p.add_argument("--t0", type=float, required=True, help="inicio [tau]")
p.add_argument("--t1", type=float, required=True, help="fin [tau]")
p.add_argument("--kemax", type=float, default=0.5,
               help="tope de la escala de color de KE por átomo [eps]")
p.add_argument("--pemax", type=float, default=0.5,
               help="rango simétrico de la escala de color de ΔPE [eps]")
p.add_argument("--dt", type=float, default=0.005, help="paso MD [tau]")
p.add_argument("--fps", type=int, default=25)
p.add_argument("--panel", default="760x600", help="tamaño de cada panel")
args = p.parse_args()

PW, PH = (int(x) for x in args.panel.split("x"))
TOP, TITLE_H, CBAR_H, DIV = 44, 34, 64, 6
FRAMES_DIR = os.path.join(
    os.environ.get("TMPDIR", "/tmp"),
    "canales_frames_" + os.path.splitext(os.path.basename(args.out))[0])
RAW_A = os.path.join(FRAMES_DIR, "a")
RAW_B = os.path.join(FRAMES_DIR, "b")
shutil.rmtree(FRAMES_DIR, ignore_errors=True)
os.makedirs(RAW_A); os.makedirs(RAW_B)
os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer
from ovito.modifiers import (ColorCodingModifier, FreezePropertyModifier,
                              ComputePropertyModifier)

BG = (0.97, 0.97, 0.985)          # fondo claro (no negro)
TXT_DARK = (0.12, 0.12, 0.16)

print(f"Cargando {args.traj} ...")
pipe_ref = import_file(args.traj, sort_particles=True)
nframes = pipe_ref.source.num_frames
print(f"  {nframes} frames, {pipe_ref.compute(0).particles.count} átomos")

times = np.array([pipe_ref.compute(f).attributes["Timestep"] * args.dt
                   for f in (0, 1)])
dt_dump = times[1] - times[0]
frames = [f for f in range(nframes)
          if args.t0 <= f * dt_dump + times[0] <= args.t1]
print(f"  dt_dump = {dt_dump:.3g} tau -> {len(frames)} frames "
      f"({len(frames)/args.fps:.1f} s de video)")

# --- Pipeline A: KE por átomo (Viridis) ------------------------------------
pipe_ke = import_file(args.traj, sort_particles=True)
pipe_ke.modifiers.append(ColorCodingModifier(
    property="c_keratom", start_value=0.0, end_value=args.kemax,
    gradient=ColorCodingModifier.Viridis()))
pipe_ke.source.data.particles_.vis.radius = 0.45
pipe_ke.source.data.cell_.vis.enabled = False

# --- Pipeline B: PE relativa al equilibrio por átomo (azul-blanco-rojo) ---
pipe_pe = import_file(args.traj, sort_particles=True)
pipe_pe.modifiers.append(FreezePropertyModifier(
    source_property="c_pe_atom", destination_property="PE0", freeze_at=0))
pipe_pe.modifiers.append(ComputePropertyModifier(
    output_property="dPE", expressions=["c_pe_atom - PE0"]))
pipe_pe.modifiers.append(ColorCodingModifier(
    property="dPE", start_value=-args.pemax, end_value=args.pemax,
    gradient=ColorCodingModifier.BlueWhiteRed()))
pipe_pe.source.data.particles_.vis.radius = 0.45
pipe_pe.source.data.cell_.vis.enabled = False

# --- Cámara común (misma para ambos paneles) -------------------------------
data_imp = pipe_ke.compute(frames[len(frames) // 2])
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

print("Renderizando panel A (KE) ...")
pipe_ke.add_to_scene()
for i, fr in enumerate(frames):
    vp.render_image(size=(PW, PH), frame=fr, renderer=renderer,
                     background=BG,
                     filename=os.path.join(RAW_A, f"f_{i:04d}.png"))
    if i % 50 == 0 or i == len(frames) - 1:
        print(f"  [{100*i//len(frames):3d}%] frame {fr:4d}", flush=True)
pipe_ke.remove_from_scene()

print("Renderizando panel B (dPE) ...")
pipe_pe.add_to_scene()
for i, fr in enumerate(frames):
    vp.render_image(size=(PW, PH), frame=fr, renderer=renderer,
                     background=BG,
                     filename=os.path.join(RAW_B, f"f_{i:04d}.png"))
    if i % 50 == 0 or i == len(frames) - 1:
        print(f"  [{100*i//len(frames):3d}%] frame {fr:4d}", flush=True)
pipe_pe.remove_from_scene()

# --- Colorbars estáticas (una vez) -----------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

def make_colorbar(cmap_name, vmin, vmax, label, width_px, height_px):
    fig, ax = plt.subplots(figsize=(width_px / 100, height_px / 100), dpi=100)
    fig.subplots_adjust(left=0.06, right=0.94, top=0.55, bottom=0.42)
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    cb = matplotlib.colorbar.ColorbarBase(
        ax, cmap=matplotlib.colormaps[cmap_name], norm=norm,
        orientation="horizontal")
    cb.set_label(label, fontsize=9, color="#1e1e28")
    cb.ax.tick_params(labelsize=8, color="#1e1e28", labelcolor="#1e1e28")
    for spine in cb.ax.spines.values():
        spine.set_edgecolor("#1e1e28")
    fig.patch.set_facecolor((0.97, 0.97, 0.985))
    ax.set_facecolor((0.97, 0.97, 0.985))
    buf = os.path.join(FRAMES_DIR, f"cbar_{cmap_name}.png")
    fig.savefig(buf, facecolor=fig.get_facecolor())
    plt.close(fig)
    return Image.open(buf).convert("RGB").resize((width_px, height_px))

cbar_ke = make_colorbar("viridis", 0.0, args.kemax,
                         "KE por átomo [ε]  (α_TB + fonones)", PW, CBAR_H)
cbar_pe = make_colorbar("bwr", -args.pemax, args.pemax,
                         "ΔPE por átomo [ε]  (α_ΔVB: azul=enlace, rojo=expansión)",
                         PW, CBAR_H)

# --- Composición final: título superior + 2 paneles + colorbars ------------
W = PW * 2 + DIV
H = TOP + TITLE_H + PH + CBAR_H
try:
    font_time = ImageFont.truetype(
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 20)
    font_title = ImageFont.truetype(
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 16)
except OSError:
    font_time = font_title = ImageFont.load_default()

print("Componiendo frames finales ...")
for i, fr in enumerate(frames):
    tau = fr * dt_dump + times[0]
    canvas = Image.new("RGB", (W, H), (int(BG[0]*255), int(BG[1]*255), int(BG[2]*255)))
    d = ImageDraw.Draw(canvas)

    d.text((W / 2, TOP / 2), f"t = {tau:6.1f} τ", font=font_time,
           fill=(30, 30, 40), anchor="mm")
    d.text((PW / 2, TOP + TITLE_H / 2), "Energía cinética por átomo",
           font=font_title, fill=(30, 30, 40), anchor="mm")
    d.text((PW + DIV + PW / 2, TOP + TITLE_H / 2), "Energía potencial relativa al equilibrio",
           font=font_title, fill=(30, 30, 40), anchor="mm")

    imA = Image.open(os.path.join(RAW_A, f"f_{i:04d}.png")).convert("RGB")
    imB = Image.open(os.path.join(RAW_B, f"f_{i:04d}.png")).convert("RGB")
    canvas.paste(imA, (0, TOP + TITLE_H))
    canvas.paste(imB, (PW + DIV, TOP + TITLE_H))
    canvas.paste(cbar_ke, (0, TOP + TITLE_H + PH))
    canvas.paste(cbar_pe, (PW + DIV, TOP + TITLE_H + PH))

    d.rectangle([PW, TOP + TITLE_H, PW + DIV, TOP + TITLE_H + PH],
                fill=(210, 210, 220))

    canvas.save(os.path.join(FRAMES_DIR, f"frame_{i:04d}.png"))
    if i % 50 == 0 or i == len(frames) - 1:
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
