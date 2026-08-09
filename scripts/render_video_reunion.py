"""
render_video_reunion.py — Renderiza un video del impacto con OVITO.
Generaliza render_impact_ovito.py: trayectoria, ventana temporal y escala
de color como argumentos; añade etiqueta de tiempo t/τ.

Uso:
    LD_LIBRARY_PATH=~/ovito_env/lib/python3.12/site-packages/PySide6/Qt/lib:$LD_LIBRARY_PATH \\
    ~/ovito_env/bin/python scripts/render_video_reunion.py \\
        --traj data/outputs/activo/viz_N79_h28/traj.lammpstrj \\
        --out docs/videos_reunion/video_N79_v116_solido.mp4 \\
        --t0 30 --t1 95 --kemax 0.5
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
               help="tope de la escala de color de KE [eps]")
p.add_argument("--dt", type=float, default=0.005, help="paso MD [tau]")
p.add_argument("--fps", type=int, default=25)
p.add_argument("--size", default="1280x720")
args = p.parse_args()

W, H = (int(x) for x in args.size.split("x"))
FRAMES_DIR = os.path.join(
    os.environ.get("TMPDIR", "/tmp"),
    "reunion_frames_" + os.path.splitext(os.path.basename(args.out))[0])
shutil.rmtree(FRAMES_DIR, ignore_errors=True)
os.makedirs(FRAMES_DIR, exist_ok=True)
os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer, TextLabelOverlay
from ovito.modifiers import ColorCodingModifier
from ovito.qt_compat import QtCore

print(f"Cargando {args.traj} ...")
pipeline = import_file(args.traj, sort_particles=True)
nframes = pipeline.source.num_frames
print(f"  {nframes} frames, {pipeline.compute(0).particles.count} átomos")

# tiempo físico de cada frame a partir de su timestep
times = np.array([
    pipeline.compute(f).attributes["Timestep"] * args.dt
    for f in (0, 1)])
dt_dump = times[1] - times[0]
frames = [f for f in range(nframes)
          if args.t0 <= f * dt_dump + times[0] <= args.t1]
print(f"  dt_dump = {dt_dump:.3g} tau -> {len(frames)} frames "
      f"({len(frames)/args.fps:.1f} s de video)")

pipeline.modifiers.append(ColorCodingModifier(
    property="c_keratom", start_value=0.0, end_value=args.kemax,
    gradient=ColorCodingModifier.Hot()))
pipeline.add_to_scene()
pipeline.source.data.particles_.vis.radius = 0.45
pipeline.source.data.cell_.vis.enabled = False

# centrar la cámara cerca del final de la ventana (bola en la superficie)
data_imp = pipeline.compute(frames[len(frames) // 2])
pos = np.array(data_imp.particles["Position"])
xc, yc, _ = pos.mean(axis=0)
zc = 12.0                       # nivel de contacto: encuadre fijo

vp = Viewport()
vp.type = Viewport.Type.Perspective
vp.camera_pos = (xc, yc - 45, zc + 8)
vp.camera_dir = (0.0, 1.0, -0.18)
vp.camera_up = (0.0, 0.0, 1.0)
vp.fov = 0.65

label = TextLabelOverlay(
    text="t = [Timestep] pasos",
    alignment=QtCore.Qt.AlignmentFlag.AlignTop
    | QtCore.Qt.AlignmentFlag.AlignLeft,
    font_size=0.045, text_color=(1.0, 1.0, 1.0), offset_y=-0.05)
vp.overlays.append(label)

renderer = TachyonRenderer()
renderer.ambient_occlusion = True
renderer.shadows = True
renderer.direct_light_intensity = 1.5
renderer.ambient_occlusion_brightness = 1.2

for i, fr in enumerate(frames):
    tau = fr * dt_dump + times[0]
    label.text = f"t = {tau:6.1f} τ"
    vp.render_image(size=(W, H), frame=fr, renderer=renderer,
                    background=(0.05, 0.05, 0.12),
                    filename=os.path.join(FRAMES_DIR, f"frame_{i:04d}.png"))
    if i % 50 == 0 or i == len(frames) - 1:
        print(f"  [{100*i//len(frames):3d}%] frame {fr:4d}  t={tau:6.1f} tau",
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
