"""
render_impact_ovito.py  —  Renderiza video del impacto usando OVITO
Lee traj.lammpstrj directamente (formato LAMMPS dump custom).

Uso (desde el directorio raíz del proyecto):
    LD_LIBRARY_PATH=~/ovito_env/lib/python3.12/site-packages/PySide6/Qt/lib:$LD_LIBRARY_PATH \\
    ~/ovito_env/bin/python scripts/render_impact_ovito.py

Salida:
    data/outputs/reprod_2012/h0_41.47/impact_video.mp4
"""

import os
import sys
import subprocess
import numpy as np

TRAJ       = "data/outputs/reprod_2012/h0_41.47/traj.lammpstrj"
VIDEO      = "data/outputs/reprod_2012/h0_41.47/impact_video.mp4"
FRAMES_DIR = "/tmp/impact_frames_ovito"

# Rango de frames: t ≈ 25τ (frame 100) a t ≈ 125τ (frame 500)
# El primer impacto ocurre ~t=56τ (frame 224), apogeo ~t=110τ (frame 440)
FRAME_START = 100
FRAME_END   = 500
FPS         = 25
WIDTH       = 1280
HEIGHT      = 720
DT_DUMP     = 0.25   # tau por frame (dump_freq=50, dt=0.005)

os.makedirs(FRAMES_DIR, exist_ok=True)

from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer
from ovito.modifiers import ColorCodingModifier

# --- Cargar trayectoria ----------------------------------------------------
print(f"Cargando {TRAJ} ...")
pipeline = import_file(TRAJ, sort_particles=True)
nframes = pipeline.source.num_frames
natoms  = pipeline.compute(0).particles.count
print(f"  {nframes} frames, {natoms} átomos")

# --- Colorear por KE atómica (escala Hot: negro→rojo→naranja→amarillo→blanco) -----
color_mod = ColorCodingModifier(
    property='c_keratom',
    start_value=0.0,
    end_value=0.5,
    gradient=ColorCodingModifier.Hot(),
)
pipeline.modifiers.append(color_mod)
pipeline.add_to_scene()

# Radio de partículas = 0.45σ, ocultar el marco de la caja de simulación
pipeline.source.data.particles_.vis.radius = 0.45
pipeline.source.data.cell_.vis.enabled = False

# --- Calcular centro de escena en el frame de impacto ----------------------
data_impact = pipeline.compute(224)
pos = np.array(data_impact.particles['Position'])
xc, yc, zc = pos.mean(axis=0)
zmax = pos[:, 2].max()

# --- Viewport: vista frontal lateral, cámara fija --------------------------
# Miramos desde -Y, ligeramente elevada, centrada en el impacto
dist = 45
vp = Viewport()
vp.type = Viewport.Type.Perspective
vp.camera_pos = (xc, yc - dist, zc + 8)
vp.camera_dir = (0.0, 1.0, -0.18)
vp.camera_up  = (0.0, 0.0, 1.0)
vp.fov = 0.65   # ~37° → campo angosto para escenas pequeñas LJ

# --- Renderizador Tachyon --------------------------------------------------
renderer = TachyonRenderer()
renderer.ambient_occlusion = True
renderer.shadows = True
renderer.direct_light_intensity = 1.5
renderer.ambient_occlusion_brightness = 1.2

# --- Renderizar frame a frame ----------------------------------------------
n_total = FRAME_END - FRAME_START + 1
print(f"Renderizando frames {FRAME_START}–{FRAME_END} ({n_total} frames @ {FPS} fps) ...")
print(f"  Duración del video: {n_total/FPS:.1f} s")

for i, fr in enumerate(range(FRAME_START, FRAME_END + 1)):
    if fr >= nframes:
        break

    out_png = os.path.join(FRAMES_DIR, f"frame_{i:04d}.png")
    vp.render_image(
        size=(WIDTH, HEIGHT),
        frame=fr,
        renderer=renderer,
        background=(0.05, 0.05, 0.12),
        filename=out_png,
    )

    if i % 50 == 0 or i == n_total - 1:
        tau = fr * DT_DUMP
        pct = 100 * i // n_total
        print(f"  [{pct:3d}%]  frame {fr:4d}  t = {tau:6.1f} τ")

print(f"\nFrames listos en {FRAMES_DIR}/")

# --- Ensamblar video con ffmpeg --------------------------------------------
print("Generando video con ffmpeg...")
cmd = [
    "ffmpeg", "-y",
    "-framerate", str(FPS),
    "-i", os.path.join(FRAMES_DIR, "frame_%04d.png"),
    "-c:v", "libx264", "-preset", "slow", "-crf", "17",
    "-pix_fmt", "yuv420p", "-movflags", "+faststart",
    VIDEO,
]
r = subprocess.run(cmd, capture_output=True, text=True)
if r.returncode != 0:
    print("ERROR ffmpeg:", r.stderr[-600:])
    sys.exit(1)

mb = os.path.getsize(VIDEO) / 1e6
print(f"\nVideo listo ({mb:.1f} MB): {VIDEO}")
