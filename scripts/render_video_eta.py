"""
render_video_eta.py — Video top-down coloreado por la anisotropía de
velocidad local eta(r,t) (Sec. 6 de onda_vs_calor_fonones.tex): rojo =
movimiento todavía coherente/radial (sonido, exergía), azul = velocidades
isotrópicas (calor, anergía). Mismo cálculo que scripts/medir_anisotropia_eta.py
pero resuelto en 2D (grilla espacial x,y, no solo bins radiales) para que se
vea la estructura real, no solo el promedio azimutal.

eta se calcula por celda de una grilla (no por átomo: un solo átomo no
tiene "anisotropía", hace falta promediar sobre varios) y se pinta cada
átomo del color de su celda. Celdas con pocos átomos (<min_atoms) se dejan
sin colorear (gris).

Uso:
    LD_LIBRARY_PATH=~/ovito_env/lib/python3.12/site-packages/PySide6/Qt/lib \\
    ~/ovito_env/bin/python scripts/render_video_eta.py \\
        --traj data/outputs/ondas_N79_libre/traj.lammpstrj \\
        --bounce data/outputs/ondas_N79_libre/bounce.dat \\
        --out docs/reunion_gonzalo/videos_reunion/eta_superficie.mp4 \\
        --cell 4.0
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
p.add_argument("--t0", type=float, default=None, help="defecto: t_imp-1")
p.add_argument("--t1", type=float, default=None, help="defecto: t_imp+15")
p.add_argument("--cell", type=float, default=4.0, help="tamaño de celda de la grilla [sigma]")
p.add_argument("--min-atoms", type=int, default=8, help="átomos mínimos por celda para calcular eta")
p.add_argument("--vmin", type=float, default=-0.5)
p.add_argument("--vmax", type=float, default=1.0)
p.add_argument("--dt", type=float, default=0.005)
p.add_argument("--fps", type=int, default=15)
p.add_argument("--panel", default="800x800")
args = p.parse_args()

PW, PH = (int(x) for x in args.panel.split("x"))
HEADER_H, FOOTER_H = 68, 56
FRAMES_DIR = os.path.join(
    os.environ.get("TMPDIR", "/tmp"),
    "eta_frames_" + os.path.splitext(os.path.basename(args.out))[0])
RAW = os.path.join(FRAMES_DIR, "raw")
shutil.rmtree(FRAMES_DIR, ignore_errors=True)
os.makedirs(RAW)
os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

d = np.loadtxt(args.bounce, comments="#")
tt, Zcm = d[:, 0], d[:, 1]
t_imp = float(tt[np.argmin(Zcm)])
t0 = args.t0 if args.t0 is not None else t_imp - 1.0
t1 = args.t1 if args.t1 is not None else t_imp + 15.0
print(f"t_imp={t_imp:.1f} tau  ventana=[{t0:.1f},{t1:.1f}] tau")

from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer
from ovito.modifiers import ExpressionSelectionModifier, DeleteSelectedModifier

print(f"Cargando {args.traj} ...")
pipe = import_file(args.traj, sort_particles=True)
nframes = pipe.source.num_frames
times = np.array([pipe.compute(f).attributes["Timestep"] * args.dt for f in (0, 1)])
dt_dump = times[1] - times[0]
frames = [f for f in range(nframes) if t0 <= f * dt_dump + times[0] <= t1]
print(f"  {len(frames)} frames ({len(frames)/args.fps:.1f} s de video)")

data0 = pipe.compute(0)
typ0 = np.asarray(data0.particles["Particle Type"])
pos0 = np.asarray(data0.particles["Position"])
is_ball0 = typ0 == 2
cx, cy = pos0[is_ball0, 0].mean(), pos0[is_ball0, 1].mean()
is_surf0 = typ0 == 1
cz = pos0[is_surf0, 2].max()
print(f"  punto de impacto (x,y,z) = ({cx:.2f}, {cy:.2f}, {cz:.2f}) sigma")

CELL = args.cell
MIN_ATOMS = args.min_atoms


def compute_eta(frame, data):
    typ = np.asarray(data.particles["Particle Type"])
    pos = np.asarray(data.particles["Position"])
    vel = np.asarray(data.particles["Velocity"])
    is_surf = typ == 1

    eta_atom = np.full(len(typ), np.nan)
    if is_surf.any():
        x, y, z = pos[is_surf, 0], pos[is_surf, 1], pos[is_surf, 2]
        vx, vy, vz = vel[is_surf, 0], vel[is_surf, 1], vel[is_surf, 2]
        dx, dy, dz = x - cx, y - cy, z - cz
        r3 = np.sqrt(dx**2 + dy**2 + dz**2)
        r3 = np.where(r3 < 1e-6, 1e-6, r3)
        v_r = (vx * dx + vy * dy + vz * dz) / r3
        v2_perp = ((vx**2 + vy**2 + vz**2) - v_r**2) / 2.0

        ix = np.floor(dx / CELL).astype(np.int64)
        iy = np.floor(dy / CELL).astype(np.int64)
        key = ix * 1_000_000 + iy   # hash de celda 2D
        uniq, inv, counts = np.unique(key, return_inverse=True, return_counts=True)

        sum_vr2 = np.bincount(inv, weights=v_r**2, minlength=len(uniq))
        sum_vp2 = np.bincount(inv, weights=v2_perp, minlength=len(uniq))
        mean_vr2 = sum_vr2 / counts
        mean_vp2 = sum_vp2 / counts
        with np.errstate(invalid="ignore", divide="ignore"):
            eta_cell = (mean_vr2 - mean_vp2) / (mean_vr2 + mean_vp2)
        eta_cell[counts < MIN_ATOMS] = np.nan

        eta_surf = eta_cell[inv]
        eta_atom[is_surf] = eta_surf

    data.particles_.create_property("Eta", data=eta_atom)
    data.particles_.create_property("HasEta", data=(~np.isnan(eta_atom)).astype(float))


GRAY = np.array([0.80, 0.80, 0.83])


def color_by_eta(frame, data):
    eta = np.asarray(data.particles["Eta"])
    has = np.asarray(data.particles["HasEta"]) > 0.5
    norm = np.clip((eta - args.vmin) / (args.vmax - args.vmin), 0.0, 1.0)
    # RdBu_r manual: azul (bajo) - blanco (medio) - rojo (alto)
    color = np.zeros((len(eta), 3))
    lo = np.array([0.13, 0.40, 0.72])   # azul
    mid = np.array([0.98, 0.98, 0.96])  # casi blanco
    hi = np.array([0.70, 0.06, 0.09])   # rojo oscuro
    below = norm < 0.5
    above = ~below
    s = np.where(below, norm / 0.5, (norm - 0.5) / 0.5)
    color[below] = (1 - s[below, None]) * lo + s[below, None] * mid
    color[above] = (1 - s[above, None]) * mid + s[above, None] * hi
    color[~has] = GRAY
    data.particles_.create_property("Color", data=color)


typ_del = ExpressionSelectionModifier(expression="ParticleType==2")
pipe.modifiers.append(typ_del)
pipe.modifiers.append(DeleteSelectedModifier())
pipe.modifiers.append(compute_eta)
pipe.modifiers.append(color_by_eta)
pipe.source.data.particles_.vis.radius = 0.55
pipe.source.data.cell_.vis.enabled = False

BG = (0.06, 0.06, 0.09)
Lx = pipe.compute(0).cell[0, 0]

vp = Viewport()
vp.type = Viewport.Type.Ortho
vp.camera_pos = (cx, cy, 40.0)
vp.camera_dir = (0.0, 0.0, -1.0)
vp.camera_up = (0.0, 1.0, 0.0)
vp.fov = Lx * 1.03

renderer = TachyonRenderer()
renderer.ambient_occlusion = False
renderer.shadows = False
renderer.direct_light_intensity = 1.6

pipe.add_to_scene()
print("Renderizando ...")
for i, fr in enumerate(frames):
    vp.render_image(size=(PW, PH), frame=fr, renderer=renderer,
                     background=BG, filename=os.path.join(RAW, f"f_{i:04d}.png"))
    if i % 20 == 0 or i == len(frames) - 1:
        print(f"  [{100*i//len(frames):3d}%] frame {fr:4d}", flush=True)
pipe.remove_from_scene()

# --- Cabecera + colorbar + composición -----------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

DARK = "#e8e8ee"
_cmap_nodes = [(0.0, (0.13, 0.40, 0.72)), (0.5, (0.98, 0.98, 0.96)), (1.0, (0.70, 0.06, 0.09))]
_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("eta_rdbu", _cmap_nodes)

fig, ax = plt.subplots(figsize=(PW * 0.5 / 100, 34 / 100), dpi=100)
fig.subplots_adjust(left=0.06, right=0.94, top=0.55, bottom=0.42)
norm = matplotlib.colors.Normalize(vmin=args.vmin, vmax=args.vmax)
cb = matplotlib.colorbar.ColorbarBase(ax, cmap=_cmap, norm=norm, orientation="horizontal")
cb.set_label(r"$\eta$  (1=sonido coherente, 0=isotrópico, calor)", fontsize=8, color=DARK)
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
    dr.text((14, 12), "Anisotropía de velocidad η (vista cenital)",
             font=font_title, fill=(230, 230, 236))
    dr.text((14, 34), "rojo = sonido (exergía)   gris = azul = calor (anergía)",
             font=font_time, fill=(180, 180, 190))
    dr.text((W - 14, 20), f"t = {tau:5.1f} τ  (t_imp = {t_imp:.1f} τ)",
             font=font_time, fill=(230, 230, 236), anchor="ra")

    im = Image.open(os.path.join(RAW, f"f_{i:04d}.png")).convert("RGB")
    canvas.paste(im, (0, HEADER_H))

    foot_y = HEADER_H + PH
    canvas.paste(cbar_img, (W // 2 - cbar_img.width // 2, foot_y + 14))

    canvas.save(os.path.join(FRAMES_DIR, f"frame_{i:04d}.png"))
    if i % 20 == 0 or i == len(frames) - 1:
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
