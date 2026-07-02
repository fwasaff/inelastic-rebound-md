"""
gen_paper_snapshots.py  —  Genera dos imágenes para el paper
  1. fig_setup.png    — bola fría sobre superficie FCC (§ Modelo)
  2. fig_impact.png   — instante de máxima compresión, KE Hot solo en la bola (§ Resultados)

Uso (desde el directorio raíz del proyecto):
    LD_LIBRARY_PATH=~/ovito_env/lib/python3.12/site-packages/PySide6/Qt/lib:$LD_LIBRARY_PATH \
    ~/ovito_env/bin/python scripts/gen_paper_snapshots.py
"""

import os
import numpy as np

TRAJ       = "data/outputs/reprod_2012/h0_41.47/traj.lammpstrj"
OUT_SETUP  = "data/outputs/reprod_2012/fig_setup.png"
OUT_IMPACT = "data/outputs/reprod_2012/fig_impact.png"

FRAME_SETUP  = 175   # t = 43.75 τ  —  bola a ~16σ sobre la superficie
FRAME_IMPACT = 224   # t = 56.0 τ   —  compresión máxima

WIDTH_SETUP  = 1400
HEIGHT_SETUP = 1050  # 4:3

WIDTH_IMPACT  = 1400
HEIGHT_IMPACT = 1050

BG = (0.92, 0.94, 0.96)   # gris claro, estándar revistas de física

from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer
from ovito.modifiers import ColorCodingModifier

print("Cargando trayectoria...")
pipeline = import_file(TRAJ, sort_particles=True)
print(f"  {pipeline.source.num_frames} frames, "
      f"{pipeline.compute(0).particles.count} átomos")

renderer = TachyonRenderer()
renderer.ambient_occlusion        = True
renderer.shadows                  = True
renderer.direct_light_intensity   = 1.6
renderer.ambient_occlusion_brightness = 1.3

# Geometría de la caja
data0 = pipeline.compute(0)
pos0  = np.array(data0.particles['Position'])
xc    = pos0[:, 0].mean()
yc    = pos0[:, 1].mean()

# ===========================================================================
# IMAGEN 1 — Setup (coloring por tipo: superficie gris-azul, bola naranja)
# ===========================================================================
print("\n[1/2] Renderizando fig_setup (frame {}, t={:.1f}τ)...".format(
      FRAME_SETUP, FRAME_SETUP * 0.25))

# Quitar cualquier modificador previo
pipeline.modifiers.clear()

# Modificador Python: asigna color por tipo de partícula
def color_by_type(frame, data):
    types  = np.array(data.particles['Particle Type'])
    colors = np.zeros((len(types), 3))
    colors[types == 1] = [0.50, 0.62, 0.75]   # superficie: azul-gris acero
    colors[types == 2] = [0.95, 0.42, 0.08]   # bola: naranja cálido
    data.particles_.create_property('Color', data=colors)

pipeline.modifiers.append(color_by_type)
pipeline.add_to_scene()
pipeline.source.data.particles_.vis.radius = 0.45
pipeline.source.data.cell_.vis.enabled     = False

# Posiciones en el frame de setup
data_s = pipeline.compute(FRAME_SETUP)
pos_s  = np.array(data_s.particles['Position'])
typ_s  = np.array(data_s.particles['Particle Type'])
z_ball_max = pos_s[typ_s == 2, 2].max()   # tope de la bola
z_surf_min = pos_s[typ_s == 1, 2].min()   # fondo de la superficie
z_mid      = (z_ball_max + z_surf_min) / 2   # centro vertical

# Cámara vista lateral: mirando en +Y desde -Y, centrada verticalmente
D_setup = 42.0
vp_setup = Viewport()
vp_setup.type       = Viewport.Type.Perspective
vp_setup.camera_pos = (xc, yc - D_setup, z_mid)
vp_setup.camera_dir = (0.0, 1.0, 0.0)   # horizontal puro
vp_setup.camera_up  = (0.0, 0.0, 1.0)
vp_setup.fov        = 0.82   # ~47° — encuadra toda la escena verticalmente

vp_setup.render_image(
    size=(WIDTH_SETUP, HEIGHT_SETUP),
    frame=FRAME_SETUP,
    renderer=renderer,
    background=BG,
    filename=OUT_SETUP,
)
print(f"  → {OUT_SETUP}  ({os.path.getsize(OUT_SETUP)//1024} KB)")

# ===========================================================================
# IMAGEN 2 — Impacto (coloring por c_keratom, mapa Hot)
# ===========================================================================
print("\n[2/2] Renderizando fig_impact (frame {}, t={:.1f}τ)...".format(
      FRAME_IMPACT, FRAME_IMPACT * 0.25))

# Reiniciar pipeline: superficie gris fija + bola coloreada por KE (mapa Hot)
pipeline.modifiers.clear()

def color_mixed(frame, data):
    types  = np.array(data.particles['Particle Type'])
    ke     = np.array(data.particles['c_keratom'])
    colors = np.zeros((len(types), 3))
    # Superficie: gris-azul frío uniforme
    colors[types == 1] = [0.55, 0.62, 0.70]
    # Bola: mapa Hot  0→negro  0.4→rojo  0.8→naranja  1→blanco  (escala 0..0.5ε)
    def hot(v):
        v = np.clip(v, 0.0, 1.0)
        r = np.clip(v / 0.4, 0, 1)
        g = np.clip((v - 0.4) / 0.4, 0, 1)
        b = np.clip((v - 0.8) / 0.2, 0, 1)
        return np.stack([r, g, b], axis=1)
    ball = types == 2
    colors[ball] = hot(ke[ball] / 0.5)
    data.particles_.create_property('Color', data=colors)

pipeline.modifiers.append(color_mixed)

# Centro de escena en el frame de impacto
data_i  = pipeline.compute(FRAME_IMPACT)
pos_i   = np.array(data_i.particles['Position'])
typ_i   = np.array(data_i.particles['Particle Type'])
xci     = pos_i[:, 0].mean()
yci     = pos_i[:, 1].mean()
z_surf_top = pos_i[typ_i == 1, 2].max()
z_ball_top = pos_i[typ_i == 2, 2].max()
z_focus    = (z_surf_top + z_ball_top) / 2

# Cámara diagonal elevada: muestra superficie FCC + bola desde arriba-frente
vp_impact = Viewport()
vp_impact.type       = Viewport.Type.Perspective
vp_impact.camera_pos = (xci, yci - 32, z_focus + 14)
vp_impact.camera_dir = (0.0, 0.85, -0.53)   # ~32° hacia abajo
vp_impact.camera_up  = (0.0, 0.0, 1.0)
vp_impact.fov        = 0.60

vp_impact.render_image(
    size=(WIDTH_IMPACT, HEIGHT_IMPACT),
    frame=FRAME_IMPACT,
    renderer=renderer,
    background=BG,
    filename=OUT_IMPACT,
)
print(f"  → {OUT_IMPACT}  ({os.path.getsize(OUT_IMPACT)//1024} KB)")

print("\nListo. Imágenes guardadas en data/outputs/reprod_2012/")
