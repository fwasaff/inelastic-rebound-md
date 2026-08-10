"""
render_corte_transversal.py — Video OVITO de un CORTE TRANSVERSAL por el punto
de impacto, coloreado por energia cinetica en escala logaritmica.

Para que sirve: es la unica vista donde se ven a la vez las tres cosas que
persigue la linea onda/calor, porque ocurren DENTRO del material y un render 3D
las esconde:

  1. la onda longitudinal bajando hacia el volumen (la "onda de choque"),
  2. la onda de Rayleigh corriendo por la superficie (mas lenta, ~0.55 c_L,
     confinada a ~1 longitud de onda de profundidad),
  3. la estela caliente que queda detras.

Por que escala LOGARITMICA: la KE por atomo va de ~1e-21 (red intacta) a ~1e-3
(frente) a ~1 (bola). En escala lineal solo se ve la bola y todo lo demas es
negro. En log se ven los tres fenomenos a la vez.

Por que un CORTE y no el volumen entero: con ~1e6 atomos, la superficie exterior
tapa el interior. El corte de unos pocos sigma de espesor por el eje del impacto
deja ver la propagacion en profundidad, que es donde esta la fisica.

Requiere el LD_LIBRARY_PATH de PySide6 (choque de Qt con el del sistema):

    LD_LIBRARY_PATH=~/ovito_env/lib/python3.12/site-packages/PySide6/Qt/lib \\
    ~/ovito_env/bin/python scripts/render_corte_transversal.py \\
        --run data/outputs/activo/ondas_N79_profunda_h28.15

Se puede correr con `nice -n 15` para no competir con una simulacion en curso.
"""

import argparse
import os
import subprocess
import glob

import numpy as np

from ovito.io import import_file
from ovito.modifiers import (SliceModifier, ComputePropertyModifier,
                             ColorCodingModifier, ExpressionSelectionModifier,
                             DeleteSelectedModifier)
from ovito.vis import (Viewport, TachyonRenderer, TextLabelOverlay,
                       ColorLegendOverlay)

p = argparse.ArgumentParser()
p.add_argument("--run", required=True, help="directorio de la corrida")
p.add_argument("--out", default=None, help="mp4 de salida")
p.add_argument("--modo", default="ke", choices=["ke", "vz", "chi"],
               help="ke  = log10 energia cinetica (cuanta energia hay)\n"
                    "vz  = velocidad vertical con mapa divergente (la ONDA se ve "
                    "como bandas suaves; el CALOR como moteado atomo a atomo)\n"
                    "chi = fraccion coherente por atomo (rojo=onda, azul=calor). "
                    "Es la vista que separa frente de estela, porque la energia "
                    "cinetica NO los distingue: ambos tienen energia, difieren en "
                    "organizacion")
p.add_argument("--ke-corte", type=float, default=1e-6,
               help="en modo chi/vz, se ocultan los atomos con KE por debajo de "
                    "esto: aun no les ha llegado nada y su 'coherencia' es ruido")
p.add_argument("--espesor", type=float, default=6.0,
               help="espesor del corte [sigma]. Fino se ve mejor la onda pero "
                    "hay menos atomos y la imagen queda rala")
p.add_argument("--ke-min", type=float, default=-7.0, help="log10(KE) minimo del color")
p.add_argument("--ke-max", type=float, default=-3.0, help="log10(KE) maximo del color")
p.add_argument("--zmax", type=float, default=72.0, help="recorte superior de la vista")
p.add_argument("--fps", type=int, default=12)
p.add_argument("--ancho", type=int, default=1280)
p.add_argument("--alto", type=int, default=720)
p.add_argument("--vmax", type=float, default=0.008,
               help="modo vz: rango simetrico del color [sigma/tau]")
p.add_argument("--cutoff", type=float, default=3.0,
               help="modo chi: radio del promedio local [sigma]")
p.add_argument("--radio", type=float, default=0.45, help="radio de particula [sigma]")
p.add_argument("--frames", type=int, default=0, help="0 = todos; N = solo N frames")
p.add_argument("--frame-inicio", type=int, default=0, help="primer frame a renderizar")
p.add_argument("--bounce", default=None, help="bounce.dat, para el instante del impacto")
p.add_argument("--dt", type=float, default=0.005)
args = p.parse_args()

TRAJ = os.path.join(args.run, "traj.lammpstrj")
OUT = args.out or os.path.join(args.run, f"corte_{args.modo}.mp4")
TMP = "/tmp/claude-1000/corte_frames"
os.makedirs(TMP, exist_ok=True)
for f in glob.glob(f"{TMP}/*.png"):
    os.remove(f)

print(f"Cargando {TRAJ} ...")
pipe = import_file(TRAJ)
n_frames = pipe.source.num_frames
print(f"  {n_frames} frames")

# --- Centro del impacto: centroide de la bola (tipo 2) en el primer frame ------
d0 = pipe.compute(0)
tipos = np.asarray(d0.particles["Particle Type"])
pos = np.asarray(d0.particles["Position"])
cx, cy = pos[tipos == 2][:, :2].mean(axis=0)
print(f"  centro del impacto: ({cx:.2f}, {cy:.2f}) sigma")

# --- Corte: losa delgada perpendicular a Y que pasa por el impacto ------------
pipe.modifiers.append(SliceModifier(normal=(0.0, 1.0, 0.0), distance=float(cy),
                                    slab_width=args.espesor))

# Recorte superior: por encima de zmax solo hay vacio (y la bola cayendo desde
# muy arriba), que comprimiria la escala de la vista sin aportar nada.
pipe.modifiers.append(ExpressionSelectionModifier(expression=f"Position.Z > {args.zmax}"))
pipe.modifiers.append(DeleteSelectedModifier())

# --- Campo de color segun el modo ----------------------------------------------
if args.modo == "ke":
    # El +1e-30 evita log(0) en los atomos que aun no se han enterado del
    # impacto (KE exactamente 0 en la red intacta a T=0).
    pipe.modifiers.append(ComputePropertyModifier(
        output_property="logKE",
        expressions=["log(c_keratom + 1e-30)/log(10)"]))
    prop, v0, v1 = "logKE", args.ke_min, args.ke_max
    grad, titulo = ColorCodingModifier.Viridis(), "log10 KE por atomo"

else:
    # Los atomos a los que aun no ha llegado nada tienen KE ~1e-21 y una
    # "coherencia" que es puro ruido numerico. Se ocultan.
    pipe.modifiers.append(ExpressionSelectionModifier(
        expression=f"c_keratom < {args.ke_corte}"))
    pipe.modifiers.append(DeleteSelectedModifier())

    if args.modo == "vz":
        # Mapa DIVERGENTE sobre la velocidad vertical. Una onda coherente se ve
        # como bandas amplias y suaves de un signo (compresion / rarefaccion);
        # el calor, como moteado que cambia de signo entre atomos vecinos. La
        # diferencia salta a la vista sin necesidad de calcular nada.
        prop, v0, v1 = "Velocity.Z", -args.vmax, args.vmax
        grad, titulo = ColorCodingModifier.BlueWhiteRed(), "v_z  [sigma/tau]"
    else:
        # chi por atomo: misma descomposicion que el analisis. Se compara la
        # velocidad de cada atomo con la MEDIA LOCAL de sus vecinos:
        #   u = media local  -> movimiento colectivo = onda
        #   v - u            -> lo que se sale de esa media = calor
        #   chi = |u|^2 / (|u|^2 + |v-u|^2)
        # chi ~ 1 (rojo): el atomo se mueve CON sus vecinos -> coherente.
        # chi ~ 0 (azul): se mueve por su cuenta -> termalizado.
        pipe.modifiers.append(ComputePropertyModifier(
            output_property="vsum", expressions=["0", "0", "0"],
            neighbor_expressions=["Velocity.X", "Velocity.Y", "Velocity.Z"],
            cutoff_radius=args.cutoff))
        pipe.modifiers.append(ComputePropertyModifier(
            output_property="nnb", expressions=["0"],
            neighbor_expressions=["1"], cutoff_radius=args.cutoff))
        UX = "((vsum.1+Velocity.X)/(nnb+1))"
        UY = "((vsum.2+Velocity.Y)/(nnb+1))"
        UZ = "((vsum.3+Velocity.Z)/(nnb+1))"
        U2 = f"({UX}^2+{UY}^2+{UZ}^2)"
        D2 = f"((Velocity.X-{UX})^2+(Velocity.Y-{UY})^2+(Velocity.Z-{UZ})^2)"
        pipe.modifiers.append(ComputePropertyModifier(
            output_property="chi", expressions=[f"{U2}/({U2}+{D2}+1e-30)"]))
        prop, v0, v1 = "chi", 0.0, 1.0
        grad, titulo = ColorCodingModifier.BlueWhiteRed(), "chi  (1=onda, 0=calor)"

pipe.modifiers.append(ColorCodingModifier(
    property=prop, start_value=v0, end_value=v1, gradient=grad))

data = pipe.compute(0)
data.particles.vis.radius = args.radio
print(f"  atomos en el corte: {data.particles.count}")

# La celda de simulacion llega hasta z=100 (espacio de caida de la bola) aunque
# ahi no quede ningun atomo. Dibujarla obliga a encuadrar 100 sigma de los que
# 30 son vacio, y aplasta la losa contra el borde inferior. Se apaga.
data.cell.vis.enabled = False

# El pipeline TIENE que estar en la escena para que el viewport lo vea; sin
# esto OVITO renderiza fondo vacio y solo avisa con un UserWarning.
pipe.add_to_scene()

# --- Camara: vista frontal ortografica (plano x-z) -----------------------------
# Encuadre FIJO y explicito, no zoom_all: asi el encuadre no cambia entre frames
# (zoom_all seguiria a la bola) y la escala espacial es comparable a lo largo de
# todo el video. fov = semialtura en unidades del mundo.
vp = Viewport(type=Viewport.Type.Front)
vp.camera_pos = (cx, cy, args.zmax / 2.0)
vp.fov = args.zmax / 2.0

# --- Reloj fisico: t - t_imp ---------------------------------------------------
# El video se anota con el tiempo RELATIVO AL IMPACTO, no con el paso de LAMMPS:
# es la unica referencia con la que se pueden leer las figuras del analisis.
bounce = args.bounce or os.path.join(args.run, "bounce.dat")
bd = np.loadtxt(bounce, comments="#")
t_imp = float(bd[np.argmin(bd[:, 1]), 0])
print(f"  t_imp = {t_imp:.2f} tau")

reloj = TextLabelOverlay(text="", alignment=0x21, offset_y=-0.04, offset_x=0.02,
                         font_size=0.035, text_color=(0, 0, 0))
vp.overlays.append(reloj)
# Leyenda abajo-izquierda (0x41 = AlignBottom|AlignLeft) y en BLANCO: ahi el
# fondo es la red intacta, que en esta escala de color es morado oscuro. Arriba
# no cabe -- la franja blanca sobre la losa mide solo ~10 sigma y el titulo de
# la leyenda se salia del encuadre.
vp.overlays.append(ColorLegendOverlay(
    modifier=pipe.modifiers[-1], alignment=0x41, orientation=0x1,
    legend_size=0.26, font_size=0.022, format_string="%.0f",
    offset_x=0.03, offset_y=0.05,
    title=titulo, text_color=(1, 1, 1)))

f0 = args.frame_inicio
f1 = n_frames - 1 if args.frames == 0 else min(f0 + args.frames - 1, n_frames - 1)
print(f"Renderizando frames {f0}..{f1} ...")
# Frame a frame (no render_anim) porque el rotulo del reloj cambia en cada uno.
rend = TachyonRenderer(shadows=False, direct_light_intensity=1.1)
for k in range(f0, f1 + 1):
    dk = pipe.compute(k)
    tk = dk.attributes["Timestep"] * args.dt - t_imp
    reloj.text = f"t - t_imp = {tk:+.2f} tau"
    vp.render_image(filename=f"{TMP}/f{k:04d}.png", size=(args.ancho, args.alto),
                    frame=k, renderer=rend, background=(1, 1, 1))
    if (k - f0) % 10 == 0:
        print(f"    frame {k}  (t-t_imp = {tk:+.2f} tau)")

pngs = sorted(glob.glob(f"{TMP}/*.png"))
print(f"  {len(pngs)} PNG generados")
if not pngs:
    raise SystemExit("OVITO no genero frames")

print(f"Codificando {OUT} ...")
subprocess.run(["ffmpeg", "-y", "-framerate", str(args.fps),
                "-pattern_type", "glob", "-i", f"{TMP}/*.png",
                "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20",
                OUT], check=True, capture_output=True)
print(f"Listo: {OUT}  ({os.path.getsize(OUT)/1e6:.1f} MB)")
