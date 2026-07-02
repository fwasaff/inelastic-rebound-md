"""
gen_energy_diagram.py  —  Diagrama de canales energéticos (versión diseño final)

Geometría calculada para figura 20×15 @ 150 dpi (3 000 × 2 250 px):
  ─ Título        y = 0.976  (13 pt → ~0.012 alto → base ≈ 0.964)
  ─ E₀            cy= 0.904, h=0.088  (top=0.948, bot=0.860)
  ─ Crossbar      y = 0.842  ((0.860+0.824)/2)
  ─ BG_TOP        y = 0.824
  ─ GHY           y = 0.806
  ─ TRUNK_Y       y = 0.784
  ─ Row 1 boxes   cy= 0.636, h=0.230  (top=0.751, bot=0.521)
  ─ Row 2 boxes   cy= 0.356, h=0.230  (top=0.471, bot=0.241)
  ─ BG_BOT        y = 0.218
  ─ LJ formula    y = 0.210
  ─ Bottom panels y = 0.010, h=0.192  (top=0.202)
"""

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

OUT = "data/outputs/reprod_2012/fig_energy_diagram.png"

# ── Paleta ────────────────────────────────────────────────────────────────────
C = dict(
    bg      = '#ECEEF2',
    e0      = '#1B2631',
    cm      = '#1A4A6B',   # azul profundo
    tb      = '#7B241C',   # rojo oscuro
    dvb     = '#A04000',   # naranja quemado
    uf      = '#7D6608',   # dorado oscuro
    ks      = '#1A6B3A',   # verde bosque
    dvs     = '#0D5C52',   # teal oscuro
    uc      = '#512E5F',   # violeta oscuro
    ball_bg = '#FEF9E7',
    ball_bd = '#922B21',
    surf_bg = '#E8F8F5',
    surf_bd = '#0D5C52',
    line    = '#5D6D7E',
    dark    = '#1A1A2E',
    mid     = '#484B4B',
    light   = '#717D7E',
)

# ── Figura ────────────────────────────────────────────────────────────────────
FW, FH = 20, 15
fig = plt.figure(figsize=(FW, FH))
ax  = fig.add_axes([0.012, 0.010, 0.976, 0.980])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')
fig.patch.set_facecolor(C['bg'])
ax.set_facecolor(C['bg'])

# ── Primitivas ────────────────────────────────────────────────────────────────
def R(x0, y0, w, h, fc, ec, lw=2.0, ls='-', r=0.010, z=2, a=1.0):
    ax.add_patch(FancyBboxPatch((x0, y0), w, h,
        boxstyle=f'round,pad={r}',
        facecolor=fc, edgecolor=ec, linewidth=lw,
        linestyle=ls, alpha=a, zorder=z))

def T(x, y, s, **kw):
    kw.setdefault('zorder', 9)
    ax.text(x, y, s, **kw)

def HL(x0, x1, y, lw=1.6):
    ax.plot([x0, x1], [y, y], color=C['line'], lw=lw,
            solid_capstyle='round', zorder=6)

def VL(x, y0, y1, lw=1.6):
    ax.plot([x, x], [y0, y1], color=C['line'], lw=lw,
            solid_capstyle='round', zorder=6)

def ARR(x, y_from, y_to, lw=1.6):
    ax.annotate('', xy=(x, y_to), xytext=(x, y_from),
        arrowprops=dict(arrowstyle='-|>', color=C['line'],
                        lw=lw, mutation_scale=10), zorder=6)

# ── Caja de canal energético ──────────────────────────────────────────────────
def canal(cx, cy, w, h, col,
          hdr,          # texto encabezado
          sym,          # símbolo LaTeX (string sin $)
          eq1,          # ecuación principal
          eq2,          # ecuación secundaria (o None)
          note,         # nota 1 línea (o None)
          pct,          # "XX %" para badge
          tsz=10.0, ssz=18.0, eqsz=9.5, nsz=8.8):
    x0, y0 = cx - w/2, cy - h/2

    # Fondo blanco + borde coloreado
    R(x0, y0, w, h, 'white', col, lw=2.8, r=0.014, z=4)

    # Encabezado coloreado (top 24 %)
    hh = h * 0.24
    R(x0, y0+h-hh, w, hh, col, col, lw=0, r=0.006, z=5)
    T(cx, y0+h-hh/2, hdr,
      ha='center', va='center', fontsize=tsz,
      fontweight='bold', color='white')

    # Símbolo grande  (66 % desde abajo)
    T(cx, y0+h*0.660, f'${sym}$',
      ha='center', va='center', fontsize=ssz,
      fontweight='bold', color=col)

    # Ecuación principal  (49 %)
    T(cx, y0+h*0.490, eq1,
      ha='center', va='center', fontsize=eqsz, color=C['dark'])

    # Ecuación secundaria  (36 %)
    if eq2:
        T(cx, y0+h*0.360, eq2,
          ha='center', va='center', fontsize=eqsz-0.5,
          color=C['mid'], style='italic')

    # Nota física  (22 %)  —  espacio amplio antes del badge
    if note:
        T(cx, y0+h*0.220, note,
          ha='center', va='center', fontsize=nsz,
          color=C['mid'])

    # Badge  (9 %)
    T(cx, y0+h*0.090, f'▸ Apogeo: {pct}',
      ha='center', va='center', fontsize=9.0,
      color='white', fontweight='bold',
      bbox=dict(facecolor=col, alpha=0.93,
                edgecolor='none', pad=3.0, boxstyle='round'))


# ═══════════════════════════════════════════════════════════════════════════════
# CONSTANTES DE GEOMETRÍA
# ═══════════════════════════════════════════════════════════════════════════════
E0_CX, E0_CY = 0.500, 0.904
E0_W,  E0_H  = 0.295, 0.088
E0_BOT        = E0_CY - E0_H/2       # 0.860

BG_TOP  = 0.824
BG_BOT  = 0.208
GHY     = 0.806                       # y-centro encabezados de grupo
BDIV    = 0.550                       # divisoria grupos
BOLA_X  = 0.275                       # x árbol bola
SURF_X  = 0.773                       # x árbol superficie

XBAR_Y  = (E0_BOT + BG_TOP) / 2      # 0.842
TRUNK_Y = GHY - 0.022                 # 0.784

BW, BH = 0.232, 0.255                 # cajas bola (más altas → más espacio interno)
SW, SH = 0.390, 0.255                 # cajas superficie
BX1, BX2 = 0.126, 0.394              # x centros col-izq y col-der bola
RY1, RY2 = 0.636, 0.353              # y centros fila-1 y fila-2
SY1, SY2 = RY1, RY2                  # alineadas con bola

BRANCH_Y1 = RY1 + BH/2 + 0.014      # 0.636+0.128+0.014 = 0.778  < TRUNK_Y ✓
BRANCH_Y2 = RY2 + BH/2 + 0.014      # 0.353+0.128+0.014 = 0.495

TRUNK_BX  = (BX1 + BX2) / 2         # 0.260  (entre columnas)

# ═══════════════════════════════════════════════════════════════════════════════
# TÍTULO
# ═══════════════════════════════════════════════════════════════════════════════
T(0.500, 0.976,
  "Canales energéticos en el impacto de un nanocúmulo LJ"
  "   —   valores al primer apogeo  (h₀ = 41 σ,  N = 87 átomos)",
  ha='center', va='top', fontsize=13, fontweight='bold', color=C['dark'])

# ═══════════════════════════════════════════════════════════════════════════════
# FONDOS DE GRUPO
# ═══════════════════════════════════════════════════════════════════════════════
R(0.012, BG_BOT, BDIV-0.014, BG_TOP-BG_BOT,
  C['ball_bg'], C['ball_bd'], lw=2.2, r=0.014, z=1)
R(BDIV+0.002, BG_BOT, 0.986-BDIV, BG_TOP-BG_BOT,
  C['surf_bg'], C['surf_bd'], lw=2.2, r=0.014, z=1)

T(BOLA_X, GHY, "BOLA  (nanocúmulo LJ)",
  ha='center', va='center', fontsize=12.5,
  fontweight='bold', color=C['ball_bd'], zorder=3)
T(SURF_X, GHY, "SUPERFICIE  (cristal FCC fijo)",
  ha='center', va='center', fontsize=12.5,
  fontweight='bold', color=C['surf_bd'], zorder=3)

ax.plot([BDIV, BDIV], [BG_BOT+0.012, BG_TOP-0.012],
        color='#B0B8C1', lw=0.9, ls='--', zorder=2)

# ═══════════════════════════════════════════════════════════════════════════════
# E₀  (raíz)
# ═══════════════════════════════════════════════════════════════════════════════
canal(E0_CX, E0_CY, E0_W, E0_H, C['e0'],
      "ENERGÍA TOTAL  E₀",
      r"E_0 = M \cdot g \cdot h_0",
      r"$E_0 = N_{ball}\cdot g\cdot(z_0 - z_{cont}) \approx 90\,\varepsilon$",
      None, None, "100 %",
      tsz=12.0, ssz=15.0, eqsz=10.0)

# ═══════════════════════════════════════════════════════════════════════════════
# ÁRBOL DE CONEXIONES
# ═══════════════════════════════════════════════════════════════════════════════
# E₀ → crossbar → grupos
VL(E0_CX, E0_BOT, XBAR_Y, lw=2.0)
HL(BOLA_X, SURF_X, XBAR_Y, lw=2.0)
ARR(BOLA_X, XBAR_Y, BG_TOP, lw=2.0)
ARR(SURF_X, XBAR_Y, BG_TOP, lw=2.0)

# BOLA: tronco vertical + ramas horizontales en T
VL(TRUNK_BX, TRUNK_Y, BRANCH_Y2, lw=1.7)
HL(BX1, BX2, BRANCH_Y1, lw=1.4)          # rama fila 1
ARR(BX1, BRANCH_Y1, RY1+BH/2, lw=1.4)
ARR(BX2, BRANCH_Y1, RY1+BH/2, lw=1.4)
HL(BX1, BX2, BRANCH_Y2, lw=1.4)          # rama fila 2
ARR(BX1, BRANCH_Y2, RY2+BH/2, lw=1.4)
ARR(BX2, BRANCH_Y2, RY2+BH/2, lw=1.4)

# SUPERFICIE: tronco + dos ramas
VL(SURF_X, TRUNK_Y, BRANCH_Y2, lw=1.7)
ARR(SURF_X, BRANCH_Y1, SY1+SH/2, lw=1.4)
ARR(SURF_X, BRANCH_Y2, SY2+SH/2, lw=1.4)

# ═══════════════════════════════════════════════════════════════════════════════
# 6 CAJAS DE CANALES
# ═══════════════════════════════════════════════════════════════════════════════
# ── BOLA fila 1 ──────────────────────────────────────────────────────────────
canal(BX1, RY1, BW, BH, C['cm'],
      "Energía cinética CM",
      r"K_{CM}",
      r"$K_{CM} = \frac{1}{2}\,M\,V_{CM}^{2}$",
      r"$M = N_{ball}\cdot m$",
      "Bola volando como proyectil",
      "0 %  (apogeo)", tsz=9.5)

canal(BX2, RY1, BW, BH, C['tb'],
      "Energía calorífica",
      r"T_B",
      r"$T_B = \sum_i\frac{1}{2}m\,|v_i-V_{CM}|^2$",
      r"$k_B\,T_{bola} = \frac{2\,T_B}{3(N-1)}$",
      "KE relativa al CM → temperatura",
      "20 %", tsz=9.5, eqsz=9.2)

# ── BOLA fila 2 ──────────────────────────────────────────────────────────────
canal(BX1, RY2, BW, BH, C['dvb'],
      "Energía configuracional",
      r"\Delta V_B",
      r"$\Delta V_B = V_B(t) - V_B(0)$",
      r"$V_B = \sum_{i<j\in bola}\varphi(r_{ij})$",
      "PE intra-bola cambia (anarmónico)",
      "25 %", tsz=9.5, eqsz=9.2)

canal(BX2, RY2, BW, BH, C['uf'],
      "Energía gravitatoria",
      r"U_F",
      r"$U_F = M \cdot g \cdot z_{CM}$",
      None,
      "Recuperable en el próximo rebote",
      "43 %", tsz=9.5)

# ── SUPERFICIE fila 1 ─────────────────────────────────────────────────────────
canal(SURF_X, SY1, SW, SH, C['ks'],
      "Energía calorífica  (fonones cin.)",
      r"K_S",
      r"$K_S = \sum_{i\in sup}\frac{1}{2}m\,v_i^{\,2}$",
      r"$k_B\,T_{sup} = \frac{2\,K_S}{3\,N_{sup}}$",
      "Sin V_CM libre en sólido → K_S = calor",
      "6 %", tsz=9.5, eqsz=9.2)

# ── SUPERFICIE fila 2 ─────────────────────────────────────────────────────────
canal(SURF_X, SY2, SW, SH, C['dvs'],
      "Energía configuracional  (fonones pot.)",
      r"\Delta V_S",
      r"$\Delta V_S = V_S(t) - V_S(0)$",
      r"$V_S = \sum_{i<j\in sup}\varphi(r_{ij})$",
      "Red cristalina deformada post-impacto",
      "5 %", tsz=9.5, eqsz=9.2)

# ═══════════════════════════════════════════════════════════════════════════════
# SEPARADOR DECORATIVO
# ═══════════════════════════════════════════════════════════════════════════════
ax.plot([0.018, 0.982], [0.205, 0.205],
        color='#C0C5CE', lw=0.8, zorder=3)

# ═══════════════════════════════════════════════════════════════════════════════
# SECCIÓN INFERIOR  —  3 paneles
# ═══════════════════════════════════════════════════════════════════════════════
PY0 = 0.008
PH  = 0.192

# ── Panel izq: respuesta sobre el calor (40 %) ───────────────────────────────
R(0.010, PY0, 0.388, PH, '#EBF5FB', '#1F618D', lw=1.6, r=0.010, z=3)

def Tp(x, dy, s, **kw):   # texto dentro de panel (dy desde tope)
    kw.setdefault('zorder', 9)
    ax.text(x, PY0+PH-dy, s, ha='center', va='top', **kw)

XI = 0.204   # centro panel izq

Tp(XI, 0.012, "¿El calor cinético = KE relativa al CM?",
   fontsize=10.8, fontweight='bold', color='#1F618D')
Tp(XI, 0.036, "SÍ — para la componente cinética del calor interno.",
   fontsize=9.2, color=C['dark'])
Tp(XI, 0.060,
   r"Bola:  $T_B = \sum\frac{1}{2}m\,|v_i-V_{CM}|^2$"
   r"   (KE en el marco del CM → temperatura)",
   fontsize=9.0, color=C['tb'])
Tp(XI, 0.083,
   r"Superficie:  $K_S = \sum\frac{1}{2}m\,v_i^2$"
   r"   (no hay V$_{CM}$ libre → toda KE es calor)",
   fontsize=9.0, color=C['ks'])
Tp(XI, 0.107,
   "Calor TOTAL (termodinámica) incluye también la PE interna:",
   fontsize=8.8, color=C['mid'], style='italic')
Tp(XI, 0.127,
   r"$E_{calor}^{bola} \approx T_B + \Delta V_B$"
   r"          "
   r"$E_{calor}^{sup} \approx K_S + \Delta V_S$  (fonones)",
   fontsize=9.0, color=C['dark'])
Tp(XI, 0.150,
   "En equilibrio (equipartición): ⟨KE⟩ ≈ ⟨PE⟩ → calor se reparte igual.",
   fontsize=8.8, color=C['mid'], style='italic')
Tp(XI, 0.172,
   r"Potencial:  $\varphi(r)=4\varepsilon\left[(\sigma/r)^{12}-c\,(\sigma/r)^{6}\right]$"
   r"   c=1 (bola–bola, sup.–sup.)   c=0 (bola–sup.)",
   fontsize=8.5, color=C['light'])

# ── Panel centro: U_c (27 %) ──────────────────────────────────────────────────
R(0.406, PY0, 0.262, PH, '#F4ECF7', C['uc'], lw=1.6, ls='--', r=0.010, z=3)

XC = 0.537

Tp(XC, 0.012, "⚡  Uᴄ — PE de colisión  (TRANSITORIA)",
   fontsize=10.5, fontweight='bold', color=C['uc'])
Tp(XC, 0.035, "Existe solo mientras hay contacto bola–sup.",
   fontsize=8.8, color=C['mid'], style='italic')
Tp(XC, 0.060,
   r"$U_c = \sum_{i\in bola,\,j\in sup}\varphi_{rep}(r_{ij})$",
   fontsize=9.5, color=C['dark'])
Tp(XC, 0.088,
   "Máxima compresión  (t ≈ 56τ):",
   fontsize=8.8, color=C['mid'])
Tp(XC, 0.107,
   "Uᴄ  ≈  60 % de E₀",
   fontsize=11.5, fontweight='bold', color=C['uc'])
Tp(XC, 0.133, "Al separarse  (t > 62τ):  Uᴄ → 0",
   fontsize=8.8, color=C['mid'])
Tp(XC, 0.152, "Se libera como  ΔV_B + ΔV_S",
   fontsize=9.2, color=C['dark'])
Tp(XC, 0.170,
   "(c_peratom divide Uᴄ/2 a bola, Uᴄ/2 a superficie)",
   fontsize=8.2, color=C['light'], style='italic')

# ── Panel der: conservación (31 %) ───────────────────────────────────────────
R(0.676, PY0, 0.314, PH, '#EAFAF1', '#1E8449', lw=1.6, r=0.010, z=3)

XR = 0.833

Tp(XR, 0.012, "✓  Conservación de energía  (NVE)",
   fontsize=10.8, fontweight='bold', color='#1E8449')
Tp(XR, 0.042,
   r"$K_{CM} + T_B + \Delta V_B + K_S + \Delta V_S + U_F$",
   fontsize=10.5, color=C['dark'])

ax.plot([0.685, 0.981], [PY0+PH-0.068, PY0+PH-0.068],
        color='#1E8449', lw=1.2, zorder=9)

Tp(XR, 0.080,
   r"$= \; E_0 \; = \; \mathrm{cte}$",
   fontsize=14, fontweight='bold', color='#1E8449')
Tp(XR, 0.108, "(durante contacto: sumar Uᴄ)",
   fontsize=8.8, color=C['light'], style='italic')
Tp(XR, 0.132,
   "Derivado de los computes LAMMPS:",
   fontsize=8.8, color=C['mid'])
Tp(XR, 0.150,
   "KB3d + TB + VB + KS + PE_surf + UF ≈ E₀",
   fontsize=9.0, color=C['dark'], family='monospace')
Tp(XR, 0.170,
   "Verificado:  suma = 100 % en cada frame",
   fontsize=9.0, fontweight='bold', color='#1E8449')

# ═══════════════════════════════════════════════════════════════════════════════
# GUARDADO
# ═══════════════════════════════════════════════════════════════════════════════
plt.savefig(OUT, dpi=150, bbox_inches='tight',
            facecolor=C['bg'], edgecolor='none')
print(f"✓  {OUT}   ({os.path.getsize(OUT)//1024} KB)")
