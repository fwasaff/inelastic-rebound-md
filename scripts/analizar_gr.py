"""
analizar_gr.py  —  Función de distribución de pares g(r) de la bola

Reproduce el análisis de la Fig. 8 de Gutiérrez et al. 2012:
muestra cómo la estructura FCC del nanocúmulo cambia (o colapsa) en función
de la energía del impacto, evaluando g(r) en cuatro instantes clave:
  t=0         — FCC perfecto inicial
  t=t_imp     — durante la compresión máxima
  t=t_apo1    — primer apogeo post-impacto
  t=t_apo2    — segundo apogeo (si la bola sigue rebotando)

Requiere (por cada h0):
    data/outputs/terminado/reprod_2012/h0_{h0}/ball_dump.lammpstrj   (dump generado por LAMMPS)
    data/outputs/terminado/reprod_2012/h0_{h0}/bounce.dat             (para identificar tiempos)

Uso (desde la raíz del proyecto):
    python3 scripts/analizar_gr.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import find_peaks

# ─── Configuración ────────────────────────────────────────────────────────────
BASE   = Path("data/outputs/terminado/reprod_2012")
OUTFIG = BASE / "fig_gr_ball.pdf"

# Casos representativos para la figura final
H0_PANELS = [11.03, 28.15, 79.53]

G          = 0.026      # gravedad [σ/τ²]
RBALL      = 2.64       # radio estimado de la bola [σ]
A0         = 1.5496     # parámetro de red FCC argón LJ [σ]
DR         = 0.025      # ancho de bin g(r) [σ]
RMAX       = 2.55       # r máximo (< cutoff 2.5σ evita artefactos)
DT         = 0.005      # paso de tiempo [τ]
DUMP_EVERY = 100        # pasos entre frames del dump → 0.5 τ/frame

# Posiciones teóricas de vecinos en FCC: a0/√2, a0, a0√(3/2), a0√2
R_FCC = [A0 / np.sqrt(2), A0, A0 * np.sqrt(1.5), A0 * np.sqrt(2)]

# Paleta y estilos para los 4 instantes de tiempo
T_COLORS = ['#C62828', '#F57F17', '#1565C0', '#2E7D32']
T_LS     = ['-',        '-',       '-',        '-'      ]
T_LW     = [1.8,        1.8,       1.8,        1.8      ]
OFFSET   = 22.0   # separación vertical entre curvas [unidades de g(r)]


# ─── Parseo del dump ──────────────────────────────────────────────────────────

def parse_ball_dump(path):
    """
    Lee ball_dump.lammpstrj (columnas: id x y z).
    Devuelve {timestep: np.array(N, 3)}.
    """
    frames = {}
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if 'ITEM: TIMESTEP' not in line:
                continue
            ts = int(f.readline().strip())
            f.readline()               # ITEM: NUMBER OF ATOMS
            n  = int(f.readline().strip())
            # Avanzar hasta la línea ITEM: ATOMS (salta BOX BOUNDS)
            while True:
                hdr = f.readline()
                if not hdr or 'ITEM: ATOMS' in hdr:
                    break
            rows = []
            for _ in range(n):
                parts = f.readline().split()
                rows.append([float(parts[1]), float(parts[2]), float(parts[3])])
            frames[ts] = np.array(rows)
    return frames


# ─── Tiempos clave ────────────────────────────────────────────────────────────

def get_key_steps(h0):
    """
    Lee bounce.dat y devuelve 4 pasos LAMMPS clave (o None si no existe):
      [paso_inicial=0, paso_impacto, paso_apo1, paso_apo2]
    Cada tiempo τ se mapea al paso de dump más cercano.
    """
    path = BASE / f"h0_{h0}" / "bounce.dat"
    data = np.loadtxt(path, comments='#')
    t_arr = data[:, 0]
    Zcm   = data[:, 1]

    # Impacto: mínimo global de Zcm (primera mitad de la simulación)
    n_search = len(Zcm) // 3
    idx_imp  = np.argmin(Zcm[:n_search * 2])
    t_imp    = t_arr[idx_imp]

    # Apogeos: máximos de Zcm después del impacto
    maxs, _ = find_peaks(Zcm, prominence=1.0, distance=100)
    maxs    = maxs[maxs > idx_imp]
    t_apo1  = t_arr[maxs[0]] if len(maxs) > 0 else None
    t_apo2  = t_arr[maxs[1]] if len(maxs) > 1 else None

    def snap(t):
        if t is None:
            return None
        step = round(t / DT)
        return (step // DUMP_EVERY) * DUMP_EVERY

    return [0, snap(t_imp), snap(t_apo1), snap(t_apo2)]


def nearest_frame(frames_dict, target_step):
    """Devuelve el timestep disponible más cercano a target_step."""
    avail = sorted(frames_dict.keys())
    return min(avail, key=lambda s: abs(s - target_step))


# ─── Cálculo de g(r) ─────────────────────────────────────────────────────────

def compute_gr(pos, rmax=RMAX, dr=DR):
    """
    Función de distribución de pares para un cluster finito.

    Normalización: g(r) → cte en la región bulk del cluster.
    Para un cristal FCC perfecto el primer pico es g(r₁) >> 1.
    Para un líquido/amorfo los picos se ensanchan y bajan.

    Parámetros
    ----------
    pos : array (N, 3)   — posiciones absolutas (CM irrelevante por ser diferencias)
    """
    N    = len(pos)
    Vol  = (4.0 / 3.0) * np.pi * RBALL**3
    rho0 = N / Vol

    r_edges   = np.arange(0.0, rmax + dr, dr)
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
    counts    = np.zeros(len(r_centers), dtype=np.int64)

    for i in range(N):
        diffs = pos[i + 1:] - pos[i]
        dists = np.sqrt((diffs**2).sum(axis=1))
        idx   = (dists / dr).astype(int)
        mask  = (dists < rmax) & (idx < len(counts))
        np.add.at(counts, idx[mask], 2)

    shell_vol = 4.0 * np.pi * r_centers**2 * dr
    shell_vol[shell_vol == 0] = np.nan
    gr = counts / (N * rho0 * shell_vol)
    return r_centers, gr


# ─── Figura ───────────────────────────────────────────────────────────────────

def make_figure(panels_data):
    """
    panels_data : list of dicts
      { 'h0': float,
        'curves': [(label_str, r_array, gr_array), ...] }
    """
    n   = len(panels_data)
    fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 6.0), sharey=False)
    if n == 1:
        axes = [axes]

    fig.patch.set_facecolor('white')
    plt.rcParams.update({'font.size': 10, 'axes.grid': False})

    # Posiciones FCC a marcar (primeros 3 vecinos)
    fcc_r     = R_FCC[:3]
    fcc_label = [r'$r_1\!=\!a_0/\!\sqrt{2}$',
                 r'$r_2\!=\!a_0$',
                 r'$r_3\!=\!a_0\sqrt{3/2}$']

    for ax, pdata in zip(axes, panels_data):
        h0     = pdata['h0']
        curves = pdata['curves']   # [(label, r, gr), ...]
        n_c    = len(curves)

        for k, (label, r, gr) in enumerate(curves):
            yoff = (n_c - 1 - k) * OFFSET   # curva t=0 arriba, última abajo
            ax.plot(r, gr + yoff,
                    color=T_COLORS[k % len(T_COLORS)],
                    lw=T_LW[k % len(T_LW)],
                    ls=T_LS[k % len(T_LS)],
                    label=label)
            # Línea base de cada curva
            ax.axhline(yoff, color='#bbbbbb', lw=0.5, ls=':')

        # Marcadores de vecinos FCC
        ymax_plot = (n_c - 1) * OFFSET + 20
        for ri, rlab in zip(fcc_r, fcc_label):
            ax.axvline(ri, color='#666666', lw=0.9, ls='--', alpha=0.6)
            ax.text(ri, ymax_plot * 1.01, rlab,
                    ha='center', va='bottom', fontsize=7.5, color='#555555')

        ax.set_xlim(0.65, RMAX)
        ax.set_ylim(-2, ymax_plot * 1.10)
        ax.set_xlabel(r'$r\;(\sigma)$', fontsize=12)
        ax.set_ylabel(r'$g(r)$ + offset', fontsize=11)

        # Título con juicio estructural
        v_imp = np.sqrt(2 * G * h0)
        if h0 <= 11.5:
            estado = 'FCC conservado'
        elif h0 <= 30:
            estado = 'FCC deformado'
        else:
            estado = 'FCC $\\rightarrow$ amorfo'

        ax.set_title(
            f'$h_0 = {h0}\\,\\sigma$   ($v_0={v_imp:.2f}\\,\\sigma/\\tau$)\n'
            f'{estado}',
            fontsize=11, fontweight='bold', pad=8)

        ax.legend(fontsize=8, loc='lower right',
                  framealpha=0.90, edgecolor='#cccccc')

    fig.suptitle(
        r'Función de distribución de pares $g(r)$ de la bola — evolución temporal',
        fontsize=12, y=1.02)

    fig.tight_layout()
    fig.savefig(OUTFIG, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"\nFigura guardada: {OUTFIG}")


# ─── Main ─────────────────────────────────────────────────────────────────────

T_LABELS = [
    r'$t=0$  (FCC inicial)',
    r'$t_\mathrm{imp}$  (compresión máx.)',
    r'$t_\mathrm{apo_1}$  (1° apogeo)',
    r'$t_\mathrm{apo_2}$  (2° apogeo)',
]


def main():
    panels_data = []

    for h0 in H0_PANELS:
        dump_path = BASE / f"h0_{h0}" / "ball_dump.lammpstrj"
        if not dump_path.exists():
            print(f"[SKIP] Sin dump para h0={h0}  ({dump_path})")
            continue

        print(f"\n── h0 = {h0} σ ──────────────────────────────────────")
        frames = parse_ball_dump(dump_path)
        print(f"  {len(frames)} frames leídos  "
              f"(N_ball={len(next(iter(frames.values())))} átomos)")

        key_steps = get_key_steps(h0)
        curves    = []

        for step_target, label in zip(key_steps, T_LABELS):
            if step_target is None:
                continue
            step   = nearest_frame(frames, step_target)
            t_real = step * DT
            pos    = frames[step]
            r, gr  = compute_gr(pos)

            t_tag  = label + f'  $t={t_real:.1f}\\,\\tau$'
            curves.append((t_tag, r, gr))

            # Diagnóstico: altura del primer pico FCC
            idx1  = int(R_FCC[0] / DR)
            peak1 = gr[max(0, idx1-1): idx1+2].max() if len(gr) > idx1 else 0
            print(f"  {label.split('$')[1]:12s}  paso={step:6d}  t={t_real:6.1f}τ  "
                  f"N={len(pos)}  g(r₁)≈{peak1:.1f}")

        if curves:
            panels_data.append({'h0': h0, 'curves': curves})

    if panels_data:
        make_figure(panels_data)
    else:
        print("\nNo se encontró ningún ball_dump.lammpstrj.")
        print("Corre primero las simulaciones con el dump activado en in.bounce_reprod.")


if __name__ == '__main__':
    main()
