# Índice de `data/outputs/`

No versionado (`.gitignore: data/outputs/*/`), salvo este índice. Es la
carpeta de trabajo: todo aquí es **regenerable** corriendo el script
correspondiente sobre el `traj.lammpstrj`/`bounce.dat` de la corrida
indicada — nada de esto es la fuente de verdad.

## Estructura (2026-08-09: reorganizado de verdad, sin symlinks)

```
data/outputs/
├── activo/       ← línea onda/calor abierta (2026-08-08/09)
├── terminado/    ← canales de energía / N_c / sweep N, cerrado
├── retirado/     ← corrida inválida (bug de units) + lo que dependía de ella
└── INDEX.md
```

Sin archivos sueltos en la raíz. Cada corrida y cada figura vive físicamente
dentro de su categoría — no son enlaces. Las **124 referencias** que citaban
estas rutas por texto fijo (scripts con constantes hardcodeadas, `CLAUDE.md`,
la bitácora, y los `.tex` del manuscrito activo, de la reunión con Gonzalo y
de los borradores `historico/`) se actualizaron con el mismo prefijo de
categoría — verificado que cada una resuelve al archivo real tras el cambio.

## `activo/` — línea onda/calor

| Nombre | Tamaño | Qué es |
|---|---|---|
| `ondas_N79_libre/` | 3,8 GB | `boundary s s f`. La única corrida de ondas válida hoy. Toda la línea de 2026-08-08/09 sale de aquí. |
| `ondas_N79_grande/` | 925 MB | `boundary p p f`, **periódica**. La corrida donde se observó el rebote por imágenes periódicas que originó esta línea. |
| `viz_N79_h28/` | 55 MB | Figuras de reunión (5 scripts la usan). |
| `fig_polarizacion_frente.*`, `fig_eco_fondo.*`, `fig_pulso_contacto.png`, `fig_chi_coherencia.*`, `fig_chi_convergencia_celda.*`, `fig_correlacion_velocidad.*` | — | Objeciones 4/5/3/8, P4 (2026-08-09). Copia git-tracked en `docs/figuras_superficie/`. |
| `fig_reflexion_borde.*` | — | P5. Copia en `docs/figuras_superficie/`. |
| `fig_frente_onda.png`, `fig_coherencia_radial.png`, `fig_ondas_superficie*.png` | — | De `ondas_N79_libre`/`grande`/`viz_N79_h28`. p≈1.00–1.06 vigente. Copia en `docs/figuras_rbef/`. |
| `fig_anisotropia_eta.png` | — | Dato válido, pero η quedó **refutada como observable** el 2026-08-09 (objeción 4). Ver `fig_polarizacion_frente.png` para el reemplazo. |

## `terminado/` — canales de energía / N_c / sweep N (cerrado)

| Nombre | Qué es |
|---|---|
| `Nc_runs/`, `reprod_2012/`, `cluster_orient/` | N_c, sweep N, equipartición. |
| `cluster_N1/`…`cluster_N87/`, `cluster_N79_long/`, `adhesion_N13/`…`adhesion_N87/` | Geometrías/corridas del sweep N y adhesión. |
| `slab_test/` | Test de tamaño finito de losa (M4). |
| `fig_alpha_TB_*.pdf`, `fig_convergencia_hstar.pdf`, `fig_equiparticion_v2.pdf`, `fig_fusion_gr.pdf`, `fig_Nc_escalas.pdf`, `fig_paso_B_*.pdf`, `fig_paso_C_*.pdf`, `fig_virial_check.pdf`, `fig_sequence*.png` | Figuras del manuscrito RBEF. Citadas por `docs/manuscrito/historico/*.tex` — rutas actualizadas y verificadas. |

## `retirado/` — corrida inválida

| Nombre | Qué es |
|---|---|
| `ondas_N79_esponja_v2/` | `traj.lammpstrj` borrado (bug de `units` en la esponja, ver `LEEME_traj_borrado.md` dentro). Solo quedan `lammps.log` + `bounce.dat`, la evidencia del bug. |
| `fig_anisotropia_eta_esponja.png`, `fig_coherencia_radial_esponja.png`, `fig_frente_onda_esponja.png`, `fig_eta_esponja_ventana_limpia.png`, `fig_ondas_esponja_v2.png`, `fig_ondas_esponja_test.png` | Figuras de la corrida inválida. Se conservan a propósito (decisión 2026-08-09): son el rastro de lo que se creyó antes de refutarlo. `sesiones/2026-08-03.md` cita `fig_eta_esponja_ventana_limpia.png` por nombre. |

## Si generas algo nuevo

Los scripts que **escriben** con ruta fija (`run_slab_test.sh`, `run_orient_sweep.sh`,
`run_Nc_dumps.sh`, `run_reprod.sh`, `run_cluster_sweep.sh`, y las constantes
`BASE`/`RUN`/`DATA` de varios `.py`) ya apuntan a la ubicación con categoría
— una corrida nueva de esas cae directo en su carpeta. Para todo lo demás
(los scripts `medir_*.py`/`analizar_*.py` que reciben `--traj`/`--out` por
línea de comandos), pasa la ruta con categoría a mano, p. ej.:

```bash
python3 scripts/medir_frente_onda.py \
    --traj data/outputs/activo/ondas_N79_libre/traj.lammpstrj \
    --out  data/outputs/activo/fig_lo_que_sea.png
```

Si es un resultado de la línea onda/calor va en `activo/`; si es del sweep N
o RBEF, en `terminado/`; si resulta inválido, se mueve a `retirado/` (no se
borra sin decidirlo).
