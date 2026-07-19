# inelastic-rebound-md

Simulación de dinámica molecular (LAMMPS, unidades LJ) del impacto de un
nanocúmulo Lennard-Jones sobre una superficie FCC bajo gravedad constante:
¿a dónde se va la energía de una pelota que rebota?

Manuscrito: `docs/manuscrito/borrador_rbef_v2.tex` (Revista Brasileira de
Ensino de Física, sección Artigos Gerais).

## Requisitos

- **LAMMPS** (probado con la versión del 7 feb 2024). En Ubuntu/Debian:
  ```bash
  sudo apt install lammps
  ```
  Esto instala el binario como `lmp` (o `lmp_serial`/`lmp_mpi` según la
  distribución). Si tu sistema no lo empaqueta, compílalo desde
  <https://docs.lammps.org/Install.html> (basta el paquete estándar, sin
  paquetes opcionales).
- **Python 3** con `numpy`, `matplotlib` y `scipy`:
  ```bash
  pip install numpy matplotlib scipy
  ```
- Opcional, para visualizar trayectorias: [OVITO](https://www.ovito.org/) o
  [ParaView](https://www.paraview.org/).

## Primer archivo a mirar

`lammps/in.bounce_cluster` — el script principal de simulación (~150 líneas
comentadas). Corre un único rebote de un cúmulo pre-minimizado
(`lammps/clusters/cluster_N79_min.dat`, el cúmulo canónico de 79 átomos)
contra una superficie FCC.

## Comando mínimo para correr un caso

```bash
mkdir -p data/outputs/mi_run && cd data/outputs/mi_run
python3 ../../../scripts/gen_tabla_repulsiva.py --output table_repulsiva.txt
lmp -in ../../../lammps/in.bounce_cluster \
    -var h0 28.15 -var cluster_file ../../../lammps/clusters/cluster_N79_min.dat \
    -log lammps.log -screen none
```

Esto produce `bounce.dat` con la serie temporal
`time Zcm Vcmz KB TB VB Uc UF KB3d` (ver definiciones en el manuscrito,
Sección 2). Un rebote de 400τ toma ~1.5 minutos en un núcleo de portátil.

## Reproducir las figuras del artículo

```bash
bash scripts/run_orient_sweep.sh      # barrido N x alturas x orientaciones (210 corridas, ~6 h CPU)
python3 scripts/medir_Nc.py           # N_c geométrico y de participación
python3 scripts/medir_fusion.py       # T_ef y g(r) para delimitar fusión
python3 scripts/fig_Nc_escalas.py     # Figura 4 del manuscrito
```

## Estructura del repositorio

```
inelastic-rebound-md/
├── lammps/            # scripts de simulación (.in) y cúmulos pre-minimizados
├── scripts/           # generación de geometrías y análisis (Python)
├── docs/              # manuscrito/, arbitraje/, reunion_gonzalo/, figuras_rbef/
└── data/outputs/      # resultados de simulación (no versionado)
```

Ver `CLAUDE.md` para el estado detallado del proyecto y la hoja de ruta.
