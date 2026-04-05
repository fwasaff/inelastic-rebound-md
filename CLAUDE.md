# CLAUDE.md – Inelastic Rebound MD

Guía de referencia para trabajar en este repositorio con Claude Code.

---

## Estado del repositorio y flujo de trabajo confirmado

**Rama principal de desarrollo:** `claude/review-lammps-project-hra1A`
**Remote:** `https://github.com/fwasaff/inelastic-rebound-md`

El repositorio local y GitHub están **sincronizados** (verificado 2026-04-05).

### Workflow estándar con Claude Code

```bash
# Siempre trabajar desde la rama de desarrollo
git checkout claude/review-lammps-project-hra1A

# Antes de empezar a editar: verificar que estamos sincronizados
git fetch origin && git status

# Después de cambios: commit + push
git add -p
git commit -m "descripción concisa"
git push
```

Claude Code opera directamente sobre este directorio. No clonar salvo que la rama local
esté desincronizada respecto a origin (verificar con `git fetch && git status`).

---

## Tiempos de simulación (medidos en este notebook)

**Hardware:** CPU único, benchmark LAMMPS 7 Feb 2024 → **~1.37 M atom-steps/s**

| Fase                    | Pasos    | Átomos | Tiempo estimado |
|-------------------------|----------|--------|-----------------|
| Equilibración (default) | 85 000   | 1 536  | ~1.5 min        |
| Un impacto (default)    | 90 000   | 1 537  | ~1.7 min        |
| Pipeline completo       | 175 000  | 1 537  | ~3.5 min        |
| Sweep 10×1×3 replicas   | 30 runs  | —      | ~13 min (4 cores) |
| Sweep 10×3×3 replicas   | 90 runs  | —      | ~38 min (4 cores) |

Para acelerar: aumentar `--jobs` hasta el número de cores disponibles.

---

## Visualización: OVITO y ParaView

### Archivos generados para visualización

| Script          | Archivo dump                | Coordenadas | Campos extra        |
|-----------------|-----------------------------|-------------|---------------------|
| in.equilibrate  | `dump.equil.lammpstrj`      | **scaled** (xs ys zs) | pe/atom |
| in.impact       | `dump.impact.lammpstrj`     | **real** (x y z)      | pe/atom |

Ambos archivos están en `data/outputs/<fase>/` tras correr el runner Python.

### OVITO (recomendado)

**Instalación:**
```bash
pip install ovito          # OVITO Python API (headless + scripting)
# o descargar OVITO Pro/Basic desde https://www.ovito.org/
```

**Abrir dump en OVITO GUI:**
1. File → Open → seleccionar `dump.impact.lammpstrj`
2. OVITO detecta automáticamente el formato LAMMPS custom
3. Para `dump.equil.lammpstrj`: las coordenadas son **scaled** — OVITO las convierte solas
4. Colorear por `pe` (energía potencial por átomo) para ver la perturbación del impacto

**Nota:** OVITO no está instalado en este sistema. ParaView sí está disponible (`/usr/bin/paraview`).

### ParaView (disponible en este sistema)

**Abrir dump:**
1. File → Open → seleccionar `dump.impact.lammpstrj`
2. En el diálogo de formato: seleccionar **"LAMMPS Dump"** (no CSV)
3. Apply → usar "Point Gaussian" representation para ver los átomos

**Colorear por energía potencial:**
- En el panel izquierdo: cambiar "Solid Color" → `c_pe_atom`
- Ajustar escala de colores con el botón de rescale

**Animación del impacto:**
- Play en la barra inferior para animar frames del dump
- `dump_freq = 100` → 1 frame cada 100 pasos = 900 frames totales para in.impact

**Filtros útiles en ParaView:**
- `Threshold` sobre `c_pe_atom` → aislar átomos perturbados
- `Glyph` → representar átomos como esferas con radio = 0.5 sigma

---

## Descripción del proyecto

Simulación de **dinámica molecular de rebote inelástico** a escala nanométrica.

**Física simulada:**
Un átomo proyectil Lennard-Jones impacta sobre un cristal FCC Argón
y rebota inelásticamente: parte de la energía cinética normal se transfiere
a fonones de la superficie (calor). El observable principal es el
**Coeficiente de Restitución (COR)**:

```
e = |v_f,z| / |v_0,z|   (0 ≤ e ≤ 1)
```

**Unidades:** LJ reducidas (eps, sigma, tau).  
Para Argón: eps ≈ 0.0104 eV, sigma ≈ 3.4 Å, tau ≈ 2.15 ps.

---

## Estructura del repositorio

```
inelastic-rebound-md/
├── CLAUDE.md                  ← este archivo
├── lammps/
│   ├── in.equilibrate         ← Fase 1: equilibrar superficie FCC
│   └── in.impact              ← Fase 2: lanzar proyectil y registrar rebote
├── python/
│   ├── __init__.py
│   ├── exceptions.py          ← jerarquía de excepciones del proyecto
│   ├── config.py              ← metaclase ConfigMeta + dataclasses de configuración
│   ├── simulation.py          ← wrappers OOP sobre subprocess LAMMPS
│   ├── analysis.py            ← análisis físico con decoradores y caché
│   ├── visualization.py       ← figuras matplotlib (context manager de estilo)
│   └── runner.py              ← CLI: equilibrate / impact / sweep / analyze / plot
├── scripts/
│   ├── run.sh                 ← pipeline completo en bash
│   └── analyze.sh             ← análisis post-simulación
├── tests/
│   ├── __init__.py
│   └── test_analysis.py       ← tests unitarios con datos sintéticos (sin LAMMPS)
├── data/
│   └── outputs/               ← salida de simulaciones (ignorado por git)
├── docs/
│   └── physics_notes.md       ← notas de física (opcional)
└── requirements.txt
```

---

## Dependencias

### LAMMPS
```
LAMMPS >= 23 Jun 2022 (stable)
Estilos requeridos: MANYBODY, GRANULAR, EXTRA-PAIR (o compilación estándar)
Variables de entorno:
  LAMMPS_EXE=lmp              # nombre del ejecutable (por defecto: lmp)
  LAMMPS_EXE=lmp_mpi          # para versión MPI
```

Verificar:
```bash
lmp -version
```

### Python
```
Python >= 3.10  (se usan union types: Path | str, match/case futuro)
numpy >= 1.24
matplotlib >= 3.7
pytest >= 7.0
```

Instalar:
```bash
pip install -r requirements.txt
```

---

## Cómo ejecutar

### Flujo completo (recomendado)

```bash
# 1. Equilibrar la superficie FCC
python -m python.runner equilibrate --T 0.1 --Nx 8 --Ny 8 --Nz 6

# 2. Simular un impacto
python -m python.runner impact --v0 2.0 --angle 0.0

# 3. Analizar resultados
python -m python.runner analyze \
    --data data/outputs/v0_2.00_ang_0.0/rebound_data.dat \
    --v0 2.0

# 4. Generar figuras
python -m python.runner plot trajectory \
    --data data/outputs/v0_2.00_ang_0.0/rebound_data.dat \
    --output-dir data/outputs/figures
```

### Script shell todo-en-uno

```bash
bash scripts/run.sh
```

### Barrido paramétrico

```bash
python -m python.runner sweep \
    --v0-min 0.5 --v0-max 5.0 --v0-steps 10 \
    --angles 0.0 30.0 60.0 \
    --replicas 3 --jobs 4
```

### Sólo LAMMPS (sin Python)

```bash
cd data/outputs/equilibration
lmp -in ../../../lammps/in.equilibrate
lmp -in ../../../lammps/in.impact -var v0 2.0 -var angle 0.0
```

---

## Variables LAMMPS sobreescribibles desde CLI

Ambos scripts aceptan `-var nombre valor` desde la línea de comandos.

| Script           | Variable         | Defecto | Unidades   |
|------------------|------------------|---------|------------|
| in.equilibrate   | T_target         | 0.1     | eps/kB     |
| in.equilibrate   | Nx, Ny, Nz       | 8,8,6   | celdas     |
| in.equilibrate   | steps_equil      | 50000   | pasos      |
| in.equilibrate   | dt               | 0.005   | tau        |
| in.impact        | v0               | 2.0     | sigma/tau  |
| in.impact        | angle            | 0.0     | grados     |
| in.impact        | T_target         | 0.1     | eps/kB     |
| in.impact        | steps_impact     | 80000   | pasos      |

---

## Tests

Los tests NO requieren LAMMPS (usan datos sintéticos).

```bash
pytest tests/ -v
```

Cobertura esperada: config, parsers, análisis físico, excepciones, análisis paramétrico.

---

## Arquitectura Python – decisiones de diseño

### Jerarquía de excepciones (`exceptions.py`)
```
ReboundMDError
├── IOError
│   ├── FileNotFoundError
│   ├── ParseError
│   └── RestartNotFoundError
├── SimulationError
│   ├── LAMMPSNotFoundError
│   ├── SimulationRuntimeError
│   └── TimeoutError
├── ConfigError
│   ├── ValidationError
│   └── MissingParameterError
└── PhysicsError
    ├── InsufficientDataError
    ├── NonPhysicalResultError
    └── ConvergenceError
```

### Metaprogramación en `config.py`
- `FieldDescriptor`: descriptor de datos con coerción de tipo y validación.
- `ConfigMeta` (metaclase): genera automáticamente `__init__`, `__repr__`, `__eq__`, `to_dict`, `from_dict`, `to_json`, `from_json`.
- `ConfigRegistry` + `__init_subclass__`: registro automático de subclases para deserialización JSON por tipo.

### Decoradores en `analysis.py`
- `@timed`: mide tiempo de ejecución de métodos.
- `@require_data`: guarda para métodos que necesitan datos cargados.
- `@cached`: caché por instancia, invalidada al recargar datos.
- `@validate_physical(qty, lo, hi)`: fábrica de decoradores paramétrica; verifica bounds físicos.

### Gestores de contexto
- `EquilibrationSimulation` / `ImpactSimulation`: crean directorio de trabajo, cleanup en `__exit__`.
- `publication_style()` en `visualization.py`: cambia `rcParams` temporalmente.

---

## Lógica física de los scripts LAMMPS

### `in.equilibrate` – evaluación de la lógica

**Correcto:**
- Lattice FCC con `a0 = 1.5496 sigma` (valor de literatura para LJ Ar a T*=0.1).
- Capa inferior fija (`setforce 0 0 0`) → absorbe momento, simula bulk semi-infinito.
- Capa intermedia con termostato Nosé-Hoover NVT → baño de calor.
- Capa superior NVE libre → superficie limpia para impacto.
- Cutoff LJ 2.5 sigma con corrección de cola (`pair_modify tail yes`).
- Rampa de temperatura al inicio para evitar inestabilidades de energía.

**Punto a vigilar:** La corrección de cola (`tail yes`) asume bulk homogéneo.
Para un slab con vacío, introduce un error pequeño pero sistemático en la presión.
Para impacto (energía relativa), el efecto es despreciable.

### `in.impact` – evaluación de la lógica

**Correcto:**
- Lee restart → la superficie ya está equilibrada.
- Proyectil tipo 2, mismo potencial LJ (homoatómico: Ar-Ar).
- Termostat NVT sólo en capas internas, proyectil en NVE → no hay sesgo artificial.
- Registro de alta frecuencia del proyectil con `fix print` (cada 10 pasos).
- Variables `vz0` y `vx0` calculadas correctamente con cos/sin del ángulo.

**Punto a vigilar:** `z_surf_top = 9.3` es una estimación (Nz=6, a0=1.5496).
El valor exacto depende de la relajación térmica. Para análisis preciso,
calcular como `max(z[type==1])` en el post-procesado Python.

---

## ¿Por qué Python avanzado + LAMMPS es una buena idea?

| Necesidad de investigación       | Herramienta Python usada            |
|----------------------------------|-------------------------------------|
| Validar parámetros de simulación | `ConfigMeta` + `FieldDescriptor`    |
| No perder resultados por errores | Jerarquía de excepciones tipada     |
| Reproducir experimentos          | `to_json` / `from_json` en configs  |
| Barrido paramétrico masivo       | `ProcessPoolExecutor` en `runner.py`|
| Evitar recalcular análisis       | `@cached` en `ReboundAnalyzer`      |
| Figuras reproducibles            | `publication_style()` context mgr   |
| Tests sin LAMMPS                 | Datos sintéticos + pytest           |

---

## Convenciones de código

- **PEP 8** estricto. Líneas ≤ 99 caracteres.
- **Type hints** en todas las firmas públicas (Python 3.10+).
- **Docstrings** en formato NumPy para clases y funciones públicas.
- **No** usar `print()` en módulos library; usar `logging`.
- **No** modificar scripts LAMMPS directamente en producción; usar `-var` overrides.
- Nuevos tests deben correr sin LAMMPS instalado.

---

## Flujo de trabajo Git

```bash
git checkout claude/review-lammps-project-hra1A   # rama de desarrollo
# ... editar ...
git add -p                                          # staging selectivo
git commit -m "descripción concisa del cambio"
git push -u origin claude/review-lammps-project-hra1A
```

---

## Extensiones futuras sugeridas

1. **Potenciales más realistas**: EAM (metales) usando `pair_style eam/fs`.
2. **Cluster proyectil**: crear con `create_atoms region` para estudiar tamaño de impacto.
3. **OVITO**: integrar `ovito` Python API para análisis estructural (CNA, DXA).
4. **Jupyter notebooks**: en `docs/notebooks/` para resultados interactivos.
5. **Comparación con teoría de Hertz**: superponer `delta ~ F^(2/3)` en curvas de fuerza.
