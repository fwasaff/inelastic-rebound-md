# CLAUDE.md – Distribución de energía en impactos de nanocúmulos LJ

Guía de referencia para el proyecto activo con Claude Code.

---

## Descripción del proyecto

Simulación de dinámica molecular (LAMMPS, unidades LJ) del impacto de un
nanocúmulo Lennard-Jones sobre una superficie FCC bajo gravedad constante.

**Pregunta central:**
> Dado un impacto con energía E₀, ¿cómo se distribuye entre los tres canales:
> CM de la bola / modos internos de la bola / fonones de la superficie?

**Observable clave:** La fracción α_int ≈ 22% que va a modos internos es
~constante para v > 0.9 σ/τ (con N=83). La hipótesis a verificar es:

```
α_int(N) = (3N - 3) / (3N - 3 + 3·N_c)
```

donde N_c es el número efectivo de átomos de superficie en la huella del impacto.

---

## Rama y remote

**Rama activa:** `claude/review-lammps-project-hra1A`
**Remote:** `https://github.com/fwasaff/inelastic-rebound-md`

```bash
git checkout claude/review-lammps-project-hra1A
git fetch origin && git status
```

---

## Estructura del repositorio

```
inelastic-rebound-md/
├── CLAUDE.md                           ← este archivo
├── lammps/
│   ├── in.bounce_reprod                ← script principal de simulación
│   └── clusters/                       ← geometrías de bola pre-minimizadas
│       ├── cluster_N1_min.dat          (1 átomo)
│       ├── cluster_N2_min.dat          (2 átomos)
│       ├── cluster_N7_min.dat          (7 átomos)
│       ├── cluster_N13_min.dat         (13 átomos)
│       └── cluster_N55_min.dat         (55 átomos)
├── scripts/
│   ├── gen_tabla_repulsiva.py          ← genera tabla LJ repulsiva (c=0)
│   ├── run_reprod.sh                   ← corre las 6 alturas del paper 2012
│   ├── analizar_reprod.py              ← extrae rebotes, alturas, COR
│   └── analizar_reprod_fisica.py       ← partición energética (4 figuras)
├── docs/
│   └── borrador_rebotes_gravedad.tex   ← borrador del artículo
└── data/
    └── outputs/
        └── reprod_2012/                ← resultados (ignorado por git)
```

---

## Hoja de ruta del estudio

### Paso 1 — Cerrar el balance energético (T_S en LAMMPS)
Agregar medición de energía de la superficie en `in.bounce_reprod`:
- `compute pe_surf_r surface reduce sum c_pe_atom`
- Sacar `c_surf_ke` y `c_pe_surf_r` en `bounce.dat`
- Verificar: E_ball + E_surf + U_grav = const rebote a rebote

### Paso 2 — Baseline N=1 (c_ps=0)
Correr las mismas 6 alturas con un átomo proyectil.
Esto mide α_surf(v) directamente (sin canal interno).

### Paso 3 — Variación sistemática N ∈ {2, 7, 13, 55} (c_ps=0)
Usar los clusters pre-minimizados en `lammps/clusters/`.
Comparar α_int(N) con la predicción analítica.

### Paso 4 — Observar h* (simulación larga)
Correr h₀=28σ, N=83, 5000τ para verificar convergencia al plateau
h* = T_B/(Mg) predicho por Gutiérrez 2012.

### Paso 5 — Adhesión c_ps=1
Repetir Pasos 2-3 con LJ completo. Medir umbral h₀*(N) donde la bola
no rebota (KE < energía de adhesión).

---

## Script LAMMPS principal: `in.bounce_reprod`

**Variables configurables vía `-var`:**

| Variable  | Defecto | Descripción |
|-----------|---------|-------------|
| h0        | 11.03   | Altura inicial del CM [sigma] |
| g         | 0.026   | Gravedad [sigma/tau²] |
| a0        | 1.5496  | Parámetro de red FCC argón LJ |
| Nx,Ny,Nz | 8,8,6   | Celdas de la superficie |
| rball     | 2.64    | Radio de la bola [sigma] |
| dt        | 0.005   | Paso de tiempo [tau] |
| t_total   | 400.0   | Tiempo de simulación [tau] |

**Nota crítica sobre `lattice fcc`:** En `units lj`, el parámetro de
`lattice fcc scale` es la densidad reducida ρ* = N/V, NO el parámetro de
red. Para a₀ = 1.5496σ se usa `$(4/(a0³)) ≈ 1.0742`.

**Observables en `bounce.dat`:**
```
time  Zcm  Vcmz  KB  TB  VB  Uc  UF  KB3d
```
- KB  = KE traslacional z del CM de la bola
- TB  = KE interna de la bola (= c_ball_ke - KB3d)
- VB  = PE intra-bola (= c_pe_ball - 0.5·c_Uc)
- Uc  = PE colisional bola-superficie (compute group/group)
- UF  = PE gravitatoria = M·g·Zcm
- KB3d= KE traslacional 3D del CM

---

## Cómo correr

### Reproducción Gutiérrez 2012 (N=83, c_ps=0)

```bash
# Generar tabla de potencial repulsivo
python3 scripts/gen_tabla_repulsiva.py \
    --output data/outputs/reprod_2012/table_repulsiva.txt

# Correr las 6 alturas
bash scripts/run_reprod.sh --jobs 4

# Analizar
python3 scripts/analizar_reprod.py --outdir data/outputs/reprod_2012
python3 scripts/analizar_reprod_fisica.py
```

### Una simulación individual

```bash
cd data/outputs/mi_run
cp ../../reprod_2012/table_repulsiva.txt .
lmp -in ../../../lammps/in.bounce_reprod \
    -var h0 28.15 -log lammps.log -screen none
```

---

## Hardware y tiempos de referencia

**CPU:** benchmark LAMMPS 7 Feb 2024 → ~1.37 M atom-steps/s

| Simulación          | Átomos | Pasos     | Tiempo estimado |
|---------------------|--------|-----------|-----------------|
| Un rebote (400τ)    | ~1620  | 80 000    | ~1.5 min        |
| 6 alturas en serie  | ~1620  | 480 000   | ~9 min          |
| Largo (5000τ)       | ~1620  | 1 000 000 | ~12 min         |

---

## Resultados de la reproducción (N=83, c_ps=0)

### Partición de energía en el primer rebote

| h₀_eff (σ) | v_imp (σ/τ) | e₁    | α_CM  | α_int | α_surf |
|------------|-------------|-------|-------|-------|--------|
| 8.53       | 0.666       | 0.906 | 0.822 | 0.068 | 0.110  |
| 16.21      | 0.918       | 0.697 | 0.486 | 0.211 | 0.303  |
| 25.77      | 1.158       | 0.662 | 0.438 | 0.178 | 0.384  |
| 38.94      | 1.423       | 0.660 | 0.435 | 0.223 | 0.342  |
| 56.41      | 1.713       | 0.618 | 0.382 | 0.219 | 0.399  |
| 77.24      | 2.004       | 0.503 | 0.253 | 0.233 | 0.514  |

**Hallazgo:** α_int ≈ 0.20–0.23 constante para v > 0.9 σ/τ.

### Nivel de contacto efectivo
z₀ ≈ 11.6–11.8 σ (media de los mínimos de Zcm). Todas las alturas h_n
se miden relativas a z₀. El parámetro LAMMPS h0 se mide desde z_box_top
(= Nz·a₀ ≈ 9.30σ); la diferencia es z₀ - z_box_top ≈ 2.3σ.

---

## Física de los scripts LAMMPS

### Potenciales
- Superficie-superficie: LJ completo, rc=2.5σ
- Bola-bola: LJ completo, rc=2.5σ
- Bola-superficie: puramente repulsivo (c=0), tabla pre-generada
  - Para c_ps=1 (Paso 5): `pair_coeff 1 2 lj/cut 1.0 1.0 2.5`

### Dinámica
- NVE global, sin termostato (reproduce Gutiérrez 2012)
- Capas inferiores fijas: `fix setforce 0 0 0` (2 capas del fondo)
- Gravedad solo en la bola: `fix gravity vector 0 0 -1`
- Pared reflectante en zlo: `fix wall/reflect` solo en `mobile_surf`

---

## Convenciones

- Unidades LJ: longitud en σ, energía en ε, tiempo en τ, masa en m
- Argón: ε ≈ 0.0104 eV, σ ≈ 3.4 Å, τ ≈ 2.15 ps, m = 39.95 uma
- Scripts de análisis: Python simple (numpy + matplotlib + scipy)
- No usar framework OOP; un script por tarea
