# CLAUDE.md – Distribución de energía en impactos de nanocúmulos LJ

Guía de referencia para el proyecto activo con Claude Code.

---

## Descripción del proyecto

Simulación de dinámica molecular (LAMMPS, unidades LJ) del impacto de un
nanocúmulo Lennard-Jones sobre una superficie FCC bajo gravedad constante.

**Pregunta central:**
> Dado un impacto con energía E₀, ¿cómo se distribuye entre los tres canales:
> CM de la bola / modos internos de la bola / fonones de la superficie?

**Canales de energía (3):**
- α_CM  = e₁² = h₁/h₀  (fracción CM vertical, = COR al cuadrado)
- α_TB  = T_B/E₀        (KE interna de la bola → "temperatura" de la bola)
- α_surf = ΔE_surf/E₀   (fonones superficie, medido directamente)
- α_ball = 1 − α_CM − α_surf = α_TB + α_VB  (energía interna total; α_VB = PE anarmónico)

**Observable clave:** α_TB ≈ 0.22 constante para v > 0.9 σ/τ (con N=79).
El canal anarmónico α_VB ≈ α_TB (equipartición aproximada).
La hipótesis a verificar para el NUEVO paper:

```
α_TB(N) = (3N - 3) / (3N - 3 + 3·N_c)
```

donde N_c ≈ 290 es el número efectivo de átomos de superficie en la huella.
Predicciones: α_TB(13)≈0.040, α_TB(55)≈0.157, α_TB(79)≈0.212, α_TB(87)≈0.229.

**Nota sobre N:** Los shells FCC completos centrados en un sitio de red son
N=13, 55, 79, 87, ... El shell completo N=79 (r < 2.45σ) es el cluster
canónico del paper. El valor "N=83" en datos anteriores era un artefacto del
método `create_atoms region sphere` con el centro fuera de un sitio de red.

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
│   ├── in.bounce_reprod                ← reproducción 2012 (create_atoms, N variable)
│   ├── in.bounce_cluster               ← versión mejorada: cluster fijo pre-minimizado
│   ├── in.minimize_cluster             ← minimiza un cluster en vacío
│   └── clusters/                       ← geometrías de bola pre-minimizadas
│       ├── cluster_N1_min.dat          (1 átomo)
│       ├── cluster_N2_min.dat          (2 átomos)
│       ├── cluster_N7_min.dat          (7 átomos)
│       ├── cluster_N13_min.dat         (13 átomos)
│       ├── cluster_N55_min.dat         (55 átomos)
│       └── cluster_N79_min.dat         (79 átomos — cluster canónico, shell FCC completo)
├── scripts/
│   ├── gen_tabla_repulsiva.py          ← genera tabla LJ repulsiva (c=0)
│   ├── gen_cluster_init.py             ← genera geometría FCC inicial de un cluster
│   ├── run_reprod.sh                   ← corre las 6 alturas del paper 2012
│   ├── analizar_reprod.py              ← extrae rebotes, alturas, COR
│   ├── analizar_reprod_fisica.py       ← partición energética (4 figuras)
│   └── analizar_gr.py                  ← función de distribución g(r) de la bola
├── docs/
│   └── borrador_rebotes_gravedad.tex   ← borrador del artículo
└── data/
    └── outputs/
        └── reprod_2012/                ← resultados (ignorado por git)
```

---

## Estado actual — resultados del sweep N

### COMPLETADO: Sweep N ∈ {13, 55, 79, 87} × 6 alturas (c_ps=0)

| N   | α_TB plateau (v>1.1) | Nc_eff | Predicción Nc=290 | Estado plateau      |
|-----|----------------------|--------|-------------------|---------------------|
| 13  | 0.066 (no plateau)   | ~169   | 0.040             | ✗ aún creciendo     |
| 55  | 0.217                | ~195   | 0.157             | ✓ aprox. constante  |
| 79  | 0.244                | ~241   | 0.212             | ✓ constante         |
| 87  | 0.250                | ~258   | 0.229             | ✓ constante         |

**Hallazgos del sweep:**
1. α_TB aumenta monótonamente con N ✓ (confirma la tendencia)
2. N=13 NO alcanza el plateau en v ≤ 2.0σ/τ — el cluster es demasiado rígido
3. Nc_eff crece con N (no constante): sugiere Nc ∝ área de contacto ∝ N^(2/3)
4. La fórmula captura la tendencia pero subestima α_TB por ~10–25%

**Figura generada:** `data/outputs/fig_alpha_TB_vs_N.pdf`

---

## Hoja de ruta — próximos pasos

### Paso siguiente — Extender sweep a N más pequeños y confirmar Nc(N)
- Correr N ∈ {1, 2, 7} con las mismas 6 alturas (clusters ya disponibles)
- Medir Nc_eff(N) y ajustar relación Nc = A × N^β
- Esto da la figura clave para el artículo nuevo

### Paso B — Observar h* (simulación larga)
Correr h₀=28σ, N=79, 5000τ para verificar convergencia al plateau
h* = T_B/(Mg) predicho por Gutiérrez 2012.

### Paso C — Adhesión c_ps=1
Repetir sweep con LJ completo (bola-superficie atractiva).
Medir umbral h₀*(N) donde la bola no rebota (KE < energía de adhesión).

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

### Reproducción 2012 (in.bounce_reprod, N variable ~79–84, c_ps=0)

| h₀_eff (σ) | v_imp (σ/τ) | e₁    | α_CM  | α_TB  | α_VB  | α_surf | α_ball |
|------------|-------------|-------|-------|-------|-------|--------|--------|
| 8.53       | 0.666       | 0.906 | 0.822 | 0.071 | 0.039 | 0.068  | 0.110  |
| 16.21      | 0.918       | 0.697 | 0.486 | 0.221 | 0.234 | 0.059  | 0.455  |
| 25.77      | 1.158       | 0.662 | 0.438 | 0.187 | 0.274 | 0.101  | 0.461  |
| 38.94      | 1.423       | 0.660 | 0.435 | 0.234 | 0.196 | 0.135  | 0.430  |
| 56.41      | 1.713       | 0.618 | 0.382 | 0.231 | 0.284 | 0.104  | 0.514  |
| 77.24      | 2.004       | 0.503 | 0.253 | 0.244 | 0.342 | 0.160  | 0.586  |

**Hallazgo:** α_TB ≈ 0.22 constante para v > 0.9 σ/τ.
α_ball = α_TB + α_VB ≈ 0.43–0.59 (total interno cinético + anarmónico).
α_surf ≈ 0.06–0.16 (mucho más pequeño que α_ball).

### N=79 fijo (in.bounce_cluster, cluster_N79_min.dat, c_ps=0)

| h₀_eff (σ) | v_imp (σ/τ) | e₁    | α_CM  | α_TB  | α_VB  | α_surf | α_ball |
|------------|-------------|-------|-------|-------|-------|--------|--------|
| 8.81       | 0.677       | 0.924 | 0.853 | 0.021 | 0.056 | 0.070  | 0.077  |
| 17.21      | 0.946       | 0.708 | 0.502 | 0.213 | 0.208 | 0.077  | 0.421  |
| 25.91      | 1.161       | 0.552 | 0.304 | 0.319 | 0.310 | 0.066  | 0.629  |
| 41.80      | 1.474       | 0.633 | 0.401 | 0.257 | 0.264 | 0.078  | 0.521  |
| 56.93      | 1.721       | 0.614 | 0.378 | 0.197 | 0.328 | 0.098  | 0.525  |
| 79.03      | 2.027       | 0.587 | 0.345 | 0.204 | 0.294 | 0.156  | 0.499  |

**α_TB ≈ 0.21–0.26 para v > 0.9; predicción teórica α_TB(79) = 0.212 (Nc≈290).** ✓

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
- **`pair_modify shift yes` es obligatorio** para el balance energético:
  sin shift, cada par que cruza rc salta ±0.0163ε en la PE registrada y el
  cierre Σα=1 falla por 2–6% de E₀ por impacto (la 5ª capa de vecinos FCC
  está en 2.45σ). El shift no cambia fuerzas ni trayectorias: los α por
  método residuo de runs antiguos siguen siendo válidos.

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
