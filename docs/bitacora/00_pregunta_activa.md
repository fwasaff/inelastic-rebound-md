# Pregunta activa

> Este archivo tiene UNA pregunta. Si hay dos, una de las dos es un desvío.
> Se actualiza cuando la pregunta se responde o cuando Gonzalo cambia el rumbo.

**Última actualización:** 2026-08-09 (P1–P5 ejecutados. η cerrado por medición.
χ medido y **no converge**. La restricción que manda resultó ser la PROFUNDIDAD
de la losa, no el tamaño lateral)

---

## La pregunta

**¿La energía que entra a la superficie es sonido coherente o calor, y dónde
está la frontera entre ambos?**

Mandato de Gonzalo tras la reunión de julio 2026: *estudiar la superficie, la
onda de sonido y el calor*. Su hipótesis explícita: **frente coherente, estela
térmica**.

**H1 (Gonzalo).** La energía entregada a la superficie se separa en un frente
balístico que conserva coherencia mecánica y una estela que termaliza detrás,
con una constante de tiempo τ_th del orden del tiempo de dispersión fonónico.
La firma que la hace falsable: **τ_th se acorta al aumentar la velocidad de
impacto**, porque la anarmonicidad escala con la amplitud.

**H0 (la nula que hay que batir).** No hay termalización en la ventana: todo es
transitorio elástico lineal + campo difuso tipo Weaver. τ_th es **independiente
de la velocidad de impacto**.

**El discriminador es τ_th(v), no el valor de τ_th en una corrida.**

### Forma operativa — EN DISPUTA

Han caído dos candidatos en dos semanas:

- **η** (anisotropía radial). Refutado por medición el 2026-08-09: mide
  polarización, no coherencia. Correlaciona a 0.9989 con la predicción de una
  onda 100% coherente que tenga la ρ_a medida. Cerrado.
- **χ = E_coh/(E_coh+E_term)** (descomposición de MD de choque). Medido el
  2026-08-09: **no converge**. Varía un factor 4.4 según el tamaño de celda,
  monótonamente y sin plateau. No es un número.

- **C(d) = ⟨v̂_i·v̂_j⟩** (correlación de dirección). El decaimiento
  τ ≈ 0.25 τ que dio el 2026-08-09 **está retirado el mismo día**: era ceguera a
  la amplitud. Ver objeción 8, cerrada abajo.

Candidato vivo: **C_amp(d) = Σ(v_i·v_j)/Σ(|v_i||v_j|)**, que pesa por amplitud y
sigue acotada en [−1,1] (Cauchy-Schwarz término a término). Es el único
estimador que ha sobrevivido a su propia objeción. Con él, en esta corrida:
**ξ = 2.12 σ** (~2 vecinos), C ≈ 0 más allá de 4 σ, y plano en el tiempo.

## Estado — 2026-08-09

**La pregunta NO está contestada. Lo que se sabe hoy es POR QUÉ no se puede
contestar con las corridas existentes.**

Tres tiempos medidos que chocan (todos de `ondas_N79_libre`):

| | valor | medido en |
|---|---|---|
| duración del forzamiento (contacto bola-superficie) | **3.20 τ** | `medir_pulso_contacto.py` |
| vuelta del eco del fondo rígido | **1.80 τ** | `medir_eco_fondo.py` |
| salida del frente de la zona de contacto (~10 σ) | **~1 τ** | `medir_polarizacion_frente.py` |

**La ventana limpia (1.8 τ) es más corta que el forzamiento (3.2 τ).** No se
mide un tiempo de termalización en una ventana más corta que la inyección. Y
dentro de ella el frente no se ha separado de la fuente: frente y estela son la
misma región.

**La restricción que manda es la profundidad de la losa (4.65 σ), no el ancho de
la caja (124 σ).** La superficie se agrandó en la dirección equivocada: el eco
lateral no vuelve hasta ~13 τ; el vertical vuelve a 1.8 τ.

### Números firmes de esta corrida

- velocidad del frente **4.76 σ/τ** (no 8.28; ese es el longitudinal de volumen
  de `equiparticion_y_virial.tex`)
- ρ_a = rms(v_z)/rms(v_ρ) = **1.71 ± 0.62** en el frente
- las dos anteriores son firmas independientes de **onda de Rayleigh**
  (c_R ≈ 0.5–0.55 c_L, ρ_a ≈ 1.5 para ν ≈ 0.25–0.3)
- τ_c (FWHM del pulso) = **2.10 τ**; λ_L = c_L τ_c = **17.4 σ**
- ξ (longitud de correlación de velocidad) = **3.38 σ** = 3.1 vecinos
- R_lateral ≈ **1**, R_fondo = **1**: ambos bordes reflejan todo
- ω_Debye = **33.1 τ⁻¹**; el dump a 0.25 τ está aliaseado por 2.6×

### Corregido en el registro

- **Cuáles corridas son periódicas y cuáles no** (el 2026-08-08 esto quedó mal
  anotado como "ninguna es periódica", y no es así):
  - `ondas_N79_grande` (jul 21, 18:12): **`boundary p p f`, SÍ periódica.** Es la
    corrida donde se observó el rebote por las imágenes periódicas, que es lo
    que motivó toda esta línea. La observación original era correcta.
  - `ondas_N79_libre` (jul 21, 19:18) y `ondas_N79_esponja_v2`: `boundary s s f`,
    header del dump `ss ss ff`. Aquí el eco lateral es de **borde libre**, no
    wraparound — y refleja igual de bien (R≈1, medido).
  O sea: se pasó de periódica a libre para matar el eco, y no lo mató, porque un
  borde libre también refleja del todo. Lo que sí cambió es *cuándo* vuelve.
- `ondas_N79_esponja_v2` sigue retirada (bug de `units`, ver 2026-08-08). El bug
  **ya está corregido** en `lammps/in.bounce_cluster` y el input lleva ahora un
  bloque de verificación probado contra la versión con bug.

## Objeciones

**Objeción 8 — CERRADA el 2026-08-09, y se llevó por delante el resultado.**
Con el estimador pesado por amplitud, el decaimiento τ ≈ 0.25 τ **no existe**:
es ruido en torno a cero desde el primer frame. Comparación en el mismo dato,
⟨C⟩ en 8–15 σ:

| t − t_imp | 0.05 | 0.30 | 0.55 | 0.80 | 1.05 | 1.30 |
|---|---|---|---|---|---|---|
| dirección | +0.500 | +0.390 | +0.188 | +0.061 | +0.022 | −0.007 |
| **amplitud** | −0.004 | −0.080 | +0.102 | +0.010 | +0.059 | +0.013 |

Lo que decaía eran los átomos **quietos** del campo lejano: sus velocidades
minúsculas están débilmente alineadas por el campo elástico cuasiestático de la
bola apoyándose, y el estimador de dirección los pesa igual que a los del
frente. Exactamente el defecto que la objeción 3 le achacaba a η, reintroducido
sin querer y detectado sólo al construir el estimador alternativo.

**Lo que queda medido, con el estimador bueno:** ξ = 2.12 σ, C ≈ 0 para d > 4 σ,
sin dependencia temporal. Dentro de la ventana accesible, la energía que entra a
la superficie **no tiene estructura coherente más allá de ~2 distancias
interatómicas**. Eso suena a "todo calor", pero **no es una respuesta a la
pregunta de Gonzalo**: la ventana no es válida (forzamiento activo, eco del
fondo a 1.8 τ, frente sin separarse de la fuente). Es una medición del campo
cercano de una fuente todavía encendida.

**Objeción 7 — ABIERTA. χ no tiene plateau.** Factor 4.4 entre ℓ = 1.55 y
6.20 σ. Sin una ℓ con n ≫ 1 y a la vez ℓ ≪ λ, χ no es reportable. Agravante: λ
no está medida — C(d) no resuelve ninguna onda periódica.

**Objeción 9 — ABIERTA, estructural.** La ventana limpia es más corta que el
forzamiento. No se arregla con análisis; requiere losa profunda con esponja de
fondo.

**Objeciones 1, 4, 5, 6 — CERRADAS el 2026-08-09.** Ver `sesiones/2026-08-09.md`.
**Objeción 2 — CERRADA el 2026-08-03.**
**Objeción 3 — CERRADA por obsolescencia** (η retirado).

## Lo siguiente

**La corrida nueva.** `ondas_N79_libre` está exprimida: con un observable que
sobrevive a sus objeciones (C_amp) y una ventana que no sirve, no queda nada más
que sacarle sin cómputo nuevo. Las tres objeciones que bloqueaban el análisis
(4, 5, 8) están cerradas; la que queda es la 9, que es de geometría.

**Y ahora sí se sabe qué volcar**, que era la decisión irreversible que faltaba:
velocidades por átomo de las **capas superiores**, en un radio que **siga al
frente** — no en un anillo fijo. Las dos lecciones de esta ronda (P1: el anillo
va sobre el frente medido, no sobre c·t; objeción 8: hay que pesar por amplitud)
apuntan a lo mismo.

**Predicción que la corrida nueva tiene que discriminar:** si en una ventana
limpia ξ_amp sigue valiendo ~2 σ en el frente, no hay frente coherente y H1 cae.
Si ξ_amp crece hasta ~λ/2 en el frente y decae hacia ~2 σ en la estela, es H1, y
el ritmo de esa caída es τ_th. El discriminador final sigue siendo **τ_th(v)**,
o sea tres velocidades.

### La corrida nueva, con las cotas ya medidas

- **Profundidad, que es lo que manda.** ~15 σ limpios + **esponja de fondo de
  ≥ 30 σ** (criterio: espesor ≥ λ_L = 17.4 σ). Total ≈ 45 σ, contra los 4.65 σ
  actuales.
- **Lateral.** Para 10 τ de ventana limpia hace falta r_limpio ≥ c_R·T = 46 σ,
  más 40 σ de esponja → semiancho 86 σ, **L ≈ 172 σ**.
- **γ_max ≈ 1–3, NO 8.** Con rampa cuadrática la profundidad óptica de ida y
  vuelta es 2γ_max W/(3c): con γ_max=3, W=40 σ da e⁻⁹·⁷ en amplitud, y hasta
  γ_max=1 da R_E ~ 10⁻³. Conviene equivocarse por abajo: γ alto no absorbe más
  y sí empina el gradiente de impedancia.
- **≥ 15 sub-anillos**, no 5: cada escalón discreto es un salto de impedancia.
- **`units box` explícito en toda región** (ya corregido, con guardia en el log).
- **dt_dump < 0.095 τ** (Nyquist sobre ω_Debye = 33.1). Usar 0.025 τ.
- **Coste de cómputo: ~20 min** (1.4 M átomos, ~4000 pasos, a los 4.74 M
  atom-steps/s medidos en el log de esponja_v2). No es el problema.
- **El problema es el dump: 48 GB** si se vuelca todo cada 0.025 τ. Hay que
  volcar sólo un grupo (capas superiores + r < 60 σ) → ~4 GB. Decidir el grupo
  ANTES de lanzar.
- **Criterio de aceptación: R_lateral < 0.01 medido** con
  `scripts/medir_reflexion_borde.py`, no supuesto.
- Soltar la bola desde cerca: los 46.7 τ de caída libre son 2.1 GB de dump inútil.

## Lo que NO es la pregunta activa

- El manuscrito RBEF (`borrador_rbef_v2.tex`): cerrado salvo la ronda 3.
- `onda_vs_calor_fonones.tex`: **no tocar**. La línea 440 dice "η no medido
  todavía"; hoy sabemos que η además no sirve. Reescribir la Sec. 6 es trabajo
  para cuando haya un observable que sobreviva.
- Canales de energía, N_c, f_ball, equipartición, sweep en N: terminado.
- Material de reunión, videos, comparaciones con 2012: producto.
- Fronteras GLE/PML exactas, SED completo, TTM: ver `ideas_aparcadas.md`.
