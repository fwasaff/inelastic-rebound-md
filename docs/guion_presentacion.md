# Guion de estudio — presentación a Gonzalo

Este documento tiene tres partes: **(A)** los conceptos que tienes que
dominar antes de la reunión, explicados desde cero; **(B)** el guion
lámina por lámina (qué decir, qué señalar en cada figura, cuánto tiempo);
**(C)** las preguntas que Gonzalo probablemente hará —es profesor de
mecánica estadística— con las respuestas preparadas.

---

## PARTE A — Conceptos que debes dominar

### A1. Unidades de Lennard-Jones (por si pregunta por números)

Todo el trabajo está en unidades reducidas: la longitud se mide en σ (el
"tamaño" de un átomo, ~3.4 Å para argón), la energía en ε (la profundidad
del pozo atractivo entre dos átomos, ~0.0104 eV) y el tiempo en
τ = √(mσ²/ε) ≈ 2.15 ps. Cuando digas "la bola cae desde 26σ" estás
diciendo "desde ~9 nanómetros"; cuando digas "el impacto dura ~2τ" es
"~4 picosegundos". La gravedad g* = 0.026 σ/τ² es **enormemente** más
grande que la gravedad real a esa escala — es una gravedad artificial
para que el experimento quepa en tiempos de simulación. Si Gonzalo
pregunta: no afecta la física del impacto, porque lo que importa es la
velocidad de impacto v_imp = √(2g·h), y esa la controlamos con h₀.

### A2. El coeficiente de restitución (COR) y por qué e < 1

- Definición: e = |v después| / |v antes| del impacto. Como la altura de
  vuelo es h = v²/2g, se cumple **h₁/h₀ = e²**.
- En los cursos de mecánica, e se trata como constante del material. El
  primer resultado del trabajo es que **e no es constante**: va de 0.92
  (impacto suave) a 0.55–0.59 (impacto fuerte). Una pelota atómica rebota
  "mejor" cuando cae despacio.
- Intuición: a baja velocidad la bola no alcanza a excitar sus
  vibraciones internas (el impacto es "adiabático" respecto de los modos
  internos); a alta velocidad el golpe es violento comparado con el
  período de vibración de los átomos y la energía se reparte.

### A3. Los cuatro canales — la ecuación más importante del paper

En el **primer apogeo** (el punto más alto después del primer rebote,
donde la velocidad del centro de masa es cero) se mide:

    α_cm + α_TB + α_ΔVB + α_surf = 1        (todo dividido por E₀)

1. **α_cm = h₁/h₀ᵉᶠ = e₁²** — la fracción que el centro de masa (CM)
   recuperó. Es la única energía "útil" para seguir rebotando.
2. **α_TB** — energía cinética de los átomos de la bola MEDIDA EN EL
   SISTEMA DEL CM (es decir, restando el movimiento de traslación).
   Es el "calor" cinético: los átomos vibran alrededor de sus posiciones.
   La T es de "temperatura-like": T_B es literalmente lo que usarías para
   definir una temperatura interna instantánea de la bola
   (T_B = (3N−3)·k_B·T_interna/2 si estuviera equilibrada).
3. **α_ΔVB** — el cambio de energía POTENCIAL interna de la bola entre
   antes del impacto y el apogeo. Ver A4: es la "expansión anarmónica".
4. **α_surf** — cuánta energía (cinética + potencial) ganó la superficie.
   Son los fonones: ondas elásticas que el golpe radia por el cristal.

**El mensaje que tienes que poder decir en una frase:** "la energía que
la pelota no recupera queda mayoritariamente DENTRO de la pelota
(α_TB + α_ΔVB ≈ 0.4–0.6), y el piso se lleva poco (0.07–0.16). La
intuición macroscópica de que 'la energía se disipa en el suelo' falla a
esta escala."

### A4. La expansión anarmónica (α_ΔVB) — el canal "nuevo"

Este es el canal conceptualmente más fino, prepáralo bien:

- El potencial LJ **no es simétrico** alrededor del mínimo: es muy duro
  al comprimir (r⁻¹²) y blando al estirar (r⁻⁶). Eso es "anarmonicidad".
- Cuando la bola se calienta internamente (T_B sube), sus átomos oscilan
  en ese pozo asimétrico y pasan **más tiempo del lado estirado**: la
  distancia media entre átomos crece. Es EXACTAMENTE el mismo mecanismo
  de la **dilatación térmica** de los sólidos — a Gonzalo le va a gustar
  esa conexión.
- Consecuencia: en el apogeo, la energía potencial interna V_B es MAYOR
  que antes del impacto (la bola está "expandida"), y esa diferencia
  ΔV_B viaja con la bola. Es energía potencial almacenada, tan real como
  la cinética.
- Números: α_ΔVB ≈ 0.21–0.33, comparable a α_TB ≈ 0.22. De hecho, que
  sean comparables ES equipartición: en un oscilador, ⟨KE⟩ ≈ ⟨PE⟩
  (virial). Si te preguntan por qué son parecidos, esa es la respuesta:
  "es el teorema del virial actuando dentro de la bola".

### A5. ¿Qué es el "plateau"? (tu duda — apréndete esto)

"Plateau" = meseta. Se refiere a esto: si graficas α_TB (la fracción
cinética interna) contra la velocidad de impacto v_imp, la curva **sube y
luego se aplana** en un valor casi constante:

- Para v_imp < 0.9 σ/τ: α_TB es pequeño y crece con v (régimen
  cuasi-elástico: el impacto apenas excita los modos internos).
- Para v_imp > 0.9–1.1 σ/τ: α_TB se queda "pegado" en ≈ 0.22 sin
  importar cuánto más fuerte golpees. ESO es el plateau: la meseta de la
  curva α_TB(v).

**Por qué existe la meseta (la física):** cuando el golpe es
suficientemente violento, el tiempo de contacto (~2τ) basta para que la
energía se reparta ("termalice") entre todos los grados de libertad
dinámicamente accesibles. Una vez que el reparto es completo, la FRACCIÓN
que le toca a cada canal ya no depende de cuánta energía entra — solo de
CUÁNTOS grados de libertad compiten. Por eso la fracción se vuelve
constante. Es el mismo motivo por el que el calor específico de un sólido
clásico no depende de la temperatura (ley de Dulong-Petit): la
equipartición reparte por igual entre modos, y la fracción de cada modo
es fija.

En cambio, si el golpe es suave (v pequeña), la termalización NO alcanza
a ocurrir y la fracción sí depende de v: por eso la curva sube antes de
aplanarse. Y en el cúmulo de N=13, que es muy rígido (una sola capa
cerrada, frecuencias propias altas), la termalización no ocurre ni
siquiera a v=2: **N=13 nunca alcanza el plateau** en el rango estudiado.
Por eso en la figura va con símbolo abierto y excluido del ajuste.

### A6. La fórmula de equipartición y N_c

    α_int(N) = (3N−3) / (3N−3 + 3N_c)

- **De dónde sale:** si en el contacto la energía E₀ se reparte por igual
  entre f_bola = 3N−3 grados de libertad internos de la bola y
  f_sup = 3N_c grados de libertad "activos" de la superficie, la fracción
  de la bola es f_bola/(f_bola+f_sup). Eso es todo: conteo de grados de
  libertad, puro espíritu de equipartición.
- **Por qué 3N−3 y no 3N−6:** un cúmulo de N átomos tiene 3N grados de
  libertad; se restan 3 de traslación del CM (esos son α_cm, no son
  "internos"). Las 3 rotaciones SÍ se cuentan como internas porque un
  impacto no perfectamente centrado excita rotación. Para N grande la
  diferencia entre 3N−3 y 3N−6 es irrelevante; para N=2 o 7 no lo es, y
  por eso esos tamaños son solo controles cualitativos.
- **N_c = número efectivo de átomos de la superficie que participan.**
  No es un parámetro microscópico medible directamente; se CALIBRA con
  los datos. El resultado interesante: si usas la fracción interna TOTAL
  (α_bola = α_TB + α_ΔVB), la calibración da **N_c ≈ N**. Lectura física:
  cada átomo de la bola "recluta" aproximadamente un átomo del piso.
  Plausible porque la bola de radio ~N^(1/3) penetra una profundidad
  comparable a su radio: volumen activado ~N, no solo el área de
  contacto ~N^(2/3).
- **La predicción estrella:** con N_c = N, la fórmula da
  α_bola → (3N)/(3N+3N) = **1/2 cuando N → ∞**. Una nanopelota grande,
  fría y sin adhesión absorbe como energía interna la MITAD de su energía
  de impacto. Límite limpio, memorizable, falsable.
- **La salvedad honesta (dilo tú antes de que la encuentre Gonzalo):**
  los datos de N=79 y 87 dan α_bola = 0.543, LEVEMENTE POR ENCIMA de la
  cota 1/2. Estrictamente eso implica N_c < N (invertir la fórmula da
  N_c/N = 1.09, 0.83, 0.83 para N=55, 79, 87). Origen probable: (i) la
  energía rotacional está incluida en α_bola pero la rotación no es un
  oscilador (no tiene PE asociada, "cuenta distinto" en equipartición), o
  (ii) α_surf está algo subestimado. Distinguirlas requiere medir el
  momento angular post-rebote (compute angmom en LAMMPS) — está declarado
  como trabajo pendiente.

### A7. h* = T_B/(Mg): qué es y qué NO es

- Idea: si al primer rebote α_cm ≈ α_TB (mira la tabla: 0.304 vs 0.319 en
  la fila h₀ᵉᶠ=25.9), entonces la altura del primer rebote
  h₁ = α_cm·h₀ᵉᶠ ≈ α_TB·h₀ᵉᶠ = T_B/(Mg). Es una COINCIDENCIA NUMÉRICA
  del primer rebote elevada a escala útil: conociendo el calentamiento
  interno predices la altura. Funciona: h₁/h* = 0.96.
- Lo que NO es: un **atractor**. "Atractor" = estado al que la dinámica
  converge y en el que se queda. Uno podría conjeturar que las alturas de
  rebote convergen a h*(T_B) porque T_B "recarga" el rebote. La
  simulación larga (5000τ, 23 rebotes) muestra lo contrario: la razón
  h_n/h*_n cae de 0.96 a ~0.18. La altura decae geométricamente
  (cada rebote pierde un factor e²) mientras que T_B, una vez depositada,
  casi no cambia (baja solo ×0.7 en 23 rebotes). Son dos relojes
  distintos: h_n decae rápido, T_B lento. Al final la bola queda
  ATRAPADA vibrando a menos de 1σ del piso — el análogo exacto del
  "prrrrt" final de una pelota real.
- El "colapso" en t≈713τ (la bola se cuela ENTRE los átomos del piso) es
  un **artefacto**: con interacción puramente repulsiva y sin atracción
  lateral, los átomos de la bola terminan filtrándose por los huecos de
  la red. Lo decimos explícitamente en el paper. Saber distinguir física
  de artefacto es parte del mensaje pedagógico.

### A8. Adhesión (c_ps = 1)

- Restauramos el término atractivo bola–piso del LJ. Resultado: en los
  24 casos (4 tamaños × 6 alturas) la bola queda pegada al primer
  contacto. No es sorpresa: la energía de adhesión medida
  (E_adh = 52–209 ε según N) es mucho mayor que la energía CM disponible
  tras el rebote (α_cm·E₀ ~ 5–50 ε).
- E_adh ∝ N^0.73 ≈ N^(2/3) = escalado del ÁREA de contacto (los átomos
  del interior no tocan el piso). Compara con la masa ∝ N. De ahí la
  predicción contraintuitiva: el umbral de caída para lograr rebote
  escala h₀* ∝ N^0.73/N = N^(−0.27) — **los cúmulos grandes escapan más
  fácil** (más energía por unidad de adhesión).
- El valor absoluto del umbral (~10²–10³ σ) es SOLO orden de magnitud:
  a esas velocidades la bola se funde al impactar (lo vio el propio
  trabajo 2012 vía g(r)), el área de contacto crecería con v, y α_cm
  está tomado del caso sin adhesión. Preséntalo siempre con esa triple
  salvedad.
- Conexión con literatura: "soft landing" de cúmulos sobre superficies
  (deposición de nanopartículas sin destruirlas) — Awasthi et al., PRB.

### A9. Los detalles técnicos que dan credibilidad

1. **Cúmulos de N exacto**: en 2012 la bola se creaba recortando una
   esfera de una red → N variable (79–84) y bola no relajada. Ahora:
   capas FCC completas (13, 55, 79, 87...) minimizadas en vacío → N
   exacto, bola inicialmente fría. Sin esto, la medición de α_TB estaría
   contaminada con energía interna inicial espuria.
2. **Potencial truncado Y DESPLAZADO (shift)**: el LJ se corta en
   r_c=2.5σ. Si solo lo cortas, U(r_c)=−0.0163ε ≠ 0, y cada par de átomos
   que cruza r_c hace saltar la energía potencial REGISTRADA en ±0.0163ε.
   Con miles de cruces por impacto, el balance de energía fallaba en
   2–6% de E₀ — lo descubrimos porque el cierre Σα=1 no cuadraba al medir
   los 4 canales de forma independiente. La solución estándar: desplazar
   el potencial (U−U(r_c)), que NO cambia las fuerzas (la derivada es la
   misma) y por tanto NO cambia ninguna trayectoria — solo arregla la
   contabilidad. Con shift: conservación a |ΔE| < 0.03ε ≈ 4×10⁻⁴E₀ por
   impacto y cierre < 10⁻⁴. Detalle importante: los resultados
   anteriores eran válidos porque α_ΔVB se calculaba por residuo.
3. **Digitalización de la Fig. 7 de 2012**: extrajimos las curvas del PDF
   del capítulo por color de píxel, calibrando con los 6 valores
   conocidos del rebote 0 (error rms de calibración: 0.08σ). Hallazgo:
   en 2012 h₀ es la ALTURA ABSOLUTA del CM (las curvas convergen al
   nivel de reposo ≈10.3σ), mientras que nuestro h₀ se mide desde el tope
   de la caja. Con la MISMA etiqueta h₀, las caídas efectivas difieren en
   ~8σ. Comparando a caída efectiva igual: acuerdo excelente a baja
   energía (e₁: 0.90 de ellos vs 0.92 nuestro) y 2012 queda ~0.06–0.09
   abajo a alta energía (explicable: bola sin minimizar + N variable).
4. **Reproducibilidad**: cada figura del paper la genera un script del
   repo desde los datos crudos. Cualquier número que Gonzalo cuestione,
   se puede regenerar delante de él.

---

## PARTE B — Guion lámina por lámina (~15 min total)

### Lámina 1 — Portada (15 s)
"Te traigo el estado del trabajo del rebote: qué reprodujimos del
capítulo 2012, qué encontramos nuevo, y una propuesta concreta de
publicación que quiero cerrar contigo hoy."

### Lámina 2 — De dónde venimos (1 min)
- Señala el bloque de la derecha: "la pregunta que ordena todo el paper
  es esta: e < 1 te dice CUÁNTA energía no vuelve; nosotros respondemos
  DÓNDE queda."
- Recorre los 5 bullets de la izquierda rápido; no te detengas en
  detalles técnicos (vuelven en la lámina 10).
- Frase puente: "el resultado central cabe en una imagen — es la
  siguiente."

### Lámina 3 — La imagen OVITO (2 min; deja que la mire)
- Deja 10 segundos de silencio; la figura se explica sola.
- Recorrido: "(a) cae — toda la energía es mecánica del CM. (b)
  compresión máxima — mira la barra: la energía está como deformación,
  54% potencial intra-bola, y el canal gris U_c es el contacto, que solo
  existe durante la colisión. (c) al separarse ya está caliente — mira
  el moteado de colores en la bola. (d) apogeo: 31% cinética interna,
  31% potencial interna, 6% superficie, 32% altura recuperada."
- Remate: "o sea: la bola se lleva ADENTRO el 62% de la energía; el piso
  se queda con el 6%. La intuición de 'se disipa en el suelo' falla."
- Detalle si pregunta por los colores: bola coloreada por KE interna
  (sin el CM), superficie por KE con escala 30 veces más sensible para
  ver los fonones.

### Lámina 4 — Los cuatro canales (2 min)
- "Esto es lo mismo cuantificado para seis alturas de caída."
- Señala la columna α_TB: "aquí está el plateau: sobre v≈0.9, la
  fracción cinética interna se clava en 0.21–0.24. La fracción deja de
  depender de cuán fuerte golpeas — señal de que el contacto alcanza a
  termalizar" (si pide detalle → tienes A5).
- Señala α_ΔVB: "y este es el canal que nadie había separado: la
  expansión anarmónica — dilatación térmica instantánea de la bola. Es
  tan grande como el cinético; es el virial actuando."
- Señala la fila 1: "impacto lento = casi elástico, 85% recuperado."
- Cierre con el bloque rojo: "suma exacta 1, cierre medido mejor que
  10⁻⁴ — la contabilidad es completa."

### Lámina 5 — Validación contra 2012 (2 min; lámina delicada, tono de colaboración)
- "Antes de extender, validamos contra ustedes. Como el capítulo no
  publica e₁ numérico, digitalicé la Fig. 7 por color, calibrando con
  las alturas iniciales conocidas."
- El hallazgo de convención: "encontré que en el capítulo h₀ es la altura
  absoluta del CM — las curvas parten en h₀ y descansan en ~10.3σ. La
  nuestra se mide desde el tope de la caja. Misma etiqueta, caída
  distinta; hay que comparar a caída efectiva igual."
- "Así comparado: 0.90 ustedes, 0.92 nosotros en el régimen suave —
  excelente. A alta energía ustedes quedan un poco abajo, consistente
  con la bola sin minimizar y el N variable del create_atoms."
- NO digas "error de ustedes"; di "diferencia de convención documentada".

### Lámina 6 — Atribución (2 min; LA lámina de la conversación)
- Léela casi literal: "esto lo necesito confirmar contigo antes de
  someter. El capítulo establece los canales, las alturas por rebote, el
  rol de U_c durante el contacto y la fusión. El borrador presenta como
  contribución nueva: el plateau cuantificado α_TB≈0.22, la fórmula de
  equipartición con N_c, y la escala h*."
- Pregunta directa: "¿esto coincide con tu recuerdo? ¿Algo de esto viene
  de notas tuyas de la época que debamos citar?"
- Escucha. Si algo era de él, la solución es simple y ya está prevista:
  se cita como comunicación privada o pasa a la introducción como idea
  previa. No es un conflicto, es limpieza de autoría.

### Lámina 7 — Equipartición (2.5 min; la lámina más "mecánica estadística")
- Panel (a): "cada curva es un tamaño de bola. Los grandes (55, 79, 87)
  suben y se aplanan — ese es el plateau. El de 13 (cuadrados abiertos
  azules) nunca se aplana: demasiado rígido para termalizar en 2τ de
  contacto. Se comporta como molécula, no como sólido."
- Panel (b): "el valor del plateau contra N. La fórmula es puro conteo
  de grados de libertad: 3N−3 de la bola contra 3N_c del piso."
- El límite: "calibrando con α_bola total sale N_c≈N — cada átomo de la
  bola recluta uno del piso — y eso da el límite α_bola → 1/2."
- La salvedad ANTES de que la encuentre él: "los datos de 79 y 87 dan
  0.543, un pelo sobre 1/2. Lo decimos en el paper: probablemente la
  rotación, que está en α_bola pero no es un oscilador. Medir el momento
  angular es el siguiente paso."

### Lámina 8 — h* (1.5 min)
- "Del plateau sale una escala: h* = T_B/Mg. Primer rebote: 0.96 —
  funciona. ¿Es un atractor? No: la razón cae a 0.18 porque la altura
  decae geométricamente y T_B casi no cambia. La bola termina atrapada
  a menos de 1σ tras 23 rebotes."
- Artefacto: "el desplome final a través de la red es un artefacto del
  potencial repulsivo puro y está documentado como tal."

### Lámina 9 — Adhesión (1.5 min)
- "Con la atracción encendida, pega siempre: 24 de 24."
- "E_adh escala como N^0.73 ≈ área de contacto; la masa como N. Por eso
  los grandes escapan más fácil: umbral ∝ N^(−0.27)."
- "El valor absoluto del umbral es solo orden de magnitud — a esas
  velocidades la bola se funde, como mostraron ustedes con g(r)."

### Lámina 10 — Rigor y límites (1.5 min)
- Columna izquierda, el shift en una frase: "descubrimos que sin
  desplazar el potencial truncado, la contabilidad energética falla en
  2–6% por impacto — puro artefacto de bookkeeping, las trayectorias no
  cambian. Con shift, cierre exacto."
- Columna derecha, sé el primero en decirlo: "T=0, una orientación, una
  trayectoria por punto: aún no hay barras de error estadísticas. Eso
  define la extensión natural del trabajo, y está declarado en el paper."

### Lámina 11 — Plan editorial (1.5 min)
- "Según lo que conversaron, apunté a la RBEF, sección Artigos Gerais:
  física computacional con punto de vista didáctico. El manuscrito ya
  cumple el protocolo: bilingüe, referencias en orden, formato de la
  revista, y agregué una sección de actividades para el aula."
- "Además existe la versión de investigación en doble columna (español e
  inglés) si después queremos un segundo paper con la estadística
  completa."
- Flujo: tu revisión → arXiv (physics.ed-ph) → submissão → repo público.

### Lámina 12 — Decisiones (2 min; cierre operativo)
Recorre las 5 y anota sus respuestas:
1. Atribución (viene de la lámina 6).
2. Idioma del envío (español está aceptado por la RBEF).
3. Su ORCID + visto bueno a arXiv.
4. Autor correspondiente y orden de autores.
5. ¿Guardamos la versión PRE para más adelante con la estadística?

---

## PARTE C — Preguntas probables de Gonzalo (y tus respuestas)

**C1. "¿Equipartición? Esto no está en equilibrio. ¿Con qué derecho usan
el teorema?"** (la pregunta más probable de un mecánico estadístico)
> Respuesta: correcto, no hay equilibrio ni ensamble — hay UNA colisión.
> Usamos equipartición como *ansatz de reparto*: si el tiempo de contacto
> basta para redistribuir la energía entre los modos dinámicamente
> accesibles, cada uno recibe en promedio lo mismo. No lo suponemos: lo
> CONTRASTAMOS. La evidencia de que el reparto ocurre es el plateau (la
> fracción deja de depender de la energía inyectada) y de que NO ocurre
> cuando no puede ocurrir: N=13, rígido, no plateau. La fórmula es
> fenomenológica con un parámetro calibrado (N_c); su mérito es predecir
> el ESCALADO con N y el límite 1/2.

**C2. "¿Y por qué la energía cinética interna equivale a la mitad de la
interna total? ¿Virial?"**
> Sí: para osciladores, ⟨KE⟩=⟨PE⟩. Por eso α_TB ≈ α_ΔVB en el plateau
> (0.22 vs 0.21–0.33; la PE fluctúa más porque incluye la deformación
> anarmónica neta, no solo la oscilación).

**C3. "¿Cómo definen temperatura de la bola si son 79 átomos?"**
> No la definimos: reportamos T_B, la KE interna, que es una cantidad
> mecánica bien definida sin hipótesis de equilibrio. Si uno QUIERE una
> temperatura instantánea, T = 2·T_B/(3N−3)k_B, pero el paper trabaja
> siempre con fracciones de energía, nunca con temperatura.

**C4. "3N−3 y no 3N−6, ¿por qué cuentan la rotación como modo interno?"**
> Porque un impacto no perfectamente central transfiere momento angular;
> esa energía no está en el CM ni en la superficie, queda 'en la bola'.
> La contamos como interna, y decimos explícitamente que la rotación no
> es un oscilador — probablemente explica el pequeño exceso sobre 1/2.
> Para N≥13 el cambio 3N−3 → 3N−6 mueve la fórmula menos que las otras
> incertezas.

**C5. "El e₁ no es monótono con h₀ (0.55 en 25.9σ pero 0.63 en 41.8σ).
¿Eso es física o ruido?"**
> No lo sabemos aún, y lo decimos: con T=0 y una sola orientación hay
> UNA trayectoria por punto, sin estadística. Puede ser resonancia con la
> fase de la vibración interna al llegar al piso, o varianza de condición
> inicial. Es la primera prioridad del trabajo futuro (3–5 réplicas con
> orientaciones aleatorias). Nota: la serie de ustedes en 2012 también
> tiene un no-monótono en esa zona (e₁=0.43 en su tercera altura).

**C6. "¿La superficie no es demasiado chica? ¿El fonón no rebota en el
borde y vuelve?"**
> La losa es 8×8×6 celdas con las 2 capas de abajo congeladas y bordes
> periódicos laterales. En el tiempo hasta el primer apogeo (~25τ) la
> onda alcanza a recorrer la caja, así que α_surf incluye energía que
> re-circula por los bordes periódicos — posible subestimación del flujo
> saliente 'verdadero' de una superficie semi-infinita. Es una limitación
> conocida de las MD de impacto; agrandar la losa y medir α_surf(L) es
> un chequeo pendiente razonable. [Sé honesto: no lo hemos barrido.]

**C7. "¿Termostato? ¿Cómo controlan T?"**
> No hay termostato: NVE puro, T inicial = 0 exacto (velocidades cero,
> cúmulo minimizado). Toda la energía viene de la caída. Por eso la
> conservación de energía es EL control de calidad, y por eso el asunto
> del shift importaba tanto.

**C8. "¿Paso de tiempo? ¿Se conserva la energía en el impacto?"**
> dt = 5×10⁻³τ, velocity-Verlet. Validado midiendo la energía total
> (incluida la PE gravitatoria) a lo largo de un impacto completo:
> |ΔE| < 0.03ε ≈ 4×10⁻⁴ de E₀. Ustedes usaron Beeman con dt 10 veces
> menor; con Verlet y potencial desplazado no hace falta.

**C9. "¿Por qué la superficie se lleva tan poco? Uno esperaría acoplo
fuerte al continuo de fonones."**
> Dos razones: (i) impedancia — el contacto dura ~2τ y excita
> preferentemente los modos de la huella, no el continuo; (ii) el canal
> dominante durante el contacto (U_c) es potencial de contacto, y esa
> energía REGRESA mayoritariamente a la bola al separarse (es lo que se
> ve en la barra gris de la lámina 3). El análisis de 2012, hecho durante
> el contacto, veía U_c grande; nuestra partición al apogeo muestra su
> destino final. Ambas descripciones son compatibles: miden en instantes
> distintos.

**C10. "¿Y esto es publicable en RBEF? Es una revista de enseñanza."**
> Sección 'Artigos Gerais': acepta física teórica/computacional/
> experimental 'bajo un punto de vista didáctico, con temas poco
> discutidos en textos clásicos, a nivel accesible a no especialistas'
> — es exactamente este manuscrito. Sin costo de publicación, SciELO,
> acepta español. Y la versión de investigación queda disponible para
> una revista de física computacional cuando tengamos la estadística de
> orientaciones (ahí sí PRE o EPJ B).

**C11. "¿Qué pasa con N intermedios? ¿19, 43?"**
> Son las capas completas que faltan entre 13 y 55, justo donde ocurre
> la transición molécula→sólido. Están en la lista de trabajo futuro;
> los cúmulos son triviales de generar con el mismo script.

**C12. "¿El N=2 con α_surf=0.18? ¿Dónde se fue el resto?"**
> e₁=0.90 → α_cm=0.81; α_surf=0.18; α_bola≈0.003. Suma ≈ 1. El dímero
> vertical casi no acopla el golpe a su único modo de vibración (que es
> perpendicular al movimiento si cae 'acostado'): todo lo que no vuelve
> al CM se lo lleva el piso. Es el único caso donde la intuición
> macroscópica ('el piso absorbe') es correcta — porque la bola no tiene
> dónde guardar energía.

---

### Chuleta de números (para tener a mano)

| Cantidad | Valor |
|---|---|
| Plateau α_TB (N=79) | 0.21–0.24, arranca en v≈0.9 σ/τ |
| α_ΔVB | 0.21–0.33 |
| α_surf | 0.07–0.16 |
| α_bola (N=79,87, plateau) | 0.543 (cota teórica 1/2) |
| N_c calibrado | ≈N (1.09 / 0.83 / 0.83 para 55/79/87) |
| h₁/h* primer rebote | 0.96; cae a ~0.18 en n=20–23 |
| Atrapamiento | 23 rebotes, t≈660τ; colapso (artefacto) t≈713τ |
| E_adh | 52–209 ε; ∝N^0.73; umbral ∝N^(−0.27) |
| Conservación con shift | \|ΔE\|<0.03ε ≈ 4×10⁻⁴E₀; cierre <10⁻⁴ |
| Validación vs 2012 | e₁ 0.90 vs 0.92 (baja E); −0.06 a −0.09 (alta E) |
| Costo | ~1.5 min/rebote en notebook, 1 núcleo |
