# Simulación de arbitraje — Revista Brasileira de Ensino de Física

**Manuscrito:** "¿A dónde se va la energía de una pelota que rebota? Partición
energética del impacto de nanocúmulos de Lennard-Jones mediante dinámica
molecular" (borrador_rbef.tex, versión de julio 2026)

**Sección:** Artigos Gerais
**Autores:** F. Wasaff, G. Gutiérrez

Este documento simula el resultado realista de una primera ronda de revisión:
carta de decisión del editor, dos informes de árbitro (el formato típico de
RBEF: un árbitro con perfil de física computacional y otro con perfil de
enseñanza), y al final una estrategia de respuesta priorizada.

---

## 1. Carta de decisión del editor

> Estimado Dr. Wasaff:
>
> Su manuscrito ha sido evaluado por dos árbitros, cuyos informes se adjuntan.
> Ambos coinciden en que el trabajo es interesante, está bien escrito y encaja
> en el perfil de la sección Artigos Gerais. Sin embargo, el Árbitro 1 plantea
> objeciones sustantivas sobre la consistencia interna del argumento de
> equipartición y sobre la ausencia de barras de error, y el Árbitro 2 solicita
> ajustes de enfoque pedagógico.
>
> **Decisión: publicable después de revisiones mayores (major revision).**
>
> Le solicitamos una versión revisada acompañada de una carta de respuesta
> punto por punto en un plazo de 60 días. La versión revisada será enviada
> nuevamente al Árbitro 1.

Esta es la decisión más probable en la vida real: el trabajo no sería
rechazado (tiene un resultado claro, código abierto y valor pedagógico
evidente), pero tampoco pasaría en primera ronda, porque el argumento central
de equipartición tiene grietas que un árbitro competente va a detectar.

---

## 2. Informe del Árbitro 1 (perfil: física computacional / materia condensada)

### Evaluación general

El manuscrito retoma el experimento computacional de Gutiérrez et al. (2012) y
lo extiende con una contabilidad energética completa al primer apogeo y un
barrido en tamaño N. La pregunta ("¿a dónde se va la energía?") está bien
planteada, la ejecución técnica es cuidadosa (el análisis del potencial
desplazado en la Sección 3b es excelente y raramente se discute), y la
disponibilidad de código y datos es ejemplar. Recomiendo publicación **solo
después de resolver los siguientes puntos mayores**, en particular M1–M3, que
afectan la conclusión principal del artículo.

### Puntos mayores

**M1. El argumento de equipartición (ec. 6) es internamente inconsistente tal
como está escrito.** La ecuación (6) afirma que α_bola = (3N−3)/(3N−3+3N_c),
es decir, que **toda** la energía E₀ se reparte por equipartición entre los
modos internos de la bola y los de la huella. Pero los propios datos de los
autores muestran que el 30–40% de E₀ permanece en el canal del centro de masa
(Tabla 1), que no aparece en el conteo de modos. Peor aún: si la ecuación (6)
es válida para α_bola, por simetría del argumento debería predecir
α_surf = 3N_c/(3N−3+3N_c) ≈ 0.5 para N_c = N — pero el valor medido es
α_surf ≈ 0.10 (N=79), cinco veces menor. No es admisible aplicar el conteo de
modos a un canal e ignorar su predicción para el otro. Los autores deben:
(i) reformular el argumento sobre la energía **disipada** (1−α_cm) en lugar de
E₀, verificando si α_bola/(α_bola+α_surf) sigue la razón de modos — con sus
datos de N=79 esa razón es 0.84, no 0.50, así que probablemente la historia
correcta es otra —; o (ii) presentar la ecuación (6) explícitamente como
heurística de escalado y no como consecuencia de equipartición; o (iii)
explicar por qué el canal del CM queda excluido del conteo (¿retroceso
coherente?). En su forma actual, el acuerdo numérico de la Tabla 2 parece una
coincidencia de calibración.

**M2. N_c aparece con tres valores distintos en el mismo artículo.** El pie de
la Figura 2(b) usa N_c = 290 y un ajuste N_c = 99·N^0.20; el texto y la
Tabla 2 usan N_c = N. Un parámetro que se recalibra tres veces no es un
parámetro predictivo, es un ajuste. La solución limpia es **medir N_c
directamente en la simulación**: contar los átomos de superficie que reciben
energía (o que están dentro del alcance del potencial de algún átomo de la
bola) durante el contacto, y reportar N_c(N) medido junto a los valores
ajustados. Esto es trivial con los dumps que los autores ya generan y
convertiría la sección más débil del artículo en la más convincente. Mientras
tanto, el "límite notable" α_bola → 1/2 (ec. 7 del texto) es consecuencia
directa de la suposición N_c = N, no de los datos, y debe presentarse como tal.

**M3. No hay barras de error: cada punto es una sola trayectoria a T=0 con
una sola orientación del cúmulo.** La propia Tabla 1 delata el problema:
e₁ no es monótono con la altura (0.552 en h=25.9σ, luego 0.633 en h=41.8σ), y
α_TB fluctúa entre 0.20 y 0.32 dentro del supuesto "plateau". Con corridas de
1.5 minutos (según los propios autores), no hay excusa para no promediar sobre
≥5 orientaciones iniciales del cúmulo (o perturbaciones aleatorias pequeñas de
las posiciones) por cada punto. Todas las fracciones α y el valor del plateau
deben reportarse como media ± desviación estándar. Sin esto, no es posible
juzgar si la diferencia entre α_bola = 0.543 y la cota 1/2 (discutida en la
Sección 5) es significativa o ruido.

**M4. Posible contaminación de α_surf por efectos de tamaño finito de la
losa.** La losa es de 8×8×6 celdas con las dos capas inferiores fijas y
condiciones periódicas laterales. Los fonones generados en el impacto se
reflejan en el fondo rígido y se reintroducen por los bordes periódicos; si el
tiempo de tránsito es comparable al tiempo de contacto (~2τ), parte de la
energía "de la superficie" puede regresar a la bola antes de la separación,
sesgando α_surf hacia abajo y α_bola hacia arriba — justamente en la dirección
de la anomalía α_bola > 1/2. Los autores deben repetir al menos un caso
representativo (N=79, h₀ intermedia) con una losa del doble de tamaño lateral
y de profundidad (p. ej. 16×16×12) y reportar la variación de α_surf y
α_bola. La caída abrupta de α_surf de 0.254 (N=55) a 0.100 (N=79) en la
Tabla 2 — un factor 2.5 entre tamaños vecinos, sin discusión — refuerza la
sospecha de que este canal no está bien controlado.

**M5. ¿En qué estado queda la bola después del impacto?** Gutiérrez et al.
(2012) reportan que a caídas grandes la bola se funde. El plateau de α_TB se
promedia para v_imp > 1.1 σ/τ e incluye impactos hasta v ≈ 2 σ/τ. Si en parte
de ese rango la bola queda líquida, el conteo de modos "3N−3 de un sólido"
pierde sentido y el plateau mezclaría dos regímenes distintos. Los autores
deben caracterizar el estado post-impacto (función g(r) o criterio de
Lindemann — noto que el repositorio ya incluye un script `analizar_gr.py`) y
delimitar el rango de validez del análisis.

### Puntos menores

**m1.** Tabla 2: para N=1 faltan e₁ y α_surf sin explicación. Para un átomo
sobre potencial puramente repulsivo, e₁ y α_surf están perfectamente
definidos y son el caso de control más limpio (α_bola = 0 exacto). Complétese.

**m2.** La afirmación de la Actividad 5 ("la partición a igual v_imp no
cambia con g") se presenta como hecho pero no se muestra dato alguno. O se
verifica con una corrida, o se reformula como conjetura a explorar.

**m3.** Sección 4, "cuatro canales cuya suma es exactamente E₀": la suma
cierra por conservación de la energía; el mérito está en que cada canal se
mide de forma independiente. Recomiendo decirlo así para no sobrevender.

**m4.** El resumen es largo para el estándar de la revista; puede recortarse
un 30% sin pérdida (por ejemplo, los valores numéricos de cada α no necesitan
estar todos en el resumen).

**m5.** Lenguaje promocional: "El resultado sorprende", "límite notable",
"predicción falsable" (dos veces). Un uso es aceptable; sistemático, cansa.

**m6.** Figura 2(b): el símbolo abierto de N=13 se explica en el pie pero el
criterio de "no alcanza el plateau" no se define cuantitativamente en el
texto. Defínase (p. ej., pendiente de α_TB(v) en el último tramo).

**m7.** La comparación con la ref. [5] vía digitalización de figuras
(Sección 3c) es valiosa; sugiero mostrarla como figura o tabla suplementaria
en lugar de solo citar dos números en el texto.

**m8.** Falta el ORCID del segundo autor; la fecha usa `\today` (fíjese la
fecha de submissão); las rutas de las figuras (`../data/outputs/...`) delatan
que el .tex no compilará en el sistema de la revista — empaquete las figuras
junto al fuente.

---

## 3. Informe del Árbitro 2 (perfil: enseñanza de la física)

### Evaluación general

Artículo bienvenido: toma un problema universalmente conocido (la pelota que
rebota), lo convierte en un laboratorio computacional genuino y entrega todo
el material para reproducirlo. La sección de actividades propuestas es de lo
mejor del manuscrito. Mi recomendación es **aceptar con revisiones
moderadas**, condicionada a mejorar el andamiaje pedagógico.

### Puntos

**P1. ¿Para quién es este artículo?** El texto oscila entre un artículo de
investigación (barrido en N, exponentes de escalado, umbral de adhesión) y
uno de enseñanza (actividades, "lecciones de simulación"). Falta un párrafo
temprano que declare el público objetivo (¿física computacional de pregrado?
¿estadística? ¿trabajo de titulación?) y los prerrequisitos concretos
(mecánica newtoniana, nociones de equipartición, uso básico de terminal
Linux). Sugiero además indicar el tiempo estimado de cada actividad de la
Sección 8, como se hace en los buenos artículos de laboratorio.

**P2. La barrera de entrada real es instalar y correr LAMMPS.** Para el
estudiante típico de habla hispana/portuguesa, el obstáculo no es la física
sino la herramienta. Recomiendo un apéndice breve (o un README destacado en el
repositorio, citado explícitamente) con: instalación en Linux/Windows,
comando mínimo para correr un caso, y qué archivo mirar primero. Dos o tres
frases en el texto principal bastan si el repositorio lo cubre.

**P3. Falta un peldaño analítico entre la ecuación (1) y la simulación.**
El salto de "e < 1" a "3N−3 modos" es grande para un estudiante. Un modelo
mínimo de juguete — por ejemplo, una cadena 1D de dos o tres masas con
resortes chocando contra una pared, que puede resolverse en papel y muestra
que parte de la energía queda vibrando en la "bola" — haría transparente el
mecanismo antes de la simulación completa. Podría ser una subsección corta o
una actividad 0. (El caso N=2 de los autores es exactamente esto; ¡úsenlo
pedagógicamente en lugar de relegarlo a "control del régimen molecular"!)

**P4.** Los detalles técnicos (a)–(c) de la Sección 3 son valiosos pero
interrumpen la narrativa física. Considérese moverlos a un apéndice
"Lecciones de simulación", al que el texto principal remita. (Entiendo que
puede haber resistencia: son parte del mensaje. Alternativa: recortarlos a la
mitad.)

**P5.** La Figura 1 es excelente para el aula. Sugiero añadir en el pie la
correspondencia explícita con los tiempos de la trayectoria (t/τ de cada
panel) y, si es posible, publicar la animación como material suplementario —
para enseñanza vale oro.

**P6.** En la ecuación (2) el interruptor c_ij dentro del potencial LJ no es
estándar; un estudiante que consulte un texto verá otra forma. Una nota al
pie aclarando que c=0 equivale a quedarse con la rama repulsiva evitaría
confusión.

**P7.** Considere citar uno o dos artículos previos de la propia RBEF sobre
dinámica molecular en la enseñanza, para anclar el trabajo en la conversación
de la revista. (Es una práctica habitual y los editores lo aprecian.)

---

## 4. Estrategia de respuesta recomendada (prioridad decreciente)

Lo que un autor experimentado haría con estos informes:

### Imprescindible (bloquea la aceptación)

1. **M1 — Reformular la equipartición.** Es el punto más peligroso porque
   ataca la tesis central. La salida honesta: aplicar el conteo de modos a la
   energía disipada (1−α_cm) y comparar α_bola/(α_bola+α_surf) con
   (3N−3)/(3N−3+3N_c) usando el N_c **medido** (ver punto 2). Si la razón de
   modos no reproduce 0.84, presentar la ec. (6) como ley de escalado
   empírica motivada por equipartición, y el límite 1/2 como extrapolación
   del ajuste. Costo: reescritura de la Sección 5, cero simulaciones nuevas.

2. **M2 — Medir N_c directamente.** Un script que cuente átomos de superficie
   con KE > umbral (o dentro de r_c de la bola) durante el contacto, corrido
   sobre los dumps existentes. Elimina de raíz la crítica de "tres
   calibraciones". Costo: ~1 día de análisis, sin simulaciones nuevas.

3. **M3 — Barras de error.** Con ~1.5 min por corrida: 5 orientaciones × 6
   alturas × 7 tamaños ≈ 210 corridas ≈ 6 horas de CPU. Totalmente factible.
   Rehace todas las tablas y figuras con media ± σ. Este punto además puede
   **resolver M4/M5 de rebote**: si α_bola = 0.54 ± 0.04, la "violación" de la
   cota 1/2 deja de ser un problema.

### Importante (el árbitro lo va a verificar en la segunda ronda)

4. **M4 — Test de tamaño de losa.** Una corrida con 16×16×12 para N=79 y una
   altura intermedia (~15 min de CPU). Si α_surf cambia poco, se reporta y
   cierra el punto; si cambia mucho, es un hallazgo que mejora el artículo.

5. **M5 — Estado de la bola.** Correr `analizar_gr.py` sobre los apogeos de
   los casos de mayor velocidad y reportar una frase + figura suplementaria.

6. **P1–P3 (Árbitro 2).** Párrafo de público objetivo, apéndice/README de
   instalación, y promover N=2 a modelo de juguete pedagógico. El punto P3 es
   una oportunidad, no una carga: el caso N=2 resuelto en papel (colisión de
   un dímero armónico con una pared) daría al artículo el peldaño analítico
   que le falta y probablemente una ecuación cerrada para α_bola(N=2).

### Negociable (conceder lo barato, argumentar lo caro)

7. **m1–m8 y P4–P7:** conceder casi todo (son mejoras baratas). Puntos donde
   cabe resistirse con argumentos:
   - **P4 (mover la Sección 3 a apéndice):** defendible mantenerla en el
     cuerpo — las "lecciones de simulación" son parte del objetivo (iii)
     declarado en la introducción. Ofrecer como compromiso recortarla.
   - **m2 (Actividad 5):** más barato verificarla con dos corridas
     (g y 10g a igual v_imp) que discutirla.

### Nota sobre el tono de la carta de respuesta

Responder punto por punto, citando el número del comentario, indicando el
cambio exacto y su ubicación (página/línea de la versión revisada). Agradecer
específicamente M1 y M2: es el tipo de crítica que mejora el artículo. Nunca
responder "el árbitro no entendió" — si no entendió, la redacción era mejorable.

---

## 5. Pronóstico

Con las revisiones 1–5 bien ejecutadas, el escenario probable en segunda
ronda es **aceptación con revisiones menores**. Los riesgos residuales:

- Que el N_c medido no sea compatible con ninguna de las calibraciones
  usadas — en cuyo caso el artículo cambia de "confirmamos equipartición" a
  "la equipartición simple falla de forma instructiva". Para una revista de
  enseñanza eso sigue siendo publicable, y hasta más honesto.
- Que las barras de error revelen que el "plateau" de α_TB tiene una deriva
  con v_imp. Mitigación: reportarlo tal cual; el mensaje pedagógico (la
  energía queda mayoritariamente en la bola) no depende de que el plateau sea
  perfectamente plano.
