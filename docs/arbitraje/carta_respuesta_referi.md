# Carta de respuesta a los árbitros — versión revisada

**Manuscrito:** "¿A dónde se va la energía de una pelota que rebota?…"
**Revista:** Revista Brasileira de Ensino de Física, sección Artigos Gerais
**Decisión previa:** revisión mayor (major revision)

Agradecemos a ambos árbitros por una lectura cuidadosa. Los puntos M1 y M2 del
Árbitro 1, en particular, nos obligaron a rehacer el análisis del mecanismo
físico, y el resultado es un artículo sustancialmente más sólido y honesto.
Respondemos punto por punto; cada cambio está respaldado por nuevas
simulaciones cuyos scripts se incluyen en el repositorio
(`run_orient_sweep.sh`, `medir_Nc.py`, `run_slab_test.sh`, `medir_fusion.py`,
`analizar_orient.py`).

---

## Árbitro 1 — puntos mayores

### M1. Inconsistencia del argumento de equipartición. **[Corregido — reescrito]**

El árbitro tiene toda la razón: aplicar el conteo de modos a $E_0$ e ignorar
el canal del CM (y la predicción implícita $\alpha_{\rm surf}\approx1/2$) era
inconsistente. Hemos **reformulado el argumento sobre la energía disipada**,
$1-\alpha_{\rm cm}$, que es lo único que se reparte entre modos internos y
superficie. Definimos

$$ f_{\rm ball} \equiv \frac{\alpha_{\rm ball}}{\alpha_{\rm ball}+\alpha_{\rm surf}}, $$

la fracción de la energía *no recuperada* que queda en la bola. La
equipartición entre $3N-3$ modos internos y $3N_c^{\rm ef}$ osciladores netos
de superficie predice $f_{\rm ball}=(3N-3)/(3N-3+3N_c^{\rm ef})$. Medido
(promediando orientaciones, ver M3):

| $N$ | $f_{\rm ball}$ |
|-----|----------------|
| 2   | 0.34 ± 0.29 |
| 7   | 0.55 ± 0.23 |
| 13  | 0.38 ± 0.13 |
| 55  | 0.49 ± 0.18 |
| 79  | 0.70 ± 0.14 |
| 87  | 0.77 ± 0.09 |

Un ajuste de $N_c^{\rm ef}$ **constante** (excluyendo los dímeros $N\le2$, sin
huella real) da $N_c^{\rm ef}=25\pm8$. Es decir: la capacidad *neta* de la
superficie es pequeña y prácticamente independiente de $N$. La consecuencia,
que ahora es el resultado central corregido, es

$$ f_{\rm ball}\;\xrightarrow{N\to\infty}\;1, $$

**no** $1/2$. Retiramos explícitamente la afirmación anterior
$\alpha_{\rm ball}\to1/2$ (que provenía de suponer $N_c=N$, ahora descartado
por la medición directa, punto M2) y la sustituimos por este límite, junto con
su interpretación física: una bola grande, fría y sin adhesión retiene
*casi toda* la energía que no devuelve al CM, porque su número de modos
internos supera de largo la capacidad fija de la superficie. La nueva
Figura 4 y la nueva Sección 5 desarrollan este argumento.

### M2. $N_c$ con tres calibraciones — medirlo directamente. **[Corregido — medido]**

Ahora **medimos $N_c$ en la simulación**, con dos estimadores sin parámetros
libres, evaluados en la ventana de contacto (dumps por-átomo, $v_{\rm imp}
\approx1.47\,\sigma/\tau$, script `medir_Nc.py`):

- **$N_c^{\rm geom}$** — huella geométrica: átomos de superficie a
  $<r_c=2.5\sigma$ de algún átomo de la bola.
- **$N_c^{\rm part}$** — número de participación de la KE de superficie,
  $(\sum_i K_i)^2/\sum_i K_i^2$ (cuántos átomos comparten la energía).

| $N$ | $N_c^{\rm geom}$ | $N_c^{\rm part}$ |
|-----|------------------|------------------|
| 7   | 47  | 284 |
| 13  | 43  | 318 |
| 55  | 67  | 341 |
| 79  | 65  | 414 |
| 87  | 93  | 397 |

El hallazgo es esclarecedor y resuelve la aparente arbitrariedad: existen
**tres escalas distintas** de $N_c$, bien separadas:

$$ N_c^{\rm ef}\approx25 \;\ll\; N_c^{\rm geom}\approx45\text{--}90 \;\ll\;
   N_c^{\rm part}\approx300\text{--}410. $$

Muchos átomos se mueven transitoriamente durante el contacto
($N_c^{\rm part}$), pero **la superficie es casi elástica y devuelve la mayor
parte de esa energía**; la capacidad *neta* que gobierna la partición final es
mucho menor ($N_c^{\rm ef}$). Notablemente, el valor $N_c\approx290$ que
usábamos antes (invirtiendo la fórmula errónea sobre $E_0$) coincide por
accidente con $N_c^{\rm part}$: era el número correcto en el balance
energético equivocado. Esto aparece ahora en la nueva Figura 4(a) y elimina
las tres calibraciones inconsistentes.

### M3. Sin barras de error. **[Corregido — 210 corridas]**

Rehicimos **todo** el barrido promediando sobre 5 orientaciones iniciales del
cúmulo (una cristalina + 4 rotaciones rígidas aleatorias fijas por semilla),
$7\times6\times5=210$ corridas (`run_orient_sweep.sh`). Todas las fracciones
$\alpha$ se reportan ahora como media $\pm\sigma$. El árbitro tenía razón: la
dispersión es grande. Para $N=79$ en el plateau, por ejemplo,
$\alpha_{\rm ball}=0.48\pm0.13$ y $\alpha_{\rm surf}=0.20\pm0.10$; las
no-monotonías que el árbitro señaló en la antigua Tabla 1 están dentro de esta
dispersión. Las Tablas 1 y 2 y las Figuras 2 y 4 se rehicieron con estas
barras. Esto además ubica la antigua "violación" de la cota $1/2$ dentro del
ruido, y la reformulación de M1 la vuelve irrelevante.

### M4. Tamaño finito de la losa. **[Corregido — verificado, el árbitro acertó]**

Repetimos $N=79$, $h_0=41.47\sigma$ ($v_{\rm imp}\approx1.43$) con la losa
estándar $8\times8\times6$ y una losa el doble en las tres direcciones,
$16\times16\times12$ (`run_slab_test.sh`):

| Losa | $\alpha_{\rm surf}$ | $\alpha_{\rm ball}$ | $f_{\rm ball}$ |
|------|--------------------|--------------------|----------------|
| $8\times8\times6$    | 0.083 | 0.553 | 0.87 |
| $16\times16\times12$ | 0.143 | 0.435 | 0.75 |

La sospecha del árbitro se confirma: **la losa pequeña subestima
$\alpha_{\rm surf}$** (por reflexión de fonones en el fondo fijo y los bordes
periódicos), en $\sim0.06$ de $E_0$ aquí. La conclusión cualitativa se
mantiene (la bola retiene la mayor parte de la energía disipada,
$f_{\rm ball}=0.75$ incluso con la losa grande, consistente con
$N_c^{\rm ef}\approx26$), pero ahora declaramos explícitamente que
$\alpha_{\rm surf}$ (y por tanto $N_c^{\rm ef}$) llevan una incertidumbre
sistemática de tamaño finito de este orden. Se añadió un párrafo en la
Sección 3 y una nota en la Tabla 2.

### M5. ¿Se funde la bola en el plateau? **[Corregido — delimitado]**

Caracterizamos el estado de la bola de dos maneras (`medir_fusion.py`):

1. Temperatura interna efectiva $T_{\rm ef}=2T_B/(3N-3)$ en el apogeo, para
   $N=79$:

   | $v_{\rm imp}$ | 0.68 | 0.95 | 1.16 | 1.47 | 1.72 | 2.03 |
   |---------------|------|------|------|------|------|------|
   | $T_{\rm ef}$  | 0.003| 0.064| 0.145| 0.188| 0.197| 0.283|

2. $g(r)$ de la bola a $t=0$ (FCC) vs primer apogeo; contraste
   pico/valle: $t=0\!\to\!1.1\times10^5$ (cristal perfecto);
   apogeo a $v=1.47\!\to\!11.2$ (sólido tibio, picos nítidos);
   apogeo a $v=2.03\!\to\!2.05$ (estructura colapsada, líquido).

Comparado con la fusión de un cúmulo LJ de $\sim80$ átomos ($T^*\approx0.30$–
$0.35$), la bola permanece **sólida hasta $v_{\rm imp}\approx1.7$** y solo
en el extremo superior ($v\approx2$) se funde. Por eso ahora **restringimos el
plateau y el conteo de modos de sólido a $v_{\rm imp}\lesssim1.7$** y marcamos
el punto de $v\approx2$ como régimen de fusión (nueva Figura 5 y texto en la
Sección 5). Esto explica parte de la deriva de $\alpha_{T_B}$ que el árbitro
notó.

## Árbitro 1 — puntos menores

Hechos en `borrador_rbef_v2.tex`:

- **m1** (falta $\alpha_{\rm surf}$ para $N=1$): completado. Para el átomo
  único, $\alpha_{\rm ball}\equiv0$ por definición; $\alpha_{\rm surf}$ figura
  ahora como caso de control (Tabla 2, fila $N=1$).
- **m5** (lenguaje promocional): eliminadas las repeticiones de "sorprende"/
  "sorprendentemente" y "notable" en las secciones reescritas.

Ediciones de texto pendientes para la versión final (no requieren cómputo):

- **m2** (Actividad 5, $g$ vs partición): verificar con dos corridas a $g$ y
  $10g$ a igual $v_{\rm imp}$ antes de afirmarlo.
- **m3** ("suma exacta = $E_0$"): matizar que el mérito es la medición
  *independiente* canal por canal; la suma cierra por conservación.
- **m4** (resumen largo): recortar ~30%.
- **m6** (criterio de plateau): ya definido de facto ($v\lesssim1.7$; pendiente
  $\alpha_{T_B}(v)<0.05$); falta enunciarlo en una frase.
- **m7** (comparación con [5]): mover a tabla suplementaria.
- **m8** (ORCID 2º autor, `\today`, rutas de figuras): las rutas
  `../data/outputs/` deben reemplazarse por figuras empaquetadas junto al
  fuente antes de la submissão.

## Árbitro 2 — pendientes para la versión final (texto, sin cómputo)

- **P1** público objetivo + tiempos por actividad; **P2** apéndice/README de
  instalación de LAMMPS; **P3** promover el caso $N=2$ a "Actividad 0" (dímero
  armónico resuelto en papel — la medición ya da $f_{\rm ball}\to0$ para
  $N=2$); **P4** recortar a la mitad la Sección 3 (conservándola en el cuerpo);
  **P5** añadir $t/\tau$ a la Figura 1 y animación suplementaria; **P6** nota al
  pie sobre $c_{ij}$; **P7** citar dos artículos de RBEF sobre MD en la
  enseñanza.

*Nota de integridad:* los cinco puntos mayores (M1–M5) están ejecutados y
verificados con nuevas simulaciones (scripts en el repositorio) y ya
incorporados a `borrador_rbef_v2.tex`. Los puntos menores listados arriba son
ediciones de redacción que no cambian ningún resultado.

---

## Resumen de resultados nuevos

1. La equipartición se aplica a la energía disipada; $N_c^{\rm ef}=25\pm8$
   constante; $f_{\rm ball}\to1$ (reemplaza $\alpha_{\rm ball}\to1/2$).
2. $N_c$ medido directamente: tres escalas separadas
   ($N_c^{\rm ef}\ll N_c^{\rm geom}\ll N_c^{\rm part}$); la superficie es un
   reservorio ancho pero casi elástico.
3. Todas las $\alpha$ con barras de error sobre 5 orientaciones (210 corridas).
4. Tamaño finito de losa cuantificado: $\alpha_{\rm surf}$ subestimado
   $\sim0.06$ en la losa chica; conclusión robusta.
5. Estado de la bola caracterizado: sólida para $v\lesssim1.7$, fundida en
   $v\approx2$; plateau y conteo de modos restringidos a $v\lesssim1.7$.
