# Simulación de arbitraje — segunda ronda (versión revisada)

**Manuscrito:** "¿A dónde se va la energía de una pelota que rebota?"
(`borrador_rbef_v2.tex`) + carta de respuesta (`carta_respuesta_referi.md`)
**Ronda anterior:** major revision (`informe_arbitraje_simulado.md`, M1–M5 + m1–m8 + P1–P7)

Metodología: releí `borrador_rbef_v2.tex` línea por línea contra la carta de
respuesta y contra los cinco puntos mayores de la ronda 1, y recalculé a mano
el ajuste de la ecuación de equipartición con los números de la Tabla 2 para
comprobar que el cierre numérico es real y no solo narrativo.

---

## 1. Veredicto por punto mayor (ronda 1 → ronda 2)

| Punto | Estado real en v2 | Comentario |
|-------|--------------------|------------|
| M1 (equipartición inconsistente) | **Resuelto** | Reformulado sobre $1-\alpha_{\rm cm}$; ec. (7)–(8) ya no predice $\alpha_{\rm surf}\approx1/2$ porque el CM queda fuera del conteo. Argumento internamente consistente. |
| M2 (tres calibraciones de $N_c$) | **Resuelto con matiz** | Ver §2 — $N_c^{\rm geom}$ y $N_c^{\rm part}$ sí son medidas directas sin parámetros libres; $N_c^{\rm ef}$ sigue siendo un parámetro *ajustado* (se obtiene invirtiendo la ec. 8), no una tercera medición independiente. Mejora real (una sola constante global en vez de tres calibraciones inconsistentes) pero la redacción sobrevende un poco el "ya no es ajuste". |
| M3 (sin barras de error) | **Resuelto** | 210 corridas, 5 orientaciones, todas las tablas con media±σ. Bien ejecutado. |
| M4 (tamaño de losa) | **Resuelto, parcial** | Un caso ($N=79$, $v\approx1.43$) muestra $\alpha_{\rm surf}$ subestimado $\sim0.06$. Correcto reportarlo como sesgo sistemático declarado, pero es un solo punto — un árbitro exigente pedirá un segundo caso (otro $N$ u otra altura) para confirmar que $\sim0.06$ no es específico de ese punto. |
| M5 (fusión de la bola) | **Resuelto** | $T_{\rm ef}$ + $g(r)$, plateau restringido a $v\lesssim1.7$. Buen uso de dos diagnósticos. Nota menor: $T_{\rm ef}=0.283$ en $v=2.03$ queda *justo debajo* del rango de fusión citado ($0.30$–$0.35$); el criterio decisivo es el colapso de $g(r)$, no $T_{\rm ef}$ — vale aclarar esto explícitamente para que no parezca contradictorio. |

**Los cinco puntos mayores de la ronda 1 están genuinamente atendidos con
simulaciones nuevas**, no solo con redacción. Verifiqué que los scripts
(`medir_Nc.py`, `run_orient_sweep.sh`, `run_slab_test.sh`, `medir_fusion.py`,
`fig_Nc_escalas.py`) y las figuras que citan existen en el repositorio.

---

## 2. Un hallazgo propio: el ajuste de $N_c^{\rm ef}=25$ tiene un residuo grande en $N=7$

Repetí a mano la predicción de la ec. (8), $f_{\rm ball}=(3N-3)/(3N-3+75)$,
contra la Tabla 2:

| $N$ | $f_{\rm ball}$ predicho | $f_{\rm ball}$ medido | residuo (en $\sigma$) |
|-----|--------------------------|------------------------|------------------------|
| 7   | 0.19 | $0.55\pm0.23$ | **1.5σ** |
| 13  | 0.32 | $0.38\pm0.13$ | 0.4σ |
| 55  | 0.68 | $0.49\pm0.18$ | 1.1σ |
| 79  | 0.76 | $0.70\pm0.14$ | 0.4σ |
| 87  | 0.78 | $0.77\pm0.09$ | 0.1σ |

El ajuste es razonable en conjunto, pero $N=7$ se desvía 1.5σ y $N=55$ 1.1σ.
El texto ya reconoce que $N=13$ no alcanza el plateau por ser demasiado
rígido/molecular; **el mismo argumento se aplica con más razón a $N=7$**, que
es más pequeño, pero el manuscrito no lo señala — solo excluye del ajuste a
$N\le2$. Esto es lo que un Árbitro 1 competente notaría en segunda ronda: o
bien se excluye $N=7$ del ajuste con la misma justificación que $N=13$
(régimen molecular, no plateau limpio), o se reporta el $\chi^2$/dof del
ajuste y se admite que el acuerdo es cualitativo, no cuantitativo punto a
punto. Es un arreglo de una frase, no de cómputo.

---

## 3. Pendientes que el propio equipo ya identificó pero **no ejecutó**

La carta de respuesta es honesta y los lista como "ediciones de texto
pendientes... no requieren cómputo". Los revisé contra `borrador_rbef_v2.tex`
y en efecto siguen sin hacerse:

**Bloqueante para el envío (no es cuestión de gusto del árbitro, es logístico):**
- **m8** — Rutas de figura siguen siendo `../data/outputs/...` (líneas 306,
  450, 499, 534, 567): el `.tex` no compilará tal cual en el sistema de la
  revista. Falta ORCID de Gutiérrez (solo aparece el de Wasaff). `\today`
  sigue sin fijar a la fecha real de envío.
- **No hay `README`** en el repositorio (`ls README*` no encuentra nada), lo
  que hace P2 del Árbitro 2 (instrucciones de instalación) doblemente
  pendiente: ni en el texto ni en el repo.

**No bloqueante pero visible en segunda ronda:**
- **m4** — el resumen (español e inglés) sigue enumerando los cuatro valores
  numéricos de α; no se recortó el 30% sugerido.
- **m6** — el criterio de plateau ($v\lesssim1.7$) se usa mecánicamente pero
  nunca se enuncia como criterio cuantitativo en una frase.
- **m7** — la comparación digitalizada con Gutiérrez 2012 sigue solo en
  texto (Sec. 3c), no como tabla/figura suplementaria.
- **P1–P7 (Árbitro 2)** — ninguno ejecutado: no hay párrafo de público
  objetivo, no se promovió $N=2$ a "Actividad 0" analítica (sigue como
  control técnico en la Sec. 2), no hay nota al pie sobre $c_{ij}$, no hay
  citas a RBEF, no hay $t/\tau$ en el pie de la Figura 1.
- **Inconsistencia nueva, menor:** la carta afirma "eliminadas las
  repeticiones de sorprende/notable" — cierto en la versión en español, pero
  el resumen en **inglés** conserva "the result is surprising" (línea 97).
  Con m5 ya eliminado en un idioma y no en el otro, un árbitro atento a la
  simetría de los resúmenes bilingües (protocolo RBEF) lo señalará.

---

## 4. Veredicto simulado de la segunda ronda

**Árbitro 1 (física computacional):** Los puntos M1–M5 están resueltos con
simulaciones nuevas y honestas, incluyendo un caso (M4) donde reconocen que
tenían razón sobre un sesgo. Mi único comentario técnico nuevo es el residuo
grande en $N=7$ (§2 arriba) y la sobreventa leve de "$N_c$ ya no es un
ajuste" cuando $N_c^{\rm ef}$ sigue siéndolo (aunque ahora es un ajuste
*consistente*, con una sola constante en vez de tres). Con esas dos frases
corregidas, **recomiendo aceptar**.

**Árbitro 2 (enseñanza):** Los puntos pedagógicos (P1–P7) de mi informe
anterior no fueron tocados. Ninguno es grave por sí solo, pero juntos siguen
dejando la pregunta abierta de mi informe anterior: ¿para quién es este
artículo? **Recomiendo aceptación condicionada** a que se incorpore al menos
el párrafo de público objetivo (P1) y se resuelva la barrera de instalación
(P2, ni siquiera hay un README).

**Decisión editorial simulada: aceptación con revisiones menores — pero no
está lista para el envío todavía.** Las revisiones mayores de fondo ya están
hechas; lo que falta es una pasada editorial corta.

---

## 5. Checklist priorizado antes de enviar

**Bloqueante (rompe el envío o la primera impresión del editor):**
1. ~~Empaquetar las figuras junto al `.tex` y corregir las 5 rutas
   `../data/outputs/...`.~~ **Hecho (2026-07-04):** figuras copiadas a
   `docs/figuras_rbef/`, rutas corregidas, recompilado sin errores (9 págs).
2. ORCID de G. Gutiérrez y fecha fija (`\today` → fecha real de envío):
   **pendiente por decisión explícita del autor** (no se sabe si Gutiérrez
   tiene ORCID; se mantiene `\today` dinámico por ahora). Queda como el
   único punto bloqueante real que falta resolver antes del envío —
   requiere confirmar el ORCID con el coautor y fijar la fecha el día de
   la subida al sistema de la revista.
3. ~~Escribir un `README.md` mínimo~~ **Hecho (2026-07-04):** creado en la
   raíz del repositorio (instalación de LAMMPS, comando mínimo verificado
   end-to-end, estructura del repo). Resuelve P2 del Árbitro 2.

**Barato, mejora la percepción del árbitro (10–20 min cada uno):**
4. Una frase sobre $N=7$: excluirlo del ajuste de $N_c^{\rm ef}$ con la
   misma lógica que $N=13$, o reportar el $\chi^2$/dof del ajuste.
5. Igualar el resumen en inglés al español (quitar "surprising").
6. Enunciar el criterio de plateau en una frase (m6).
7. Recortar ~30% el resumen (m4).
8. Nota al pie sobre $c_{ij}$ en la ec. (2) (P6, dos líneas).
9. Un párrafo de público objetivo + prerrequisitos al final de la
   introducción (P1).

**Opcional / negociable:**
10. Mover la comparación con Gutiérrez 2012 a tabla suplementaria (m7).
11. Promover $N=2$ a "Actividad 0" analítica (P3) — valioso pero es una
    subsección nueva, no una edición de una frase; se puede diferir a una
    ronda de "minor revision" real si el tiempo aprieta.
12. Citar 1–2 artículos de RBEF sobre MD en enseñanza (P7).

Con 1–3 resuelto el manuscrito ya es técnicamente enviable; con 4–9 también
resuelto, llega a la revista sin flancos abiertos evidentes.
