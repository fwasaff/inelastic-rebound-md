# Plan de reunión con Gonzalo — presentación de resultados y cierre pre-envío

**Objetivo:** que Gonzalo evalúe el trabajo completo (física, análisis,
conclusiones), discutir los resultados punto por punto, y salir de la reunión
con las 4 decisiones que bloquean el envío a RBEF.

**Material:** `presentacion_gonzalo.pdf` (Beamer, ~15 láminas) + 2 videos
OVITO + manuscrito impreso (`borrador_rbef_v2.pdf`) + repo abierto en el
notebook por si pide ver un script.

---

## 1. Estado del proyecto (para tener claro antes de entrar)

**El manuscrito v2 está científicamente listo.** Los cinco puntos mayores del
arbitraje simulado (equipartición inconsistente, N_c con tres calibraciones,
sin barras de error, tamaño de losa, fusión) están resueltos con simulaciones
nuevas, no con redacción. Verificado contra los datos:

- `tabla_orient.txt` (210 corridas) coincide con las Tablas 1 y 2 del tex.
- `Nc_medido.npy` coincide con las columnas N_c^geom / N_c^part.
- Conservación de energía re-verificada hoy: |ΔE| = 0.023 ε = 4.3×10⁻⁴ E₀
  en 400 τ (figura nueva `figuras_reunion/fig_conservacion.pdf`).

**Lo que falta para enviar (nada es cómputo):**

1. **Bloqueante:** ORCID de Gonzalo + fijar fecha de envío → *sale de esta
   reunión*.
2. Menores (≈1 h total): frase excluyendo N=7 del ajuste (misma lógica que
   N=13) o reportar χ²/dof; quitar "surprising" del abstract inglés; enunciar
   el criterio de plateau; recortar el resumen ~30%; nota al pie sobre c_ij;
   párrafo de público objetivo (P1 del árbitro 2).

---

## 2. Agenda propuesta (60 min)

| min | Tema | Lámina(s) |
|-----|------|-----------|
| 0–5   | Contexto: qué retomamos de 2012 y qué preguntamos ahora | 2 |
| 5–15  | **Video 1** (rebote sólido) + los 4 canales + cierre verificado | 3–5 |
| 15–25 | Control de calidad: conservación, shift, termalización, validación contra su Fig. 7 | 6–8 |
| 25–40 | **El resultado central:** f_ball → 1, N_c^ef ≈ 25, las tres escalas de N_c | 9–11 |
| 40–47 | **Video 2** (fusión a v≈2) + límites: fusión, losa, N=7/13 | 12–13 |
| 47–55 | Rebotes sucesivos (h*), adhesión | 14 |
| 55–60 | Decisiones: autoría/atribución, ORCID, idioma, flujo arXiv→RBEF | 15 |

**Estrategia narrativa:** contar la historia *incluyendo el error corregido*.
El argumento viejo (α_ball → 1/2 con N_c = N) era inconsistente; el arbitraje
simulado lo detectó; la medición directa de N_c lo descartó y la reformulación
sobre la energía disipada da un resultado más limpio y más interesante
(f_ball → 1). Presentarlo así tiene dos ventajas con Gonzalo: demuestra
control metodológico, y se anticipa a la objeción que él mismo haría (es la
objeción natural de cualquier físico estadístico).

**Regla práctica:** los videos primero enganchan, pero cada afirmación
cuantitativa debe volver a una figura con barras de error. Si Gonzalo quiere
profundizar en un punto, dejar la agenda y seguirlo — las láminas de respaldo
y este documento cubren las derivadas.

---

## 3. Preguntas que Gonzalo probablemente hará — respuestas preparadas

### Sobre la comparación con su trabajo de 2012

**P: "¿Por qué sus e₁ no coinciden con mi Fig. 7?"**
R: Convención de h₀. En 2012, h₀ es la altura *absoluta* del CM
(z_rest ≈ 10.3 σ); nosotros medimos desde el tope de la losa. A caída
efectiva igual, coincidimos en el régimen cuasi-elástico (0.90 vs 0.92) y
dentro de ~0.1 a alta energía; el resto se explica por cúmulo no minimizado
(N variable 79–84) e integrador distinto (Beeman). Esto está en la Sec. 3c
del manuscrito como "lección de convenciones".

**P: "El 2012 ya mostraba que U_c domina durante el contacto. ¿Qué agregan?"**
R: El 2012 describe el balance *durante* el contacto; nosotros medimos la
partición *al primer apogeo*, cuando los canales ya se separaron. U_c es
transitorio y regresa mayoritariamente a la bola. Lo nuevo: los 4 α con
cierre verificado, el plateau α_TB ≈ 0.22, la reformulación de equipartición
sobre energía disipada, N_c medido, h*, y el caso con adhesión.

**P: "¿La fórmula de equipartición y h* vienen de mis notas?"**
R: Pregunta abierta de autoría — es una de las decisiones de la reunión
(lámina final). Si algo viene de notas suyas no publicadas, se cita como
corresponde.

### Sobre la física del resultado central

**P: "¿Cómo puede N_c^ef ser ~25 y constante si la huella crece con N?"**
R: Es la distinción clave entre las tres escalas. Cientos de átomos se mueven
*transitoriamente* durante el contacto (N_c^part ≈ 300–410, medido), la
huella geométrica es 45–90 (medido), pero la superficie es casi elástica y
*devuelve* la mayor parte: el sumidero *neto* equivale a solo ~25 osciladores
(ajustado de f_ball(N)). El viejo N_c ≈ 290 coincidía por accidente con
N_c^part: era el número correcto en el balance equivocado.

**P: "¿N_c^ef es medido o ajustado?"**
R: Ajustado — un solo parámetro global sobre f_ball(N), con N_c^geom y
N_c^part como mediciones independientes que lo acotan. Ser explícito aquí:
el arbitraje simulado ya nos marcó no sobrevender el "ya no es ajuste".

**P: "¿Por qué equipartición si esto no es un sistema en equilibrio?"**
R: Es equipartición *aproximada* en el subsistema que sí termaliza durante el
contacto (~10 τ, varias decenas de períodos vibracionales). Dos evidencias:
T_B ≈ ΔV_B (virial de un sólido casi armónico, visible en
fig_conservacion b), y el éxito del conteo de modos en f_ball(N) para
N ≥ 55. Los fallos (N=7, N=13) son justamente los casos que *no* alcanzan a
termalizar — la excepción confirma el mecanismo.

**P: "¿Y por qué N=7 sí parece seguir la fórmula y N=13 no?"** *(o al revés)*
R: N=7 tiene residuo 1.5σ (el mayor del ajuste) y N=13 está fuera de
tendencia; ambos son régimen molecular. La v3 los tratará simétricamente
(excluirlos del ajuste o reportar χ²/dof). Con N ≥ 55 el acuerdo es
cuantitativo; por debajo, cualitativo.

**P: "f_ball → 1 implica que la superficie no importa para bolas grandes.
¿No es contraintuitivo? ¿Vale en el límite macroscópico?"**
R: Vale para *este* régimen: bola fría, sin adhesión, superficie cristalina
casi elástica, un solo impacto. No extrapolamos a pelotas macroscópicas
(plasticidad, defectos, adhesión cambian el sumidero). El punto físico es que
el reservorio interno crece como 3N−3 y el neto de la superficie no crece.

### Sobre el rigor numérico

**P: "¿Se conserva la energía?"**
R: |ΔE| < 0.023 ε = 4.3×10⁻⁴ E₀ por impacto (dt = 0.005 τ), figura dedicada.
α_ΔVB se reporta por cierre; la verificación independiente con V_B medido
muestra que lo único ausente es la KE lateral del CM (mediana 5×10⁻⁴ E₀,
< 0.02 E₀ para N≥55; declarada en el manuscrito). En el caso canónico
N=79/h₀=28.15 el cierre directo es < 10⁻⁴.

**P: "¿El shift del potencial no cambia la dinámica?"**
R: No: desplaza U por una constante por par dentro de r_c; las fuerzas son
idénticas. Sin shift, cada par que cruza r_c salta ±0.016 ε y el cierre
falla 2–6% de E₀. Trayectorias iguales, contabilidad limpia.

**P: "¿La losa 8×8×6 no es demasiado chica? Los fonones se reflejan."**
R: Correcto, y está cuantificado: con losa 16×16×12, α_surf sube de 0.083 a
0.143 (~0.06 E₀ de sesgo) y f_ball baja 0.87 → 0.75. La conclusión no cambia;
el sesgo está declarado en el manuscrito junto a los datos. Punto débil
conocido: es *un* caso; si lo pide, un segundo caso (otro N u otra altura)
son ~30 min de cómputo.

**P: "¿La bola sigue sólida cuando aplican el conteo de modos?"**
R: Sí hasta v ≈ 1.7 (T_ef ≤ 0.20 vs fusión ~0.30–0.35; g(r) con picos FCC,
contraste ~11). A v ≈ 2, g(r) colapsa (contraste ~2): fundida — coincide con
la fusión que el 2012 reportó. El plateau se restringe a v ≲ 1.7. Nota: a
v=2.03, T_ef = 0.283 queda bajo el rango de fusión; el criterio decisivo es
g(r), no T_ef (aclararlo si pregunta).

**P: "¿Barras de error?"**
R: 5 orientaciones × 6 alturas × 7 tamaños = 210 corridas; todo se reporta
como media ± σ. La dispersión entre orientaciones es grande (es el mensaje:
sin promediar, las no-monotonías parecen física y son ruido).

**P: "¿Momento angular post-rebote? ¿Rotación de la bola?"**
R: No medido — declarado como trabajo futuro (`compute angmom` es una
línea). Parte de la dispersión entre orientaciones puede ser canal
rotacional; no afecta el cierre (está dentro de T_B).

**P: "¿Qué pasa con temperatura inicial finita?"**
R: Futuro. T=0 aísla el mecanismo (toda la energía viene de la caída).
Con T>0 el plateau debería sobrevivir mientras T << T_fusión.

**P: "¿Y ese colapso de la bola a través de la red en la simulación larga?"**
R: Artefacto del potencial puramente repulsivo (sin barrera lateral),
documentado como tal en el texto y usado como ejercicio de análisis crítico.

### Sobre la publicación

**P: "¿Por qué RBEF y no una revista de investigación?"**
R: El ángulo pedagógico es genuino (todo corre en un notebook, actividades
incluidas, repo público) y RBEF tiene sección para esto, sin costos, acceso
abierto. La versión de investigación (doble columna, estadística completa)
existe y puede ir después a una revista de física computacional — decisión 5.

---

## 4. Decisiones que deben salir de la reunión

1. **Atribución:** plateau α_TB, fórmula de equipartición y h* —
   ¿contribución nueva o notas previas suyas? (afecta redacción de la intro)
2. **ORCID de Gonzalo** — único bloqueante logístico del envío.
3. **Idioma** del envío (español es aceptado por RBEF).
4. **Orden de autores / autor de correspondencia** (hoy: Wasaff¹, Gutiérrez²).
5. **Flujo:** ¿arXiv (physics.ed-ph) simultáneo con la submissão? ¿Se guarda
   la versión extendida para una segunda revista?

## 5. Checklist del día anterior

- [ ] Compilar `presentacion_gonzalo.pdf` y probarla en el proyector/notebook.
- [ ] Probar los 2 videos (`docs/videos_reunion/`) con sonido/reproductor.
- [ ] Imprimir `borrador_rbef_v2.pdf` (9 págs) y esta pauta.
- [ ] `git push` — repo al día por si revisa desde su oficina.
- [ ] Tener a mano: `tabla_orient.txt`, `fig_conservacion.pdf`,
      `fig_termalizacion.pdf`, informe de arbitraje simulado (por si pregunta
      cómo se estresó el manuscrito).
