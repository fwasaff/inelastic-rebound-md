# Pregunta activa

> Este archivo tiene UNA pregunta. Si hay dos, una de las dos es un desvío.
> Se actualiza cuando la pregunta se responde o cuando Gonzalo cambia el rumbo.

**Última actualización:** 2026-08-03

---

## La pregunta

**¿La energía que entra a la superficie es sonido coherente o calor, y dónde
está la frontera entre ambos?**

Mandato de Gonzalo tras la reunión de julio 2026: *estudiar la superficie, la
onda de sonido y el calor*. Su hipótesis explícita: **frente coherente, estela
térmica**.

Criterio de decisión (Sec. 6 de `docs/manuscrito/onda_vs_calor_fonones.tex`):
la anisotropía de velocidad

    η(r,t) = (<v_r²> − <v_perp²>) / (<v_r²> + <v_perp²>)

η=1 movimiento puramente radial (exergía, sonido) · η=0 isotropía (calor).

## Estado

**Medido, sin escribir.** `data/outputs/fig_anisotropia_eta.png` (2026-07-21) y
la versión con esponja (2026-07-22) dan:

- η ≈ 0.6–0.7 en el frente (0.75–1.0 × r_frente)
- η ≈ −0.2 en la estela (r < 0.30 × r_frente)

La hipótesis de Gonzalo se sostiene en primera lectura. **Pero
`onda_vs_calor_fonones.tex` §Aplicación todavía dice "no medido todavía"** — el
documento quedó atrás de la medición.

Contexto ya establecido en el mismo documento:
- r_frente ∝ t^p con p = 1.00–1.06 (balístico, umbral)
- r_rms ∝ t^p con p ≈ 0.70–0.75 (ponderado por energía)
- Kn = ℓ_fonón/L_limpio ~ 1–5 (Kaburaki 1998): frontera balística
- Tres regímenes propuestos: balístico → sonido difuso (Weaver 1984) → calor

## Las dos objeciones abiertas sobre el resultado

1. **η ≈ −0.2 no es η ≈ 0.** Negativo significa v_perp > v_r: dominio
   transversal, no isotropía térmica. Un campo realmente termalizado da η = 0.
   ¿Es modos de Rayleigh/transversales atrapados en la superficie, o es un
   artefacto de proyectar el radial 3D sobre una capa cuasi-2D?
2. **El frente colapsa a t ≈ 14τ** en el panel derecho. Hay que verificar que
   la ventana limpia declarada (~7.5τ) es la que se usó, y que el colapso es el
   eco del borde y no física.

## Lo siguiente

Responder (1) y (2). Si sobreviven, escribir §Aplicación con el resultado real.
Ninguna de las dos requiere una simulación nueva.

## Lo que NO es la pregunta activa

- El manuscrito RBEF (`borrador_rbef_v2.tex`): cerrado salvo la lista de la
  ronda 3 de arbitraje. Editarlo no es investigar.
- Canales de energía, N_c, f_ball, equipartición, sweep en N: terminado.
- Material de reunión, videos, comparaciones con 2012: producto de
  comunicación, no resultado.
