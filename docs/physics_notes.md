# Notas de Física – Inelastic Rebound MD

## 1. Potencial Lennard-Jones

$$V(r) = 4\varepsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right]$$

En unidades reducidas LJ: $\varepsilon = \sigma = m = 1$.

| Argón (unidades físicas)  | Valor                  |
|---------------------------|------------------------|
| $\varepsilon$             | 0.0104 eV (~120 K·kB)  |
| $\sigma$                  | 3.4 Å                  |
| $\tau = \sigma\sqrt{m/\varepsilon}$ | ~2.15 ps   |
| $T^* = k_BT/\varepsilon$  | 0.1 ≈ 12 K             |

**Cutoff:** $r_c = 2.5\,\sigma$. Error de truncamiento en energía de par:
$\Delta V = 4\varepsilon[(\sigma/r_c)^{12} - (\sigma/r_c)^6] \approx -0.0163\,\varepsilon$ por par.

---

## 2. Red FCC Argón a $T^* = 0.1$

La constante de red de equilibrio para LJ-FCC en el punto de referencia del artículo de Verlet (1967):

$$a_0^* = 1.5496 \quad [\sigma]$$

Densidad de número reducida:
$$\rho^* = \frac{N\sigma^3}{V} = \frac{4}{a_0^{*3}} \approx 1.070$$

Con $\sigma = 3.4$ Å:  $a_0 = 1.5496 \times 3.4 = 5.27$ Å (muy próximo al valor experimental de Ar sólido).

---

## 3. Coeficiente de Restitución (COR)

Para impacto normal ($\theta = 0$):

$$e = \frac{|v_{f,z}|}{|v_{0,z}|} \in [0, 1]$$

Para impacto oblicuo (ángulo $\theta$ respecto a la normal):

$$v_{0,z} = -v_0 \cos\theta, \quad v_{0,x} = v_0 \sin\theta$$

La componente tangencial se conserva aproximadamente para superficies lisas:
$$v_{f,x} \approx v_{0,x}  \quad \text{(sin fricción)}$$

---

## 4. Balance de energía en el impacto

Energía cinética incidente (normal):
$$E_{k,0} = \frac{1}{2}m v_{0,z}^2$$

Energía cinética tras rebote:
$$E_{k,f} = \frac{1}{2}m v_{f,z}^2 = e^2 E_{k,0}$$

Energía transferida a la superficie (fonones):
$$\Delta E = E_{k,0}(1 - e^2)$$

Esta energía se distribuye entre:
- Fonones acústicos de largo alcance (ondas elásticas, rápidas).
- Fonones ópticos locales (vibración del punto de impacto, lentas).

La temperatura de superficie sube:
$$\Delta T^* \approx \frac{\Delta E}{(3/2) N_{\text{calor}} k_B}$$

donde $N_{\text{calor}}$ son los átomos de la capa termostatada.

---

## 5. Ley de potencias: COR vs. velocidad

Para impactos LJ a $T^* \approx 0$, la literatura reporta:

$$e(v_0) \approx A \cdot v_0^{-\beta}$$

con exponentes $\beta \in [0.1, 0.3]$ dependiendo del rango de velocidades
(régimen plástico vs. elástico). A baja velocidad ($v_0 \lesssim 1$) el rebote
es casi elástico ($e \to 1$); a alta velocidad domina la disipación ($e \to 0$).

**Referencias:**
- Glosli & McClelland, *PRL* **70**, 1960 (1993).
- Harrison *et al.*, *J. Chem. Phys.* **97**, 6573 (1992).
- Tomsic *et al.*, *J. Chem. Phys.* **119**, 7805 (2003).

---

## 6. Criterio de contacto

El proyectil "está en contacto" con la superficie cuando la distancia al átomo
de superficie más cercano es menor que el cutoff del potencial:

$$d < r_c = 2.5\,\sigma$$

En la práctica, el script LAMMPS registra `d_contact = z_{proj} - z_{surf,top}`,
y el contacto ocurre para `d_contact < 2.5`.

---

## 7. Selección del timestep

Regla de oro para MD con LJ:
$$\Delta t \lesssim 0.005\,\tau  \quad \text{(equilibración)}$$
$$\Delta t \lesssim 0.002\,\tau  \quad \text{(impacto, gradientes altos)}$$

La frecuencia de vibración máxima del sólido LJ es:
$$\omega_{\max} \approx 10-12 \quad [\tau^{-1}]$$

Para muestrear correctamente: $\Delta t < 1/(10\,\omega_{\max}) \approx 0.008\,\tau$.
El timestep de impacto $\Delta t = 0.002\,\tau$ es conservadoramente correcto.

---

## 8. Corrección de cola en slab geometry

El comando `pair_modify tail yes` aplica la corrección de Lustig (1994) para
energía y presión de largo alcance, asumiendo función de distribución radial
$g(r) = 1$ para $r > r_c$. Esta suposición es válida para sistemas bulk.

En un slab con vacío, $g(r)$ no es 1 más allá del slab. El error en la presión
normal es:
$$|\Delta P_{zz}| \sim \frac{8\pi\rho^2\varepsilon\sigma^3}{3} \cdot \left[
\frac{2}{5}\left(\frac{\sigma}{r_c}\right)^9 - \left(\frac{\sigma}{r_c}\right)^3
\right] \cdot \frac{L_z^{\text{slab}}}{L_z^{\text{box}}}$$

Para nuestros parámetros esto es $\lesssim 1\%$ del impacto energético. **Aceptable.**

---

## 9. Termostatización selectiva

| Grupo       | Integrador | Propósito                                      |
|-------------|------------|------------------------------------------------|
| Bottom 2 capas | `setforce 0 0 0` | Pared rígida; absorbe impulso lineal    |
| Capas medias   | NVT Nosé-Hoover | Baño de calor; drena fonones         |
| Capa superior  | NVE       | Superficie libre; interacción directa con proyectil |
| Proyectil      | NVE       | Sin sesgo térmico; conserva energía total       |

Esta zonificación garantiza:
1. Que la energía que entra por el impacto pueda disiparse hacia el bulk.
2. Que el proyectil no intercambie calor artificialmente con el termostato.
3. Que la superficie libre muestre la dinámica correcta sin fricción de Langevin.
