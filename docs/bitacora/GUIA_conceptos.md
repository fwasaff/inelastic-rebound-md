# Guía de conceptos: qué medimos, con qué palabras, y qué significa

> Para leer antes de la bitácora. La bitácora usa estos símbolos dando por
> supuesto que se entienden; aquí se explican de cero.
> Escrita el 2026-08-09.

---

## 1. La pregunta, en una página

La bola cae, golpea la superficie y rebota. Parte de la energía se queda en la
superficie. **¿En qué forma se queda?**

Hay dos formas extremas, y son físicamente distintas:

- **Sonido (movimiento coherente).** Los átomos se mueven *juntos*, de forma
  organizada, como una ola. Esa energía todavía puede hacer trabajo: si pusieras
  algo en el camino de la ola, lo empujaría. En termodinámica esto se llama
  **exergía** — energía con capacidad de hacer trabajo.
- **Calor (movimiento térmico).** Los átomos se agitan cada uno por su lado, sin
  organización. Esa energía ya no puede hacer trabajo neto: los empujones se
  cancelan entre sí.

La diferencia **no está en la cantidad de energía**, sino en su *organización*.
La misma cantidad de energía puede estar en una u otra forma, y el paso de la
primera a la segunda es irreversible (segunda ley).

**La hipótesis de Gonzalo:** el impacto produce un **frente coherente** (una ola
que viaja hacia afuera) y detrás de él una **estela que se termaliza** (la ola se
desorganiza y se vuelve calor). Y la pregunta cuantitativa es **cuánto tarda**
esa desorganización.

**Por qué es difícil:** ver "energía" es fácil (basta sumar ½mv²). Ver
*organización* es lo difícil, porque hay que decidir, con un número, si los
átomos se mueven juntos o cada uno por su lado. Toda la dificultad técnica de
estos dos meses está ahí.

---

## 2. Las dos "τ" que se confunden

Esto genera confusión constante y conviene fijarlo:

| Uso | Qué es | Ejemplo |
|---|---|---|
| **τ como unidad de tiempo** | La unidad natural de tiempo del sistema LJ. Para argón, **1 τ ≈ 2.15 ps**. | "el eco vuelve a 1.8 τ" = a 3.9 picosegundos |
| **τ como constante de decaimiento** | El tiempo característico que tarda algo en caer a 1/e (≈37%) de su valor. Se escribe con subíndice: **τ_th**, τ_decaimiento. | "τ_th ≈ 0.25 τ" = decae en 0.25 unidades de tiempo |

Cuando lees "τ_th ≈ 0.25 τ", el primero es *la cantidad que se mide* y el segundo
es *la unidad en que se mide*. Igual que decir "la altura es 1.8 metros".

### Tabla de conversión de los números que aparecen

| En unidades LJ | En unidades físicas | Qué es |
|---|---|---|
| 1.8 τ | 3.9 ps | vuelta del eco del fondo |
| 2.1 τ | 4.5 ps | duración del pulso de contacto (FWHM) |
| 3.2 τ | 6.9 ps | contacto bola-superficie completo |
| 13.3 τ | 28.6 ps | ventana lateral (el frente llega al borde) |
| 1.10 σ | 0.37 nm | distancia al primer vecino |
| 4.65 σ | 1.58 nm | espesor de la losa **actual** |
| 45 σ | 15.3 nm | espesor de la losa **propuesta** |
| 8.28 σ/τ | 1309 m/s | velocidad longitudinal (volumen) |
| 4.60 σ/τ | **727 m/s** | velocidad del frente **medida** |

---

## 3. "Plateau" = **meseta**. Y por qué es *el* concepto central

**Meseta** = un tramo donde la curva se aplana, donde el resultado deja de
cambiar aunque cambies el parámetro.

Parece un detalle técnico. **No lo es: es la diferencia entre un número físico y
un número inventado.**

### La idea de "perilla"

Muchas mediciones requieren elegir algo arbitrario. Por ejemplo, para preguntar
"¿los átomos se mueven juntos?" hay que elegir **en qué región** los comparas: en
cubitos de 2 σ de lado, o de 5 σ. Eso es una **perilla**: un parámetro que tú
eliges y que la naturaleza no fija.

El problema: si giras la perilla y el resultado cambia, **el resultado es tuyo,
no de la naturaleza**.

La solución: buscar un rango de la perilla donde el resultado *no cambie*. Ese
rango plano es la **meseta**. Si existe, el valor de la meseta es física. Si no
existe, no tienes un número, tienes una opinión.

### Los dos ejemplos de este proyecto

**Ejemplo bueno — α_TB tiene meseta.** Al aumentar la velocidad de impacto,
α_TB sube y luego se queda quieto en ≈0.22 para v > 0.9 σ/τ. Ese aplanamiento es
lo que permite decir "α_TB vale 0.22": el número no depende de qué velocidad
exacta elijas dentro de la meseta.

**Ejemplo malo — χ no tiene meseta.** Al cambiar el tamaño del cubito ℓ:

| ℓ (σ) | 1.55 | 2.07 | 2.58 | 3.10 | 3.87 | 4.65 | 6.20 |
|---|---|---|---|---|---|---|---|
| χ | 0.73 | 0.62 | 0.54 | 0.42 | 0.26 | 0.19 | 0.17 |

Baja sin parar, un factor 4.4. **En ningún tramo se aplana.** Por eso χ se
descartó: no es que valga 0.42 con incertidumbre; es que no vale nada, porque el
valor lo fija tu elección de ℓ.

> **La regla, y es la que hemos aplicado todo el tiempo:** si el resultado
> depende de una perilla arbitraria, no es física, es la perilla. Antes de creer
> un número, hay que girar la perilla y ver si aguanta.

---

## 4. Los símbolos, en una tabla

| Símbolo | Se lee | Qué mide | Estado |
|---|---|---|---|
| **η** | eta | Anisotropía: ¿el movimiento apunta hacia afuera o en todas direcciones? | ❌ **Refutado.** Mide polarización, no coherencia |
| **ρ_a** | rho-sub-a | Razón de la elipse: cuánto sube y baja un átomo vs cuánto va y viene | ✅ Medido: **1.71 ± 0.62** |
| **χ** | chi (ji) | Fracción de energía que es movimiento coherente | ❌ **Sin meseta**, no reportable |
| **C(d)** | ce-de-de | ¿Se mueven igual dos átomos separados una distancia d? | ✅ Candidato vivo (pesado por amplitud) |
| **ξ** | xi (csi) | Longitud de correlación: hasta qué distancia se mueven juntos | ✅ Medido: **2.12 σ** (≈2 vecinos) |
| **τ_th** | tau-th | Tiempo de termalización: cuánto tarda la ola en volverse calor | ⬜ **Es lo que buscamos** |
| **c_L, c_R** | ce-ele, ce-erre | Velocidad del sonido longitudinal / de Rayleigh | ✅ c_R medida = 4.6 σ/τ |
| **λ** | lambda | Longitud de onda: el "largo" de la ola | ✅ λ_L ≈ 17.4 σ |
| **α** | alfa | Fracción de la energía inicial que va a cada canal | ✅ Cerrado (línea anterior) |

---

## 5. Las ondas: qué esperamos ver, físicamente

Cuando golpeas un sólido se generan **varios tipos de onda a la vez**, y viajan a
velocidades distintas. Los tres que importan aquí:

### (a) Onda longitudinal — la rápida, hacia el volumen
Los átomos se comprimen y descomprimen en la dirección en que viaja la onda
(como un resorte). Es la más rápida: **c_L ≈ 1309 m/s**. Se va hacia abajo, al
interior del material.

### (b) Onda de Rayleigh — la de superficie, la que nos importa
Vive **en la superficie** y decae rápido con la profundidad (a una longitud de
onda ya casi no queda nada). Los átomos describen **elipses**: suben-bajan a la
vez que van-vienen. Es más lenta: **c_R ≈ 0.5–0.55 × c_L**.

**Y esto es exactamente lo que medimos**, por dos vías independientes:

| Firma | Teoría (Rayleigh) | Medido |
|---|---|---|
| velocidad c_R/c_L | 0.50 – 0.55 | **0.56** |
| razón de elipse ρ_a | ≈1.5 (para ν≈0.3) | **1.71 ± 0.62** |

Dos números independientes que coinciden. **Sabemos que hay una onda de Rayleigh
y sabemos a qué velocidad va.** Eso está sólido.

### (c) Calor — que *no* es una onda
El calor no viaja como frente ordenado: se difunde. Los átomos se agitan
localmente y contagian a sus vecinos, sin dirección privilegiada. Por eso su
firma es una **longitud de correlación corta** (ξ ~ 1-2 vecinos): más allá de un
par de átomos, el movimiento ya no está relacionado.

**Ahí está el criterio, y es la idea más simple de todo esto:**

> - Onda coherente → los átomos se mueven juntos a lo largo de **media longitud
>   de onda** (aquí serían ~8 σ). ξ grande.
> - Calor → los átomos se mueven juntos solo con sus **vecinos inmediatos**
>   (~1-2 σ). ξ pequeña.
>
> Medir ξ es medir directamente qué fracción es ola y qué fracción es calor.

---

## 6. Los cuatro observables que probamos (y por qué fallaron tres)

Vale la pena verlo como secuencia, porque cada fallo enseñó algo:

### 1️⃣ η — anisotropía. **Refutado el 2026-08-09**
*La idea:* si el movimiento apunta radialmente hacia afuera, es la ola; si apunta
en todas direcciones por igual, es calor.

*Por qué falla:* eso vale para una onda longitudinal en un medio infinito. Pero
en una **superficie** el modo coherente es una **elipse**, no un movimiento
radial. Una onda de Rayleigh perfectamente coherente da η ≈ −0.06, que el
criterio lee como "calor".

*La prueba:* η medido y la predicción de una onda 100% coherente con la elipse
medida correlacionan a **0.9989**. η no contiene ninguna información más allá de
la forma de la elipse.

### 2️⃣ χ — fracción coherente. **Sin meseta**
*La idea (estándar en MD de choque):* en cada cubito, separa la velocidad en
"velocidad media del cubito" (= movimiento colectivo = ola) más "lo que se sale
de esa media" (= agitación = calor).

*Por qué falla:* necesita un tamaño de cubito ℓ que sea mucho mayor que la
distancia entre átomos y mucho menor que la longitud de onda. Ese rango no
existe aquí (ver §3).

### 3️⃣ C(d) solo dirección. **Ciego a la amplitud**
*La idea:* mirar si dos átomos separados d se mueven en la misma dirección.

*Por qué falla:* al normalizar cada velocidad a longitud 1, **un átomo quieto
pesa igual que uno del frente**. Daba un decaimiento precioso de 0.50 a 0 con
τ ≈ 0.25 τ... que eran los átomos *lejanos y quietos*, débilmente alineados
porque la bola al apoyarse deforma la superficie elásticamente. Nada que ver con
termalización.

*Y es el mismo error que le achacábamos a η*, reintroducido sin darnos cuenta.

### 4️⃣ C(d) pesada por amplitud. **✅ El candidato vivo**
*La corrección:* en vez de normalizar cada velocidad, se normaliza la suma:

```
C(d) = Σ(v_i · v_j) / Σ(|v_i| · |v_j|)
```

Así los átomos energéticos pesan más (que es lo que queremos: nos interesa dónde
está la energía), pero el resultado sigue estando acotado entre −1 y +1 — porque
la desigualdad de Cauchy-Schwarz garantiza |v_i·v_j| ≤ |v_i||v_j| **término a
término**.

Es el único que ha sobrevivido a su propia objeción.

---

## 7. Qué sabemos hoy, de verdad

### ✅ Sólido
- Hay una onda de superficie tipo **Rayleigh**, a **727 m/s**, con polarización
  elíptica ρ_a ≈ 1.7. Dos firmas independientes que concuerdan con la teoría.
- El pulso de contacto dura **2.1 τ** (FWHM) y el contacto total **3.2 τ**.
- Ambos bordes de la caja (lateral libre y fondo rígido) **reflejan todo**
  (R ≈ 1). No absorben nada.
- El eco del fondo vuelve a **1.8 τ**, medido directamente.

### ❌ No sabemos (y es la pregunta)
- **Si el frente es coherente.** No está medido de forma válida.
- **Cuánto tarda en termalizarse** (τ_th). No está medido.
- Por tanto: **no, todavía no hemos visto las ondas que predijo Gonzalo** en el
  sentido de haberlas *clasificado*. Vemos **una** onda; no sabemos aún si es la
  "coherente" de su hipótesis o ya es calor propagándose.

### 🚫 Por qué la corrida actual no puede contestarlo

Tres tiempos medidos que chocan:

```
forzamiento (la bola en contacto):  ████████████████ 3.2 τ
ventana limpia (eco del fondo):     █████████ 1.8 τ
```

**La ventana limpia es más corta que el forzamiento.** Estás intentando medir
cuánto tarda algo en relajarse mientras todavía lo estás empujando. Y en esos
1.8 τ el frente aún no se ha separado de la bola: "frente" y "estela" son la
misma región.

**Y el culpable es la profundidad, no el ancho:**

| | ventana limpia |
|---|---|
| lateral (borde a 61 σ) | **13.3 τ** |
| fondo (losa de 4.65 σ) | **1.8 τ** |

El fondo es **7.4× más restrictivo**. La superficie se agrandó a lo ancho, que
era la dirección que ya no molestaba.

---

## 8. El plan, en orden

| # | Qué | Por qué | Coste |
|---|---|---|---|
| **1** | Losa profunda con **esponja de fondo** | Es lo único que abre la ventana. Sin esto no hay medición posible | implementar + validar |
| **2** | Validar que la esponja absorbe (**R < 1%**) | Sin medirlo, "la esponja funciona" es fe. Ya nos costó una corrida entera | ~5 min |
| **3** | Una corrida de producción | Ver si en ventana limpia el frente separa | ~20 min |
| **4** | Medir **ξ(r,t)** con C(d) pesada por amplitud | El observable que sobrevivió | ~10 min |
| **5** | Repetir a **3 velocidades** | τ_th(v) es el discriminador final entre H1 y H0 | ~1 h |

### Qué decide cada hipótesis en el paso 4

- **Si ξ ≈ 2 σ en el frente** (igual que ahora, en ventana limpia): no hay frente
  coherente. La hipótesis de Gonzalo cae.
- **Si ξ crece a ~8 σ en el frente y decae a ~2 σ en la estela**: eso *es* la
  hipótesis de Gonzalo, y el ritmo de esa caída es **τ_th**.
- **Y el paso 5 es el que manda:** si τ_th se acorta al aumentar la velocidad, es
  termalización de verdad (anarmónica). Si no cambia con la velocidad, es solo
  dispersión elástica lineal disfrazada.

---

## 9. Anexo: ¿2D o 3D? — la pregunta está mal planteada, y conviene desarmarla

Hay que separar **visualizar** de **medir**. No es lo mismo.

### Para medir: (r, z) es lo correcto, y no pierde nada

El impacto es un punto sobre un plano. Eso da **simetría cilíndrica**: lo que
pasa a 10 σ hacia el norte es lo mismo que a 10 σ hacia el este. Entonces las
coordenadas naturales no son (x, y, z) sino:

- **r** = distancia horizontal al impacto
- **z** = profundidad
- (el ángulo azimutal φ es dirección de simetría)

Promediar sobre φ **no pierde información: gana estadística**, porque estás
promediando muchas copias del mismo fenómeno. Es lo contrario de perder una
dimensión.

Y para tu problema concreto —el eco del fondo— **necesitas la vista (z, t)**,
que es justo la que muestra el pulso bajar, rebotar y volver. En un render 3D eso
no se ve: ocurre *dentro* del material.

### La excepción, y es real

La simetría cilíndrica es un **supuesto no verificado**. La superficie es un
cristal FCC(001), que tiene simetría de 4 pliegues, no infinita: la onda puede ir
más rápido en la dirección [100] que en la [110]. Eso **sí** requiere el mapa
completo (x, y) y **no se ha medido nunca en este proyecto**.

Es barato comprobarlo: partir el mapa en sectores azimutales y ver si la
velocidad del frente cambia. Si no cambia, promediar sobre φ queda justificado
por dato y no por fe.

### Para visualizar: 3D sí, pero para comunicar

Un render 3D (OVITO) es **excelente para mostrarle a Gonzalo** qué está pasando,
y para tu propia intuición. Es malo para medir: no puedes leer un número de él.

### Resumen

| Quiero... | Vista |
|---|---|
| ver el eco del fondo | **(z, t)** ← el problema actual |
| ver el frente separarse de la estela | **(r, t)** |
| ver la onda de Rayleigh confinada arriba | **(r, z)** ← la clave de la corrida nueva |
| ver si el cristal es anisótropo | **(x, y)** ← nunca medido |
| enseñárselo a alguien | **render 3D** |

Tu intuición de "en 2D pierdo la dimensión extra" es correcta **para la vista
(x,y) de arriba**, que efectivamente aplasta la profundidad — y la profundidad es
justo donde está el problema. Pero (r, z) no aplasta la profundidad: la pone en
el eje vertical. Ésa es la vista que hay que mirar.
