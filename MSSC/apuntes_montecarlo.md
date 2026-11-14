# Manual Aplicado de Métodos Monte Carlo y Computación Científica en Julia

> **Enfoque del manual:** Guía universitaria aplicada a la Física Computacional y Estadística. Prioriza la intuición geométrica y física, la deducción rigurosa sin formalismos innecesarios de teoría de la medida, la estabilidad numérica y la implementación idiomática de alto rendimiento en **Julia**.

---

## Índice

1. [La filosofía Monte Carlo en Física](#1-la-filosofía-monte-carlo-en-física)
2. [Fundamentos de Muestreo Monte Carlo](#2-fundamentos-de-muestreo-monte-carlo)
3. [Integración Monte Carlo: *Hit-or-Miss* y Crudo](#3-integración-monte-carlo-hit-or-miss-y-crudo)
4. [Muestreo Estratificado (*Stratified Sampling*)](#4-muestreo-estratificado-stratified-sampling)
5. [Método de Aceptación-Rechazo (*Rejection Sampling*)](#5-método-de-aceptación-rechazo-rejection-sampling)
6. [Muestreo por Importancia (*Importance Sampling*)](#6-muestreo-por-importancia-importance-sampling)
7. [Cadenas de Markov Monte Carlo (MCMC y Metropolis-Hastings)](#7-cadenas-de-markov-monte-carlo-mcmc-y-metropolis-hastings)
8. [Guía Maestra de Julia para Simulación Científica](#8-guía-maestra-de-julia-para-simulación-científica)
9. [Implementaciones Completas de Referencia en Julia](#9-implementaciones-completas-de-referencia-en-julia)
10. [Tabla de Decisión y Checklist Anti-Errores](#10-tabla-de-decisión-y-checklist-anti-errores)

---

## 1. La filosofía Monte Carlo en Física

En Física nos enfrentamos continuamente a dos problemas que comparten exactamente la misma raíz matemática:

1. **Calcular integrales multidimensionales:** Por ejemplo, la energía potencial total de una distribución continua de carga, una sección eficaz de dispersión o el volumen de un cuerpo geométrico en $d$ dimensiones.
2. **Calcular valores esperados estadísticos (observables):** Por ejemplo, la magnetización o la energía media de un sistema de $10^{23}$ partículas a temperatura $T$:
   $$\langle O \rangle = \frac{1}{Z} \int O(\mathbf{x})\, e^{-\beta E(\mathbf{x})}\, d^d\mathbf{x}, \qquad \beta = \frac{1}{k_B T}$$

Cuando la dimensión $d$ es pequeña ($d = 1, 2, 3$), las cuadraturas deterministas (regla del trapecio, Simpson, cuadratura gaussiana) son extremadamente precisas. Pero cuando $d$ crece ($d = 5, 20, 1000\dots$), aparece la **maldición de la dimensionalidad**: una rejilla con solo 10 puntos por eje requiere $10^{1000}$ evaluaciones.

La idea de Monte Carlo (Ulam, Von Neumann, Metropolis) consiste en sustituir la suma exhaustiva o la rejilla determinista por un **muestreo estadístico inteligente**: en lugar de recorrer todo el espacio, lanzamos sondas aleatorias distribuidas según las leyes de la probabilidad y estimamos los promedios físicos como promedios muestrales.

### Mapa mental de las técnicas

| Método | ¿Qué problema resuelve? | ¿Genera muestras de $\rho(x)$? | ¿Rechaza muestras? |
| :--- | :--- | :--- | :--- |
| **Transformada Inversa** | Muestrear $\rho(x)$ con integral analítica invertible | Sí (exactas, i.i.d.) | No |
| **Integración MC Cruda** | Calcular $\int_\Omega f(x)\,dx$ sobre dominio acotado | No (usa uniforme) | No |
| **Muestreo Estratificado** | Reducir la varianza de una integral partiendo el dominio | No (uniforme por celdas) | No |
| **Aceptación-Rechazo** | Muestrear $\rho(x)$ difícil en baja/media dimensión | Sí (exactas, i.i.d.) | **Sí** |
| **Importance Sampling** | Calcular $\int f(x)\,dx$ o $\langle O \rangle_\rho$ reduciendo varianza | No (da **pesos** $w_i$) | No |
| **MCMC (Metropolis-Hastings)** | Muestrear $\rho(\mathbf{x})$ en **muy alta dimensión** o sin conocer $Z$ | Sí (asintóticamente, correlacionadas) | Sí (repite estado actual) |

---

## 2. Fundamentos de Muestreo Monte Carlo

### 2.1. Densidad de probabilidad (PDF) y función acumulada (CDF)

Considera una variable aleatoria continua $X$ definida sobre un dominio o **soporte** $\Omega \subseteq \mathbb{R}$.

* **Función de Densidad de Probabilidad (PDF), $\rho(x)$:** Tiene el mismo significado físico que una **densidad lineal de masa** normalizada a masa total $1$ (o que $|\psi(x)|^2$ en Mecánica Cuántica). Sus unidades son $[x]^{-1}$, y no es una probabilidad sino una densidad: la probabilidad de encontrar el sistema en un intervalo infinitesimal $[x, x+dx]$ es $\rho(x)\,dx$. Debe cumplir:
  $$\rho(x) \ge 0 \quad \forall x \in \Omega, \qquad \int_\Omega \rho(x)\,dx = 1$$

* **Función de Distribución Acumulada (CDF), $F(x)$:** Es la probabilidad total acumulada desde el extremo inferior hasta $x$ (adimensional):
  $$F(x) = P(X \le x) = \int_{-\infty}^x \rho(t)\,dt \in [0, 1]$$
  Es una función monótona creciente que satisface $F'(x) = \rho(x)$.

* **Valor esperado de un observable físico $O(x)$:**
  $$\langle O \rangle_\rho = \mathbb{E}_\rho[O(X)] = \int_\Omega O(x)\,\rho(x)\,dx$$

> **Regla fundamental de consistencia:** Toda integral que contenga un factor $\rho(x)$ normalizado a $1$ es un valor esperado respecto a $\rho$. Recíprocamente, si un código genera muestras $x_i$ distribuidas según $\rho(x)$, el promedio aritmético $\frac{1}{N}\sum_{i=1}^N O(x_i)$ converge a $\int_\Omega O(x)\,\rho(x)\,dx$, **no** a $\int_\Omega O(x)\,dx$.

---

### 2.2. Ley de los Grandes Números, Teorema Central del Límite y el escalado $N^{-1/2}$

Supongamos que extraemos $N$ muestras **independientes e idénticamente distribuidas (i.i.d.)** $x_1, x_2, \dots, x_N \sim \rho(x)$ y definimos el estimador muestral:

$$\hat{O}_N = \frac{1}{N}\sum_{i=1}^N O(x_i)$$

1. **Insesgadez:** La esperanza de nuestro estimador coincide exactamente con el valor real:
   $$\mathbb{E}[\hat{O}_N] = \frac{1}{N}\sum_{i=1}^N \mathbb{E}_\rho[O(x_i)] = \langle O \rangle_\rho$$

2. **Varianza del estimador (por independencia):** Sea $\sigma_O^2 = \text{Var}_\rho(O(X)) = \langle O^2 \rangle_\rho - \langle O \rangle_\rho^2$ la varianza intrínseca del observable. Como las muestras son independientes, la varianza de la suma es la suma de las varianzas:
   $$\text{Var}(\hat{O}_N) = \text{Var}\left(\frac{1}{N}\sum_{i=1}^N O(x_i)\right) = \frac{1}{N^2}\sum_{i=1}^N \text{Var}(O(x_i)) = \frac{1}{N^2}\,(N \sigma_O^2) = \frac{\sigma_O^2}{N}$$

3. **Teorema Central del Límite (TCL):** Para $N$ suficientemente grande (y siempre que $\sigma_O^2 < \infty$), el error $\hat{O}_N - \langle O \rangle_\rho$ sigue una distribución Gaussiana centrada en cero con desviación estándar:
   $$\sigma(\hat{O}_N) = \frac{\sigma_O}{\sqrt{N}} \propto N^{-1/2}$$

#### ¿Por qué Monte Carlo gana en dimensiones altas?

En un espacio de dimensión $d$ con $N$ evaluaciones totales de la función:

| Método | Error en 1D | Error en $d$ dimensiones | ¿A partir de qué $d$ gana Monte Carlo? |
| :--- | :--- | :--- | :--- |
| Regla del Trapecio | $\mathcal{O}(N^{-2})$ | $\mathcal{O}(N^{-2/d})$ | $d > 4$ |
| Regla de Simpson | $\mathcal{O}(N^{-4})$ | $\mathcal{O}(N^{-4/d})$ | $d > 8$ |
| **Monte Carlo** | $\mathcal{O}(N^{-1/2})$ | **$\mathcal{O}(N^{-1/2})$** | **Independiente de $d$** |

Por esta razón, en un problema 5D (como estimar el volumen de una hiperesfera 5D) o en física estadística con cientos de grados de libertad, el error de Monte Carlo al graficar $\log(\sigma)$ frente a $\log(N)$ presenta siempre una recta de **pendiente $-0.5$**.

#### Cálculo práctico de la barra de error en código

En una única simulación con $N$ muestras, no conocemos el valor real $\sigma_O$, pero podemos estimarlo con la varianza muestral (con corrección de Bessel $N-1$):

$$s_O^2 = \frac{1}{N-1}\sum_{i=1}^N \left(O(x_i) - \hat{O}_N\right)^2, \qquad \text{Error estándar} = \frac{s_O}{\sqrt{N}}$$

Y si repetimos el experimento completo $R$ veces independientes obteniendo estimaciones $\hat{O}^{(1)}_N, \dots, \hat{O}^{(R)}_N$, la desviación estándar empírica entre esas $R$ realizaciones (`std(buffer)`) mide directamente $\sigma(\hat{O}_N) \approx \sigma_O / \sqrt{N}$.

---

### 2.3. Método de la Transformada Inversa

Todo generador pseudoaleatorio moderno (como `Xoshiro256++` en Julia) entrega números uniformes $u \sim U(0, 1)$. ¿Cómo transformarlos en muestras de una densidad $\rho(x)$?

**Idea intuitiva:** Si $F(x)$ comprime todo el eje real en el intervalo $[0, 1]$ asignando tramos más anchos de $[0, 1]$ allí donde $\rho(x)$ es más grande, entonces recorrer el camino inverso $x = F^{-1}(u)$ proyecta puntos uniformes de $[0, 1]$ acumulándolos precisamente donde $\rho(x)$ es alta.

**Demostración en tres líneas:** Sea $U \sim U(0, 1)$ y definamos $X = F^{-1}(U)$. Como $F$ es monótona creciente:
$$P(X \le x) = P(F^{-1}(U) \le x) = P(U \le F(x)) = F(x)$$
Luego $X$ sigue exactamente la distribución deseada.

#### Ejemplo físico: Tiempos de desintegración radiactiva

Un núcleo inestable tiene una tasa de decaimiento $\lambda$ y la densidad de probabilidad del tiempo de vida es $\rho(t) = \lambda e^{-\lambda t}$ para $t \ge 0$:

1. Calculamos la acumulada:
   $$F(t) = \int_0^t \lambda e^{-\lambda \tau}\,d\tau = 1 - e^{-\lambda t}$$
2. Igualamos a $u \sim U(0, 1)$ y despejamos $t$:
   $$u = 1 - e^{-\lambda t} \implies e^{-\lambda t} = 1 - u \implies t = -\frac{1}{\lambda}\ln(1 - u)$$

> **Detalle numérico crucial en Julia:** `rand(rng)` genera valores en el intervalo semiabierto $[0, 1)$ (puede devolver `0.0`, nunca `1.0`). Si haces `log(rand(rng))` corres el riesgo de evaluar `log(0.0) == -Inf`. Al usar `log(1.0 - rand(rng))`, el argumento pertenece a $(0, 1]$ y nunca puede valer cero.

**¿Cuándo falla este método?** Cuando la integral $F(x)$ no tiene forma cerrada (ej. la Gaussiana $\int e^{-t^2/2}dt$) o no puede invertirse analíticamente sin métodos iterativos costosos (como una distribución polinómica de grado alto o la distribución de Planck). Para esos casos necesitamos Aceptación-Rechazo o MCMC.

---

## 3. Integración Monte Carlo: *Hit-or-Miss* y Crudo

Queremos calcular numéricamente una integral definida en $d$ dimensiones:
$$I = \int_\Omega f(\mathbf{x})\,d^d\mathbf{x}$$

### 3.1. Método de *Hit-or-Miss* (Dar en el blanco)

Es la analogía de lanzar dardos al azar sobre un rectángulo (o hipercubo) de volumen conocido $V_{\text{caja}}$ que contiene por completo la región $\Omega$ cuya área o volumen queremos medir:

1. Genera $N$ puntos uniformes $\mathbf{x}_i$ dentro de la caja de volumen $V_{\text{caja}}$.
2. Cuenta cuántos puntos ($N_{\text{dentro}}$) satisfacen la condición geométrica de pertenecer a $\Omega$.
3. Estima el volumen como:
   $$\hat{V}_\Omega = V_{\text{caja}} \cdot \frac{N_{\text{dentro}}}{N}$$

#### Ejemplo físico-geométrico: Estimación de $\pi$ mediante la hiperesfera 5D

Una hiperesfera de radio $L$ en $d=5$ dimensiones tiene volumen exacto:
$$V_5(L) = \frac{\pi^{5/2}}{\Gamma(5/2 + 1)}\, L^5 = \frac{\pi^{5/2}}{\Gamma(3.5)}\, L^5$$

Si la encerramos en el hipercubo $[-L, L]^5$, el volumen de la caja es $V_{\text{caja}} = (2L)^5 = 32 L^5$. Lanzando puntos $\mathbf{x} \in [-L, L]^5$, aceptamos el punto si $r^2 = \sum_{k=1}^5 x_k^2 < L^2$. Como la razón teórica es $p = V_5(L)/V_{\text{caja}}$:
$$p = \frac{\pi^{5/2}}{32\,\Gamma(3.5)} \implies \pi = \left(32 \cdot p \cdot \Gamma(3.5)\right)^{2/5}$$

#### Limitación de *Hit-or-Miss*

Cada lanzamiento es una variable de Bernoulli con probabilidad de éxito $p$: tiene varianza $\sigma_{\text{HM}}^2 = p(1 - p)$. Si $p$ es muy pequeño (en dimensiones altas la esfera ocupa una fracción minúscula del hipercubo: en $d=10$, $p \approx 0.0025$; en $d=20$, $p \approx 2.5 \times 10^{-8}$), la mayoría de los puntos fallan y el error relativo explota como $\sim 1/\sqrt{Np}$.

---

### 3.2. Integración Monte Carlo Cruda (Teorema del Valor Medio)

Si queremos integrar una función $f(\mathbf{x})$ sobre una caja/intervalo de volumen conocido $V$ (en 1D, $V = b - a$), **nunca uses *Hit-or-Miss* en una dimensión extra**: usa el valor medio directamente.

Multiplicamos y dividimos la integral por $V$:
$$I = \int_\Omega f(\mathbf{x})\,d^d\mathbf{x} = V \int_\Omega f(\mathbf{x}) \underbrace{\frac{1}{V}}_{\rho_U(\mathbf{x})}\,d^d\mathbf{x} = V \cdot \mathbb{E}_U[f(\mathbf{X})]$$

donde $\rho_U(\mathbf{x}) = 1/V$ es la **densidad uniforme sobre $\Omega$**. El estimador Monte Carlo crudo es:

$$\hat{I}_{\text{crudo}} = \frac{V}{N} \sum_{i=1}^N f(\mathbf{x}_i), \qquad \mathbf{x}_i \sim \text{Uniforme}(\Omega)$$

En 1D sobre $[a, b]$, generamos $x_i = a + (b - a) u_i$ con $u_i \sim U(0, 1)$, y calculamos:
$$\hat{I}_{\text{crudo}} = \frac{b - a}{N} \sum_{i=1}^N f(x_i)$$

#### ¿Por qué Monte Carlo Crudo es mejor que *Hit-or-Miss* para integrar $f(x)$?

Si $0 \le f(x) \le H$ en $[a, b]$, se demuestra fácilmente que la varianza de *Hit-or-Miss* es:
$$\text{Var}(\hat{I}_{\text{HM}}) = \frac{1}{N}\left[ (b-a)H \cdot I - I^2 \right]$$
mientras que la varianza de Monte Carlo crudo es:
$$\text{Var}(\hat{I}_{\text{crudo}}) = \frac{1}{N}\left[ (b-a)\int_a^b f(x)^2\,dx - I^2 \right]$$
Como $f(x) \le H$, tenemos $f(x)^2 \le H f(x)$, de modo que **siempre** $\text{Var}(\hat{I}_{\text{crudo}}) \le \text{Var}(\hat{I}_{\text{HM}})$.

---

## 4. Muestreo Estratificado (*Stratified Sampling*)

### 4.1. El problema del "agrupamiento" aleatorio

Cuando generas puntos puramente aleatorios $x_i \sim U(a, b)$, el azar no reparte los puntos de manera homogénea: siempre aparecen zonas con cúmulos (*clusters*) y zonas vacías. Si $f(x)$ varía significativamente a lo largo de $[a, b]$, que caigan demasiados o muy pocos puntos en una región introduce fluctuaciones innecesarias en $\hat{I}$.

### 4.2. Formulación matemática

Dividimos el dominio de integración $\Omega$ en $K$ subregiones (**estratos**) disjuntas que cubren todo el espacio:
$$\Omega = \bigcup_{k=1}^K \Omega_k, \qquad \Omega_i \cap \Omega_j = \emptyset \; (i \ne j)$$

Sea $V_k = \text{Vol}(\Omega_k)$ el volumen de cada estrato, $V = \sum_{k=1}^K V_k$ el volumen total y $p_k = V_k / V$ la fracción de volumen que ocupa el estrato $k$. Por aditividad de la integral:

$$I = \sum_{k=1}^K \int_{\Omega_k} f(\mathbf{x})\,d\mathbf{x} = \sum_{k=1}^K V_k \,\mathbb{E}_{U(\Omega_k)}[f(\mathbf{X})]$$

Si dentro de cada estrato $k$ lanzamos $N_k$ muestras uniformes $\mathbf{x}_{k,1}, \dots, \mathbf{x}_{k,N_k} \sim U(\Omega_k)$ (con presupuesto total $\sum_{k=1}^K N_k = N$), el estimador estratificado es:

$$\hat{I}_{\text{est}} = \sum_{k=1}^K \frac{V_k}{N_k} \sum_{j=1}^{N_k} f(\mathbf{x}_{k,j})$$

### 4.3. Reducción de varianza: ¿por qué funciona?

Dentro de cada estrato $k$, sea $\mu_k = \mathbb{E}_{U(\Omega_k)}[f(X)]$ su media local y $\sigma_k^2 = \text{Var}_{U(\Omega_k)}(f(X))$ su varianza local. La varianza total de la función $f(X)$ sobre $\Omega$ (la que sufre Monte Carlo crudo) puede descomponerse exactamente como:

$$\sigma_{\text{total}}^2 = \underbrace{\sum_{k=1}^K p_k \sigma_k^2}_{\text{Varianza dentro de estratos}} + \underbrace{\sum_{k=1}^K p_k (\mu_k - \mu)^2}_{\text{Varianza entre medias de estratos}}$$

Si repartimos las muestras proporcionalmente al tamaño de cada estrato ($N_k = p_k N$, asignación proporcional):

$$\text{Var}(\hat{I}_{\text{est, prop}}) = \sum_{k=1}^K V_k^2 \frac{\sigma_k^2}{N_k} = \frac{V^2}{N}\sum_{k=1}^K p_k \sigma_k^2 \le \frac{V^2}{N}\sigma_{\text{total}}^2 = \text{Var}(\hat{I}_{\text{crudo}})$$

**Interpretación física:** Al fijar el número de sondas que caen en cada celda, el muestreo estratificado **anula por completo** la fluctuación debida a que caigan más puntos en la zona alta que en la zona baja de la función (elimina la varianza *entre* estratos).

> **Caso extremo en 1D (*Jittered Sampling*):** Si divides $[a, b]$ en $K = N$ subintervalos iguales y tomas 1 solo punto aleatorio dentro de cada celda, para una función diferenciable la varianza dentro de cada celda de ancho $h = (b-a)/N$ cae como $\sigma_k^2 \sim h^2 \propto N^{-2}$. Entonces la varianza total escala como $N^{-3}$ y el error escala como **$\mathcal{O}(N^{-3/2})$** en vez de $\mathcal{O}(N^{-1/2})$.

### 4.4. Asignación Óptima de Neyman

¿Podemos hacerlo aún mejor si no todos los estratos tienen la misma variabilidad? Supón que en una región $f(x)$ es casi constante ($\sigma_k \approx 0$) y en otra oscila bruscamente ($\sigma_k$ grande). No tiene sentido gastar la mitad de las muestras donde la función es plana.

Queremos minimizar la varianza del estimador estratificado:
$$\text{Var}(\hat{I}_{\text{est}}) = V^2 \sum_{k=1}^K \frac{p_k^2 \sigma_k^2}{N_k}$$
sujeto a la restricción presupuestaria $\sum_{k=1}^K N_k = N$. Aplicando multiplicadores de Lagrange:
$$\frac{\partial}{\partial N_k}\left[ V^2\sum_{j=1}^K \frac{p_j^2 \sigma_j^2}{N_j} + \lambda \left(\sum_{j=1}^K N_j - N\right) \right] = 0 \implies -\frac{V^2 p_k^2 \sigma_k^2}{N_k^2} + \lambda = 0 \implies N_k \propto p_k \sigma_k$$

Imponiendo $\sum_{k=1}^K N_k = N$, obtenemos la **asignación óptima de Neyman**:

$$N_k^* = N \cdot \frac{p_k \sigma_k}{\sum_{j=1}^K p_j \sigma_j}$$

En la práctica, se hace un pequeño sondeo piloto (por ejemplo, $5\%$ del presupuesto) para estimar $\sigma_k$ en cada estrato y luego se distribuye el $95\%$ restante según $N_k^*$.

**Limitación en alta dimensión:** Dividir un hipercubo de dimensión $d$ en $m$ trozos por eje genera $K = m^d$ estratos. Para $d = 20$ y $m = 2$ ya tienes más de un millón de estratos. En dimensiones altas se recurre al Muestreo por Hipercubo Latino (*Latin Hypercube*) o se cambia a *Importance Sampling* / MCMC.

---

## 5. Método de Aceptación-Rechazo (*Rejection Sampling*)

### 5.1. El problema y la intuición geométrica

Queremos generar muestras $x_i$ distribuidas **exactamente** según una densidad $\rho(x)$ (sin pesos, todas las muestras aceptadas valen igual), pero $F(x)$ no es invertible.

**Idea geométrica (lanzar dardos bajo la curva $\rho(x)$):**
Imagina el recinto plano limitado entre el eje horizontal y la curva $y = \rho(x)$. Si pudieras elegir puntos $(x, y)$ **uniformemente al azar sobre toda esa área**:
* La densidad de puntos cuya coordenada horizontal está en $[x, x+dx]$ es proporcional a la altura de la curva en $x$, es decir, a $\rho(x)\,dx$.
* Por tanto, descartando la coordenada $y$, la coordenada $x$ sigue exactamente $\rho(x)$.

Como no sabemos muestrear directamente debajo de la curva complicada $\rho(x)$, construimos una función techo (**envolvente**) más sencilla que sepamos muestrear y que la cubra por completo:

$$M \cdot g(x) \ge \rho(x) \quad \forall x \in \Omega$$

donde $g(x)$ es una densidad de probabilidad fácil de muestrear (por ejemplo, Uniforme en $[a, b]$, Exponencial o Gaussiana) y $M \ge 1$ es una constante positiva.

### 5.2. Algoritmo paso a paso

1. Genera una coordenada horizontal candidata $x_0 \sim g(x)$.
2. Genera una coordenada vertical aleatoria $u \sim U(0, 1)$ de forma independiente (lo que equivale a elegir una altura uniforme $y = u \cdot M g(x_0)$ entre $0$ y el techo $M g(x_0)$).
3. **Condición de aceptación:** Si el punto $(x_0, y)$ cae debajo de la curva objetivo $\rho(x_0)$, es decir:
   $$u \le \frac{\rho(x_0)}{M\,g(x_0)}$$
   entonces **aceptamos** $x_0$. En caso contrario, lo **rechazamos** y repetimos desde el paso 1.

### 5.3. Demostración de exactitud (Teorema de Bayes)

¿Por qué la distribución de los puntos aceptados es exactamente $\rho(x)$ y no $g(x)$?

Calculemos primero la probabilidad total de aceptar una propuesta en un intento cualquiera:
$$P(\text{Aceptación}) = \int_\Omega P(\text{Aceptar} \mid X = x)\,g(x)\,dx = \int_\Omega \frac{\rho(x)}{M\,g(x)}\,g(x)\,dx = \frac{1}{M}\int_\Omega \rho(x)\,dx = \frac{1}{M}$$

*(Nota física inmediata: la **eficiencia del algoritmo** es $1/M$. Si $M=4$, aceptas el $25\%$ de las propuestas; si $M=1.2$, aceptas el $83.3\%$.)*

Ahora aplicamos Bayes para hallar la densidad condicionada a que el punto haya sido aceptado:
$$\rho(x \mid \text{Aceptado}) = \frac{P(\text{Aceptar} \mid X = x)\,g(x)}{P(\text{Aceptación})} = \frac{\frac{\rho(x)}{M\,g(x)}\,g(x)}{1/M} = \rho(x)$$

La función auxiliar $g(x)$ y la constante $M$ **se cancelan exactamente**.

### 5.4. Ventaja oculta: Densidades SIN normalizar

Supón que solo conoces $\tilde{\rho}(x)$ (por ejemplo, el factor de Boltzmann $e^{-\beta E(x)}$ sin calcular la función de partición $Z$, o $\tilde{\rho}(x) = 1 - 1.2x^3$ sin calcular su integral), tal que $\rho(x) = \tilde{\rho}(x) / Z$.

Basta elegir $\tilde{M}$ tal que $\tilde{M}\,g(x) \ge \tilde{\rho}(x)$ y aceptar con probabilidad:
$$p_{\text{acc}}(x_0) = \frac{\tilde{\rho}(x_0)}{\tilde{M}\,g(x_0)}$$

Al repetir la demostración anterior, el factor $Z$ desaparece del cociente: **Aceptación-Rechazo no necesita que $\rho(x)$ esté normalizada** (la eficiencia pasa a ser $Z / \tilde{M}$).

---

### 5.5. Los 4 errores fatales en Aceptación-Rechazo (y cómo evitarlos)

#### Error 1: Truncar la propuesta con un `if` sin renormalizar $g(x)$ (El error de la Gaussiana truncada)
Supón que $\rho(x)$ solo vive en el intervalo $[a, b] = [0, 0.8]$, y decides proponer candidatos usando una normal estándar $x_0 \sim \mathcal{N}(0, 1)$, pero metes un filtro:
```julia
# ❌ INCORRECTO:
x0 = randn(rng)
if 0.0 <= x0 <= 0.8
    prob = rho(x0) / (M * g_normal(x0))
    if rand(rng) <= prob
        ...
    end
end
```
**¿Por qué está mal?**
1. La distribución de probabilidad de los $x_0$ que entran dentro del `if` ya **NO es** $g_{\text{normal}}(x)$, porque has tirado el $\sim 71\%$ de las muestras fuera de $[0, 0.8]$. Su densidad real es la **normal truncada**:
   $$g_{\text{trunc}}(x) = \frac{g_{\text{normal}}(x)}{\int_{0}^{0.8} g_{\text{normal}}(t)\,dt} = \frac{g_{\text{normal}}(x)}{\Phi(0.8) - \Phi(0)} \approx \frac{g_{\text{normal}}(x)}{0.2881}$$
2. Si aun así divides por $M \cdot g_{\text{normal}}(x_0)$ usando `if` solo como atajo, **ojo al análisis fino**:
   * Si en un intento donde $x_0 \notin [0, 0.8]$ lo cuentas como **rechazado** (equivalente a definir $\rho(x_0) = 0$ fuera de $[0, 0.8]$), entonces el algoritmo **sí** muestrea $\rho(x)$, **siempre que** $M \cdot g_{\text{normal}}(x) \ge \rho(x)$ en todo $[0, 0.8]$ y cuentes todos los intentos (incluidos los fuera de rango) en el denominador de la eficiencia ($\sim 1/M$).
   * **Sin embargo**, el fallo grave aparece en cuanto olvidas que al filtrar con un `if` estás condicionando a $[0, 0.8]$: si reutilizas esa misma lógica en *Importance Sampling* (dividiendo entre `g_normal(x0)` sin renormalizar por $0.2881$), **destruyes el resultado** multiplicándolo por $\sim 3.47$; y en Aceptación-Rechazo estás desperdiciando de entrada el $71\%$ de los números generados antes siquiera de evaluar $\rho$.
   * **La solución limpia en un intervalo finito $[a, b]$:** Usa una propuesta **Uniforme** sobre el mismo soporte $[a, b]$:
     $$g(x) = \frac{1}{b - a}, \qquad x_0 = a + (b - a)\cdot \text{rand(rng)}$$
     No desperdicias ni una sola muestra fuera del intervalo y necesitas el mínimo $M = (b - a) \max_{x \in [a,b]} \rho(x)$. Para $\rho(x) = A(1 - 1.2x^3)$ en $[0, 0.8]$ ($A = 1/0.67712 \approx 1.4768$), la envolvente uniforme óptima tiene $M = 0.8 \cdot A \approx 1.1815$, ¡con una eficiencia del **$84.6\%$** en lugar de $\sim 25\%$!

#### Error 2: Elegir $M$ demasiado pequeño ($M g(x) < \rho(x)$ en alguna zona)
Si en algún subintervalo $\rho(x) > M g(x)$, entonces $p_{\text{acc}} > 1$. Como una probabilidad no puede superar $1$, al comparar `rand() <= prob` aceptas con probabilidad $1$ (en lugar de la necesaria para compensar la falta de propuestas), recortando o "achatando" la densidad generada por debajo de la curva real.

#### Error 3: Reemplazar el generador aleatorio dentro del bucle
```julia
# ❌ CATASTRÓFICO:
for i in 1:n_samples
    rng = Xoshiro(seed) # ¡Reinicia el estado en cada iteración!
    x0 = rand(rng)      # Genera el MISMO número 1,000,000 de veces
end
```
El objeto `rng` guarda el estado interno de la secuencia pseudoaleatoria. Instanciarlo dentro del bucle destruye la aleatoriedad por completo.

#### Error 4: Confundir el cálculo de una integral usando muestras de Aceptación-Rechazo
Si ya tienes muestras aceptadas $x_1, \dots, x_{N_{\text{acc}}} \sim \rho(x)$:
* Para hallar $\langle O \rangle_\rho = \int_a^b O(x)\rho(x)\,dx$: promedias directamente $\frac{1}{N_{\text{acc}}} \sum O(x_i)$.
* Para hallar la integral pura $\int_a^b f(x)\,dx$: como $\int f(x)\,dx = \int \frac{f(x)}{\rho(x)}\rho(x)\,dx$, el estimador es $\frac{1}{N_{\text{acc}}} \sum \frac{f(x_i)}{\rho(x_i)}$ (**sin multiplicar** por $(b-a)$, porque $\rho(x)$ ya está normalizada a $1$).
* **Crítica conceptual:** Si tu objetivo final era solo calcular $\int_a^b f(x)\,dx$ y ya conoces $\rho(x)$, usar Aceptación-Rechazo para luego dividir por $\rho(x)$ es dar un rodeo ineficiente: tiras muestras a la basura y obtienes el mismo estimador que **Importance Sampling** directo usando $\rho(x)$ desde el principio, pero con menos muestras útiles.

---

## 6. Muestreo por Importancia (*Importance Sampling*)

### 6.1. Filosofía: No tires muestras, asigna pesos

En Aceptación-Rechazo descartábamos muchas muestras caras de evaluar para que las supervivientes tuvieran todas el mismo peso ($1$).
En *Importance Sampling* (IS) le damos la vuelta al problema: **aprovechamos el 100% de las muestras** generadas desde una distribución fácil $g(x)$ (llamada distribución de importancia o propuesta), y corregimos el sesgo multiplicando cada punto por un **peso estadístico**.

### 6.2. Dos formulaciones fundamentales

#### Caso A: Calcular la integral $\int_\Omega f(x)\,dx$

Sea $g(x)$ una densidad normalizada ($\int_\Omega g(x)\,dx = 1$) tal que $g(x) > 0$ en todo punto donde $f(x) \ne 0$. Multiplicamos y dividimos por $g(x)$ dentro de la integral:

$$I = \int_\Omega f(x)\,dx = \int_\Omega \frac{f(x)}{g(x)}\, g(x)\,dx = \mathbb{E}_g\left[\frac{f(X)}{g(X)}\right]$$

Generamos $N$ muestras $x_1, \dots, x_N \sim g(x)$ (¡sin rechazar ninguna!) y calculamos:

$$\hat{I}_{\text{IS}} = \frac{1}{N}\sum_{i=1}^N \frac{f(x_i)}{g(x_i)}$$

Observa qué pasa si $g(x) = \frac{1}{b-a}$ (Uniforme en $[a, b]$): la fórmula se reduce a $\frac{b-a}{N}\sum f(x_i)$, ¡recuperando Monte Carlo crudo como un caso particular!

#### Caso B: Calcular un promedio físico $\langle O \rangle_\rho$ respecto a una $\rho(x)$ difícil

Queremos calcular $\langle O \rangle_\rho = \int_\Omega O(x)\,\rho(x)\,dx$, pero no sabemos muestrear $\rho(x)$:

$$\langle O \rangle_\rho = \int_\Omega O(x)\,\frac{\rho(x)}{g(x)}\,g(x)\,dx = \mathbb{E}_g\big[O(X)\,w(X)\big], \qquad w(x) = \frac{\rho(x)}{g(x)}$$

El factor $w(x_i) = \rho(x_i)/g(x_i)$ es el **peso de importancia** de la muestra $x_i$:
* Si muestreamos en una zona donde $g(x_i) > \rho(x_i)$ (hemos sobre-muestreado esa zona), su peso $w(x_i) < 1$ compensa el exceso.
* Si muestreamos en una zona donde $g(x_i) < \rho(x_i)$ (la hemos visitado poco), su peso $w(x_i) > 1$ compensa la escasez.

#### Caso C (*Self-Normalized Importance Sampling*, SNIS): cuando tampoco conoces la constante de normalización $Z$ de $\rho(x) = \tilde{\rho}(x)/Z$

Podemos estimar $Z = \int \tilde{\rho}(x)\,dx = \int \frac{\tilde{\rho}(x)}{g(x)} g(x)\,dx \approx \frac{1}{N}\sum_{i=1}^N \tilde{w}(x_i)$ con las mismas muestras, donde $\tilde{w}(x_i) = \tilde{\rho}(x_i)/g(x_i)$. Al dividir numerador y denominador, la $N$ se cancela y obtenemos:

$$\langle O \rangle_\rho \approx \frac{\sum_{i=1}^N \tilde{w}(x_i)\,O(x_i)}{\sum_{i=1}^N \tilde{w}(x_i)}$$

---

### 6.3. Elección de $g(x)$: El secreto de la varianza mínima

La varianza de un solo término en $\hat{I}_{\text{IS}}$ es:
$$\text{Var}_g\left(\frac{f(X)}{g(X)}\right) = \int_\Omega \frac{f(x)^2}{g(x)}\,dx - I^2$$

¿Qué densidad $g^*(x)$ minimiza esta varianza? Usando la desigualdad de Cauchy-Schwarz sobre $\int |f| = \int \frac{|f|}{\sqrt{g}}\sqrt{g}$:
$$\left(\int_\Omega |f(x)|\,dx\right)^2 \le \left(\int_\Omega \frac{f(x)^2}{g(x)}\,dx\right) \underbrace{\left(\int_\Omega g(x)\,dx\right)}_{=1}$$

El límite inferior se alcanza cuando $\frac{|f(x)|}{\sqrt{g(x)}} \propto \sqrt{g(x)}$, es decir:

$$g^*(x) = \frac{|f(x)|}{\int_\Omega |f(t)|\,dt}$$

**Resultado espectacular:** Si $f(x) \ge 0$ y elegimos $g(x) \propto f(x)$, entonces el cociente $f(x)/g(x)$ es **una constante exacta igual a $I$ en todo el dominio**. La varianza es **CERO** (con 1 sola muestra ya darías el resultado exacto, si conocieras la constante de normalización!).

En la práctica real no conocemos $\int |f|$, pero la receta física es clara: **elige $g(x)$ con una forma parecida a $|f(x)|$** (por ejemplo, una recta o polinomio simple que imite la pendiente de $\cos(\frac{\pi}{2} x)$ en $[0, 0.8]$, como $1 - 1.2x^3$ o $1 - x^2$, que hace casi constante el cociente $f(x)/g(x)$ y reduce drásticamente la varianza respecto a la uniforme!).

### 6.4. ¡Cuidado con las colas ligeras! (Peligro de varianza infinita)

Existe una trampa mortal en *Importance Sampling*: elegir una $g(x)$ que decaiga a cero **más rápido** que $|f(x)|$ en los extremos (colas ligeras).

* **Ejemplo desastroso:** Intentas integrar $f(x) = \frac{1}{1+x^2}$ (colas lentas/pesadas de Lorentz) en $(-\infty, \infty)$ usando $g(x) = \mathcal{N}(0, 1) \propto e^{-x^2/2}$ (colas gaussianas ultrarrápidas).
* Casi siempre las muestras caen en $[-3, 3]$ y parece que converge, hasta que una muestra rara cae en $x = 6$. En ese punto:
  $$\frac{f(6)}{g(6)} \propto \frac{1/37}{e^{-18}} \approx 1.8 \times 10^6$$
  De repente, un único punto pega un latigazo numérico que arruina la media. Matemáticamente, $\int \frac{f(x)^2}{g(x)}\,dx = \infty$: **el estimador tiene varianza infinita y el TCL deja de cumplirse**.
* **Regla de oro:** $g(x)$ debe tener **colas más pesadas** (decaer más lento) que $f(x)$ o $\rho(x)$, o dicho al revés, el cociente $f(x)/g(x)$ debe permanecer acotado en todo $\Omega$.

### 6.5. Comparativa directa: Aceptación-Rechazo vs. Importance Sampling

| Aspecto | Aceptación-Rechazo | Importance Sampling |
| :--- | :--- | :--- |
| **Propósito principal** | Obtener **muestras no ponderadas** $x_i \sim \rho(x)$ | Estimar **valores esperados o integrales** con baja varianza |
| **¿Se tiran muestras?** | Sí (proporción $1 - 1/M$) | **No**, se aprovecha el $100\%$ |
| **Relación clave con $g(x)$** | Exige $M g(x) \ge \rho(x)$ (envolvente superior) | Exige $g(x) > 0$ donde $f(x) \ne 0$ y colas pesadas |
| **Peso de cada muestra** | $1$ (todas valen igual) | $w_i = f(x_i)/g(x_i)$ o $\rho(x_i)/g(x_i)$ |
| **Error típico de mezcla** | Calcular $\int f\,dx$ aceptando muestras de $\rho$ para luego dividir por $\rho$ | Creer que IS devuelve muestras distribuidas según $\rho$ (solo devuelve pares $(x_i, w_i)$) |

---

## 7. Cadenas de Markov Monte Carlo (MCMC y Metropolis-Hastings)

### 7.1. ¿Por qué colapsan Transformada Inversa, Aceptación-Rechazo e IS en alta dimensión?

Considera un sistema clásico de Física Estadística: un modelo de Ising con $L \times L = 10 \times 10$ espines ($d = 100$ variables binarias), o un fluido de $N = 100$ átomos en 3D ($d = 300$ coordenadas continuas). La densidad de equilibrio es Boltzmann:
$$\pi(\mathbf{x}) = \frac{1}{Z} e^{-\beta E(\mathbf{x})}$$

1. **No podemos calcular $Z$:** Implicaría integrar en $d = 300$ dimensiones.
2. **Aceptación-Rechazo es imposible:** Aunque construyas una propuesta factorizada que por cada coordenada solo necesite un factor $m_1 = 1.05$ ($95\%$ de eficiencia en 1D), en $d = 300$ dimensiones la constante global es:
   $$M = (1.05)^{300} \approx 2.27 \times 10^6 \implies \text{Eficiencia} \approx 0.00004\%$$
   ¡Prácticamente todas las configuraciones propuestas tienen energías enormes (átomos solapados) y son rechazadas!
3. **Importance Sampling también colapsa:** Una muestra uniforme en 300D tiene $e^{-\beta E} \approx 0$ el $99.99999\%$ de las veces; un solo punto afortunado acumula el $100\%$ del peso (*degeneración de pesos*).

**La solución genial de Metropolis, Rosenbluth, Rosenbluth, Teller y Teller (1953):**
En lugar de intentar generar configuraciones independientes "desde la nada", partimos de una configuración $\mathbf{x}^{(0)}$ y vamos dando **pequeños pasos** locales (moviendo una partícula o volteando un espín cada vez), diseñando las reglas de transición para que el paseo aleatorio visite cada estado $\mathbf{x}$ con frecuencia proporcional a $\pi(\mathbf{x})$.

---

### 7.2. Cadenas de Markov y Distribución Estacionaria

Una **Cadena de Markov** es una sucesión aleatoria de estados $X^{(0)}, X^{(1)}, X^{(2)}, \dots$ que satisface la **propiedad de Markov** (ausencia de memoria):
$$P\left(X^{(t+1)} = y \;\middle|\; X^{(t)} = x,\, X^{(t-1)} = x_{t-1}, \dots\right) = P(x \to y)$$

El núcleo o probabilidad de transición $P(x \to y)$ cumple $\int P(x \to y)\,dy = 1$ (desde $x$, el sistema debe terminar en algún estado $y$, que puede ser el mismo $x$).

Queremos diseñar $P(x \to y)$ de modo que tras muchos pasos, la densidad de probabilidad de encontrar el sistema en un estado $x$ converja a nuestra densidad objetivo $\pi(x)$, llamada **distribución estacionaria** (o de equilibrio):

$$\int \pi(x)\, P(x \to y)\,dx = \pi(y)$$

Para que la cadena converja de manera única a $\pi(x)$ independientemente del estado inicial $X^{(0)}$, se necesitan dos propiedades físicas básicas (**ergodicidad**):
1. **Irreducibilidad:** Desde cualquier estado posible ($\pi(x)>0$) se puede alcanzar cualquier otro estado en un número finito de pasos (el sistema no queda encerrado para siempre en un sector aislado del espacio de fases).
2. **Aperiodicidad:** La cadena no queda atrapada en ciclos rígidos deterministas (ej. $A \to B \to A \to B$).

---

### 7.3. Balance Detallado (Reversibilidad Microscópica)

Resolver directamente la ecuación integral de estacionariedad $\int \pi(x) P(x \to y)\,dx = \pi(y)$ es complicado. En su lugar imponemos una condición más fuerte y muy natural en física estadística: el **Balance Detallado**:

$$\pi(x)\,P(x \to y) = \pi(y)\,P(y \to x) \qquad \forall x, y$$

**Significado físico:** En equilibrio térmico, el flujo de probabilidad que va del microestado $x$ al microestado $y$ es exactamente igual al flujo inverso que vuelve de $y$ a $x$ (no hay corrientes de probabilidad circulares netas).

**Prueba de que el Balance Detallado implica estacionariedad:** Integra ambos lados respecto a $x$:
$$\int \pi(x)\,P(x \to y)\,dx = \int \pi(y)\,P(y \to x)\,dx = \pi(y) \underbrace{\int P(y \to x)\,dx}_{=1} = \pi(y)$$

---

### 7.4. El Algoritmo de Metropolis-Hastings paso a paso

Dividimos cada transición $x \to y$ ($x \ne y$) en dos etapas:
1. **Propuesta:** Desde el estado actual $x$, proponemos un candidato $y$ usando una distribución de propuesta sencilla $q(x \to y)$ (ej. un pequeño desplazamiento gaussiano $y = x + \Delta x$, con $\Delta x \sim \mathcal{N}(0, \delta^2)$).
2. **Aceptación:** Aceptamos el salto $x \to y$ con probabilidad $A(x \to y)$; si se rechaza, **el sistema permanece en el estado actual** $X^{(t+1)} = x$. *(¡Ojo! A diferencia de Aceptación-Rechazo, en MCMC el estado actual se vuelve a contar en la cadena cuando se rechaza un salto, ¡no se descarta la iteración!)*

La probabilidad de transición para $x \ne y$ es $P(x \to y) = q(x \to y)\,A(x \to y)$. Al sustituir en el balance detallado:

$$\pi(x)\,q(x \to y)\,A(x \to y) = \pi(y)\,q(y \to x)\,A(y \to x) \implies \frac{A(x \to y)}{A(y \to x)} = \frac{\pi(y)\,q(y \to x)}{\pi(x)\,q(x \to y)}$$

Para que la cadena explore el espacio lo más rápido posible, queremos que las probabilidades de aceptación $A$ sean tan grandes como sea posible (sin superar $1$). La elección estándar de **Metropolis-Hastings (1970)** es:

$$A(x \to y) = \min\left(1,\; \frac{\pi(y)\,q(y \to x)}{\pi(x)\,q(x \to y)}\right)$$

#### El caso estrella de la Física: Metropolis simétrico (1953)

Si la propuesta es simétrica ($q(x \to y) = q(y \to x)$, como un salto gaussiano o uniforme centrado en $x$), el factor $q$ se cancela:

$$A(x \to y) = \min\left(1,\; \frac{\pi(y)}{\pi(x)}\right)$$

Y si $\pi(x) = \frac{1}{Z} e^{-\beta E(x)}$ es la distribución canónica de Boltzmann:

$$\frac{\pi(y)}{\pi(x)} = \frac{\frac{1}{Z} e^{-\beta E(y)}}{\frac{1}{Z} e^{-\beta E(x)}} = e^{-\beta (E(y) - E(x))} = e^{-\beta \Delta E}$$

* **¡La función de partición $Z$ desaparece por completo!**
* Si el movimiento **baja la energía** ($\Delta E \le 0$), $e^{-\beta \Delta E} \ge 1 \implies$ se acepta con probabilidad $1$.
* Si el movimiento **sube la energía** ($\Delta E > 0$), se acepta con probabilidad $e^{-\beta \Delta E}$ (activación térmica que permite escapar de mínimos locales!).

---

### 7.5. Aspectos prácticos críticos en MCMC

1. **Termalización (*Burn-in*):**
   La cadena empieza en un punto arbitrario $x^{(0)}$ que probablemente está muy lejos de las zonas de equilibrio de $\pi(x)$. Los primeros $N_{\text{burn}}$ pasos sirven para que el sistema "olvide" su condición inicial y alcance el equilibrio térmico. **Esas muestras iniciales se descartan antes de medir promedios.**

2. **Autocorrelación y Tamaño Efectivo de Muestra ($N_{\text{eff}}$):**
   Las muestras sucesivas $x^{(t)}, x^{(t+1)}$ están fuertemente correlacionadas (están muy cerca en el espacio). Por eso, la fórmula ingenua del error $\sigma / \sqrt{N}$ **subestima brutalmente el error real** en MCMC!
   Si las muestras tienen un tiempo integrado de autocorrelación $\tau_{\text{int}} \ge 1/2$, el error estándar real escala como:
   $$\sigma(\hat{O}_N) \approx \sqrt{\frac{2\tau_{\text{int}}}{N}}\,\sigma_O = \frac{\sigma_O}{\sqrt{N_{\text{eff}}}}, \qquad N_{\text{eff}} = \frac{N}{2\tau_{\text{int}}}$$
   En la práctica se estima el error mediante *binning* (*blocking analysis*): agrupar los datos en bloques de tamaño creciente $B$ y observar cuándo la desviación estándar entre las medias de los bloques satura en una meseta.

3. **Calibración del tamaño de salto ($\delta$):**
   * $\delta$ **demasiado pequeño:** Casi todos los saltos se aceptan ($\sim 99\%$), pero el caminante avanza milimétricamente (movimiento difusivo lento) y $\tau_{\text{int}}$ es enorme.
   * $\delta$ **demasiado grande:** Saltas constantemente fuera de los pozos de potencial hacia zonas de energía altísima; la tasa de aceptación cae cerca de $0\%$, te quedas atascado en el mismo sitio y $\tau_{\text{int}}$ vuelve a dispararse.
   * **Regla empírica:** Ajustar $\delta$ durante la fase previa (no mientras mides!) para conseguir tasas de aceptación en torno al **$30\%$ – $50\%$** en bajas dimensiones (y $\approx 23.4\%$ en muy alta dimensión para propuestas gaussianas).

4. **Estabilidad numérica en logaritmos:**
   En física y estadística nunca calcules $\pi(y)/\pi(x)$ directamente si $\pi(x) \propto e^{-\beta E(x)}$ puede desbordar a `0.0` (*underflow*) o `Inf` (*overflow*). Trabaja siempre con el **logaritmo de la densidad** $\ell(x) = \ln \tilde{\pi}(x) = -\beta E(x)$:
   $$\ln A(x \to y) = \min\big(0,\; \ell(y) - \ell(x) + \ln q(y \to x) - \ln q(x \to y)\big)$$
   Y acepta si `log(rand(rng)) <= log_alpha` (o `rand(rng) <= exp(log_alpha)` cuando $\Delta \ell < 0$).

---

## 8. Guía Maestra de Julia para Simulación Científica

Julia está diseñado para tener la expresividad de Python/MATLAB y la velocidad de C/Fortran. Sin embargo, para alcanzar esa velocidad y evitar errores sutiles, hay que entender cómo compila Julia el código (**JIT basado en inferencia de tipos**).

### 8.1. *Scope* (Ámbito de variables) y la regla de oro: "Mete todo en funciones"

1. **Nunca hagas cálculos intensivos en el ámbito global (*global scope*):**
   Las variables declaradas sueltas en un script (fuera de toda función) son variables globales cuyo tipo podría cambiar en cualquier momento. Por eso el compilador de Julia **no puede optimizarlas** y genera código decenas de veces más lento y con miles de asignaciones de memoria.
   * **Regla práctica:** Define funciones para todo (`function simular(...) ... end`) y en el nivel global solo llama a esas funciones.

2. **Los bloques `for`, `while`, `let` crean un ámbito local:**
   Si declaras `limit_inf = 0.0` *dentro* de un bucle `for`, esa variable **deja de existir** al llegar al `end` del bucle:
   ```julia
   # ❌ ERROR (UndefVarError: `limit_inf` not defined):
   for i in 1:10
       limit_inf = 0.0
   end
   println(limit_inf)

   # ✅ CORRECTO: declárala antes del bucle o dentro de la función
   limit_inf = 0.0
   for i in 1:10
       ...
   end
   ```

3. **No olvides cerrar cada `function`, `for`, `if` con su `end`:** Julia no usa llaves `{}` ni indentación significativa (como Python); la indentación es solo estética, quien manda es `end`.

---

### 8.2. Estabilidad de Tipos (*Type Stability*)

Cuando Julia compila una función para unos argumentos concretos (por ejemplo, `Float64` e `Int64`), intenta deducir el tipo exacto de cada variable local para usar registros de CPU nativos. Si una variable cambia de tipo en mitad de la función, rompe la optimización.

#### Trampa 1: El array vacío sin tipo (`[]`)
```julia
# ❌ LENTÍSIMO: [] crea un Vector{Any}
samples = []
push!(samples, 0.42) # Cada número queda "empaquetado" (boxed) en el heap

# ✅ RÁPIDO: especifica el tipo concreto
samples = Float64[]
# Aún mejor si sabes el tamaño máximo: prealoca o reserva memoria
sizehint!(samples, n_samples)
```

#### Trampa 2: Inicializar acumuladores como enteros (`0`) en lugar de flotantes (`0.0`)
```julia
# ❌ INESTABILIDAD DE TIPO: cum empieza como Int64 y luego pasa a Float64
function mala_suma(n::Int)
    cum = 0
    for i in 1:n
        cum += cos(0.1 * i) # Aquí cambia de Int64 a Float64!
    end
    return cum
end

# ✅ ESTABLE EN TIPOS: cum es Float64 desde el principio
function buena_suma(n::Int)
    cum = 0.0
    for i in 1:n
        cum += cos(0.1 * i)
    end
    return cum
end
```

#### Cómo auditar tus funciones como un profesional: `@code_warntype`
En el REPL de Julia puedes escribir:
```julia
@code_warntype buena_suma(1000)
```
Si las variables salen marcadas en **azul/normal**, tu función tiene estabilidad de tipos perfecta (velocidad tipo C). Si ves tipos en **rojo** (`Any`, `Union{Int64, Float64}`), tienes una fuga de tipos que frena el código.

---

### 8.3. Qué son las Macros (`@...`) y cuándo usarlas

Las macros transforman tu código antes de que se compile. Sus reglas exactas son:

| Macro | Para qué sirve | Uso correcto | Error común |
| :--- | :--- | :--- | :--- |
| `@inline` | Sugiere al compilador incrustar el cuerpo de una **función** pequeña para evitar el coste de llamada | `@inline mi_func(x) = x^2` o `@inline mi_func(x)` | **NO** se pone delante de un `for` (`@inline for ...` es incorrecto) |
| `@inbounds` | Desactiva la comprobación de índices dentro de arrays | `@inbounds for i in eachindex(v) ... end` | Usarlo cuando el índice puede salirse de rango (puede corromper memoria/crashear) |
| `@simd` | Permite vectorización SIMD en bucles simples sin dependencias entre iteraciones | `@simd for i in 1:n s += a[i] end` | Usarlo en bucles con `if` complejos, `push!` o dependencias secuenciales (`rng`) |
| `@time` | Mide tiempo de ejecución y **asignaciones de memoria** | `@time simular(1_000_000)` | Medir la 1.ª ejecución (incluye tiempo de compilación JIT; ejecuta siempre 2 veces!) |
| `Threads.@threads` | Paraleliza las iteraciones de un bucle `for` entre los hilos disponibles | Ver sección 8.5 | Compartir un `rng` o hacer `push!` concurrente al mismo array |

> **Observación clave sobre `@time`:** Si un bucle Monte Carlo de $10^6$ iteraciones en Julia te muestra muchas megabytes/gigabytes en `allocations` al usar `@time`, casi siempre indica un error de diseño (inestabilidad de tipos, variable global o creación innecesaria de vectores dentro del bucle). Un buen bucle Monte Carlo con acumuladores escalares hace **`0 allocations`**.

---

### 8.4. Números aleatorios (`Random`) y precisión numérica

1. **Instanciación del generador (`Xoshiro`):**
   `Xoshiro` (Xoshiro256++) es el generador rápido y de altísima calidad por defecto en Julia moderno ($\ge 1.7$). Instáncialo **una sola vez por secuencia independiente**:
   ```julia
   using Random
   rng = Xoshiro(42) # Semilla fija => experimento 100% reproducible
   ```
2. **Muestrear en un intervalo $[a, b]$:**
   ```julia
   x = a + (b - a) * rand(rng)   # Uniforme en [a, b)
   z = randn(rng)                # Normal estándar N(0, 1)
   ```
3. **Evitar *Catastrophic Cancellation* al calcular la varianza:**
   La fórmula algebraica $\text{Var}(X) = \mathbb{E}[X^2] - (\mathbb{E}[X])^2$ parece cómoda en un solo bucle (`sum_x += x; sum_x2 += x^2`), pero es **numéricamente inestable**: si la media es $10^8$ y la desviación típica es $0.1$, restar dos números cercanos a $10^{16}$ pierde casi todos los dígitos significativos de `Float64` (incluso puede dar varianza negativa!).
   **Solución estándar:** El **algoritmo online de Welford** calcula media y varianza en una sola pasada, con `0` memoria extra y máxima precisión:
   ```julia
   mean_x = 0.0
   m2     = 0.0
   for i in 1:n
       x = ...
       delta   = x - mean_x
       mean_x += delta / i
       delta2  = x - mean_x
       m2     += delta * delta2
   end
   var_x = m2 / (n - 1) # Varianza muestral insesgada
   ```

---

### 8.5. Paralelización segura con `Threads.@threads`

Para aprovechar varios núcleos (arrancando Julia con `julia --threads auto`), recuerda dos reglas inviolables en computación paralela:

1. **Nunca compartas un objeto `rng` entre hilos:** Dos hilos modificando el estado interno del mismo `Xoshiro` a la vez producen una *race condition* (corrupción de secuencia y pérdida de reproducibilidad).
   * **Patrón seguro por realización (ideal cuando haces `reps` experimentos independientes):**
     ```julia
     buffer = zeros(Float64, reps) # Prealocación! Cada hilo escribe en su índice único `i`
     Threads.@threads for i in 1:reps
         rng_local = Xoshiro(seed + i) # Semilla única y determinista por hilo/tarea
         buffer[i] = ejecutar_experimento(N, rng_local)
     end
     ```
   * **Nota de rigor estadístico:** Al usar `Xoshiro(seed + i)`, los generadores producen flujos independientes gracias al excelente salto de estado de `Xoshiro256++`. Otra opción limpia para paralelizar dentro de un mismo experimento grande de $N$ muestras es dividir $N$ en $T =$ `Threads.nthreads()` bloques, dar a cada bloque su propio `rng` y sus acumuladores locales, y al final combinar las sumas de cada bloque.

2. **Nunca hagas `cum += ...` o `push!(v, ...)` sobre la misma variable desde varios hilos:** Cada hilo debe trabajar en sus propias variables locales o escribir en su casilla exclusiva `resultados[i]`.

---

### 8.6. Vectorización (*Broadcasting* `.`) vs. Bucles `for` en Julia

A diferencia de Python/NumPy o MATLAB, donde los bucles `for` son lentos y estás obligado a vectorizar todo en arrays enormes, en **Julia los bucles `for` bien tipados dentro de funciones son igual de rápidos que C y no consumen RAM**.

* **En Monte Carlo:** Escribir un bucle `for` con acumuladores escalares es **mejor y más rápido** que hacer `mean(f.(a .+ (b-a) .* rand(10_000_000)))`, porque la versión vectorizada reserva cientos de megabytes de RAM para los vectores temporales, mientras que el bucle `for` usa `0 bytes` de RAM y vive entero en la caché L1 de la CPU!
* Cuando sí operes sobre vectores ya existentes (por ejemplo, calcular una pendiente o graficar curvas de referencia), usa la sintaxis punto (`.`) para fusionar operaciones sin crear arrays intermedios:
  ```julia
  # Fusión automática de bucle en Julia (solo reserva el array final `y`)
  y = @. exp(ordenada + pendiente * log_N)
  ```

---

## 9. Implementaciones Completas de Referencia en Julia

Todas las funciones siguientes cumplen estabilidad de tipos, cero asignaciones en el bucle crítico, manejo correcto de `Xoshiro` y retorno estructurado (sin mezclar cálculo puro con `println` dentro del núcleo numérico).

### 9.1. Integración Monte Carlo Cruda (con error estándar vía Welford)

```julia
using Random

"""
    mc_integrate_crude(f, a::Float64, b::Float64, n_samples::Int; seed::Int=1)

Estima ∫_a^b f(x) dx mediante Monte Carlo crudo usando n_samples muestras uniformes.
Devuelve una tupla con nombre: `(integral, std_error)`.
"""
function mc_integrate_crude(f, a::Float64, b::Float64, n_samples::Int; seed::Int=1)
    @assert b > a "El límite superior debe ser mayor que el inferior"
    @assert n_samples > 1 "Se necesitan al menos 2 muestras para estimar la varianza"

    rng = Xoshiro(seed)
    ancho = b - a
    media_f = 0.0
    m2 = 0.0

    for i in 1:n_samples
        x = a + ancho * rand(rng)
        fx = Float64(f(x))

        # Algoritmo de Welford para media y varianza estables
        delta = fx - media_f
        media_f += delta / i
        m2 += delta * (fx - media_f)
    end

    integral = ancho * media_f
    var_f = m2 / (n_samples - 1)
    std_error = ancho * sqrt(var_f / n_samples)

    return (integral = integral, std_error = std_error)
end
```

---

### 9.2. Muestreo Estratificado 1D (Uniforme por estratos)

```julia
using Random

"""
    mc_integrate_stratified(f, a::Float64, b::Float64, n_strata::Int, samples_per_stratum::Int; seed::Int=1)

Divide [a, b] en `n_strata` subintervalos iguales y muestrea `samples_per_stratum` puntos
en cada uno. Reduce la varianza eliminando el agrupamiento aleatorio entre celdas.
"""
function mc_integrate_stratified(
    f, a::Float64, b::Float64, n_strata::Int, samples_per_stratum::Int; seed::Int=1
)
    @assert b > a && n_strata >= 1 && samples_per_stratum >= 1
    rng = Xoshiro(seed)

    ancho_total = b - a
    h = ancho_total / n_strata # Tamaño de cada estrato
    suma_medias_estratos = 0.0
    suma_var_estratos = 0.0

    for k in 1:n_strata
        a_k = a + (k - 1) * h
        media_k = 0.0
        m2_k = 0.0

        for j in 1:samples_per_stratum
            x = a_k + h * rand(rng)
            fx = Float64(f(x))
            delta = fx - media_k
            media_k += delta / j
            m2_k += delta * (fx - media_k)
        end

        suma_medias_estratos += media_k
        if samples_per_stratum > 1
            suma_var_estratos += m2_k / (samples_per_stratum - 1)
        end
    end

    integral = h * suma_medias_estratos
    # Varianza del estimador estratificado: ∑ h^2 * (σ_k^2 / N_k)
    std_error = samples_per_stratum > 1 ? h * sqrt(suma_var_estratos / samples_per_stratum) : NaN

    return (integral = integral, std_error = std_error)
end
```

---

### 9.3. Aceptación-Rechazo (con soporte exacto $[a, b]$ y propuesta uniforme)

Este patrón resuelve limpiamente el caso de $\tilde{\rho}(x) = 1 - 1.2x^3$ en $[0, 0.8]$ usando su soporte exacto sin truncamientos erróneos:

```julia
using Random

"""
    sample_rejection_uniform(rho_unnorm, a::Float64, b::Float64, M_unnorm::Float64, n_proposals::Int; seed::Int=1)

Genera muestras exactas distribuidas según rho(x) ∝ rho_unnorm(x) en [a, b]
usando como propuesta una Uniforme(a, b), cuya densidad es g(x) = 1 / (b - a).
Requiere que M_unnorm * (1/(b-a)) >= rho_unnorm(x) para todo x en [a, b]
(equivale a elegir M_unnorm >= (b - a) * max(rho_unnorm)).
"""
function sample_rejection_uniform(
    rho_unnorm, a::Float64, b::Float64, M_unnorm::Float64, n_proposals::Int; seed::Int=1
)
    rng = Xoshiro(seed)
    ancho = b - a
    g_val = 1.0 / ancho
    techo = M_unnorm * g_val

    samples = Float64[] # Tipo concreto Float64[]! Nunca []
    sizehint!(samples, n_proposals) # Evita realocaciones continuas al hacer push!

    for _ in 1:n_proposals
        x0 = a + ancho * rand(rng)
        rho_val = Float64(rho_unnorm(x0))

        # Verificación defensiva de envolvente válida
        @assert 0.0 <= rho_val <= techo "Envolvente M insuficiente en x = $x0 (rho=$rho_val > techo=$techo)"

        prob = rho_val / techo
        if rand(rng) <= prob
            push!(samples, x0)
        end
    end

    efficiency = length(samples) / n_proposals
    return (samples = samples, accepted = length(samples), efficiency = efficiency)
end
```

> **Consejo:** Si tu objetivo final no es guardar los puntos `samples` sino calcular $\int_a^b f(x)\,dx$, mira el bloque siguiente (**Importance Sampling**): es mucho más directo y eficiente!

---

### 9.4. *Importance Sampling* (Integral pura y SNIS para observables físicos)

```julia
using Random

"""
    mc_importance_integral(f, sample_g, pdf_g, n_samples::Int; seed::Int=1)

Estima I = ∫ f(x) dx mediante Importance Sampling: E_g[f(X) / g(X)].
- `sample_g(rng)`: función que devuelve una muestra x ~ g.
- `pdf_g(x)`: función que evalúa la densidad NORMALIZADA g(x).
"""
function mc_importance_integral(f, sample_g, pdf_g, n_samples::Int; seed::Int=1)
    @assert n_samples > 1
    rng = Xoshiro(seed)
    media_w = 0.0
    m2 = 0.0

    for i in 1:n_samples
        x = Float64(sample_g(rng))
        gx = Float64(pdf_g(x))
        @assert gx > 0.0 "La densidad propuesta g(x) debe ser estrictamente positiva en x = $x"

        val = Float64(f(x)) / gx
        delta = val - media_w
        media_w += delta / i
        m2 += delta * (val - media_w)
    end

    std_error = sqrt((m2 / (n_samples - 1)) / n_samples)
    return (integral = media_w, std_error = std_error)
end

"""
    mc_snis_expectation(observable, rho_unnorm, sample_g, pdf_g, n_samples::Int; seed::Int=1)

Self-Normalized Importance Sampling (SNIS): estima ⟨O⟩_ρ = ∫ O(x) ρ(x) dx
cuando ρ(x) = rho_unnorm(x) / Z tiene constante de normalización Z desconocida.
"""
function mc_snis_expectation(observable, rho_unnorm, sample_g, pdf_g, n_samples::Int; seed::Int=1)
    rng = Xoshiro(seed)
    suma_pesos = 0.0
    suma_ponderada = 0.0

    for _ in 1:n_samples
        x = Float64(sample_g(rng))
        w = Float64(rho_unnorm(x)) / Float64(pdf_g(x))
        suma_pesos += w
        suma_ponderada += w * Float64(observable(x))
    end

    return suma_ponderada / suma_pesos
end
```

---

### 9.5. MCMC: Metropolis-Hastings 1D (Potencial físico de doble pozo)

Ejemplo realista de física estadística: muestrear la distribución de Boltzmann $\pi(x) \propto e^{-\beta V(x)}$ de una partícula en un potencial de doble pozo $V(x) = x^4 - 2x^2$, trabajando en logaritmos para máxima estabilidad numérica:

```julia
using Random

"""
    metropolis_1d(log_pi, x0::Float64, step_std::Float64, n_steps::Int, n_burn::Int; seed::Int=1)

Ejecuta Metropolis simétrico con propuesta Gaussiana N(x, step_std^2).
- `log_pi(x)`: logaritmo de la densidad sin normalizar (ej. -beta * V(x)).
- `n_burn`: pasos iniciales de termalización que se descartan.
Devuelve las muestras tras burn-in y la tasa de aceptación.
"""
function metropolis_1d(
    log_pi, x0::Float64, step_std::Float64, n_steps::Int, n_burn::Int; seed::Int=1
)
    @assert n_steps > n_burn >= 0
    rng = Xoshiro(seed)

    x_curr = x0
    log_pi_curr = Float64(log_pi(x_curr))

    n_keep = n_steps - n_burn
    chain = Vector{Float64}(undef, n_keep) # Prealocación exacta sin push!
    accepted = 0

    for t in 1:n_steps
        # Propuesta simétrica: q(x -> y) = q(y -> x)
        x_prop = x_curr + step_std * randn(rng)
        log_pi_prop = Float64(log_pi(x_prop))

        # Razón logarítmica de aceptación: log(π(y)/π(x))
        log_alpha = log_pi_prop - log_pi_curr

        if log_alpha >= 0.0 || log(1.0 - rand(rng)) <= log_alpha
            x_curr = x_prop
            log_pi_curr = log_pi_prop
            accepted += 1
        end

        # Guardar solo después de termalizar; si se rechaza, se repite x_curr!
        if t > n_burn
            @inbounds chain[t - n_burn] = x_curr
        end
    end

    return (chain = chain, acceptance_rate = accepted / n_steps)
end
```

---

### 9.6. Paralelización multihilo segura (`Threads.@threads`)

Patrón de referencia para repetir `reps` experimentos Monte Carlo independientes en paralelo sin *race conditions* ni asignaciones compartidas:

```julia
using Random, Statistics

"""
    mc_parallel_experiments(f, a::Float64, b::Float64, n_samples::Int, reps::Int; seed::Int=1)

Ejecuta `reps` estimaciones Monte Carlo independientes en paralelo sobre múltiples hilos.
Cada realización crea su propio `Xoshiro(seed + i)` y escribe en su índice único `buffer[i]`.
Devuelve `(buffer, media, std_empirica)`.
"""
function mc_parallel_experiments(f, a::Float64, b::Float64, n_samples::Int, reps::Int; seed::Int=1)
    @assert b > a && n_samples >= 1 && reps >= 2
    ancho = b - a
    buffer = zeros(Float64, reps) # Prealocación thread-safe

    Threads.@threads for i in 1:reps
        rng_local = Xoshiro(seed + i) # RNG privado por iteración
        cum = 0.0
        for _ in 1:n_samples
            x = a + ancho * rand(rng_local)
            cum += Float64(f(x))
        end
        buffer[i] = ancho * (cum / n_samples)
    end

    return (buffer = buffer, media = mean(buffer), std_empirica = std(buffer))
end
```

---

## 10. Tabla de Decisión y Checklist Anti-Errores

### ¿Qué método elijo para mi problema?

1. **¿Tu integral es 1D–3D y suave?** $\rightarrow$ Usa cuadratura determinista (Simpson / Gauss-Legendre) o **Muestreo Estratificado**.
2. **¿Tu integral es de dimensión media/alta ($d \ge 4$) sobre una caja acotada?** $\rightarrow$ Usa **Monte Carlo Crudo**, o **Importance Sampling** si la función tiene picos concentrados.
3. **¿Necesitas muestras exactas i.i.d. de $\rho(x)$ en 1D–2D y sabes invertir $F(x)$?** $\rightarrow$ **Transformada Inversa** (100% de eficiencia).
4. **¿Necesitas muestras exactas i.i.d. de $\rho(x)$ en baja dimensión y no puedes invertir $F(x)$?** $\rightarrow$ **Aceptación-Rechazo** con una envolvente $M g(x) \ge \rho(x)$ de soporte coincidente.
5. **¿Quieres simular un sistema de Física Estadística en alta dimensión o sin conocer $Z$?** $\rightarrow$ **MCMC (Metropolis-Hastings)** calibrando $\delta$ para aceptar $\sim 25\%\text{--}50\%$ de pasos y descartando el *burn-in*.

### Checklist antes de ejecutar cualquier código Monte Carlo en Julia

* [ ] **¿Está todo el cálculo metido dentro de una `function`?** (Nunca ejecutes bucles pesados en el ámbito global).
* [ ] **¿El `rng = Xoshiro(seed)` se instancia FUERA del bucle?**
* [ ] **¿Todas las variables y acumuladores mantienen su tipo?** (`cum = 0.0`, `Float64[]`, nada de `[]`).
* [ ] **¿He evitado poner `@inline` delante de bucles `for`?**
* [ ] **¿Coinciden el soporte de $g(x)$ y el del dominio de integración?** (Si usas una propuesta en $[a, b]$, no dividas por la densidad de una Gaussiana sin truncar).
* [ ] **¿Está normalizada la densidad $g(x)$ en *Importance Sampling*?** ($\int_\Omega g(x)\,dx = 1$, y con colas más pesadas que $f$).
* [ ] **En *Importance Sampling*, ¿estoy promediando $f(x_i)/g(x_i)$ sin volver a multiplicar por $(b - a)$?** (Solo se multiplica por $(b - a)$ en Monte Carlo crudo, donde $g(x) = 1/(b - a)$!).
* [ ] **En `Threads.@threads`, ¿tiene cada hilo su propio `rng` y su propio índice de escritura sin `push!` compartido?**
