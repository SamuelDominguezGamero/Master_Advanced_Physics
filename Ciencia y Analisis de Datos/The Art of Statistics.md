# The Art of Statistics: How to Learn from Data
**Autor:** David Spiegelhalter  
**Edición base:** *El arte de la estadística: Cómo aprender de los datos*  
**Nivel:** Máster Universitario en Ciencia y Análisis de Datos (UNED)

---

## Capítulo 1: Cómo entender las cosas como proporciones: datos categóricos y porcentajes

El punto de partida del análisis de datos consiste en clasificar, contar y expresar resultados mediante proporciones. Aunque calcular porcentajes parece una operación matemática elemental, el proceso de transformar la realidad en datos cuantitativos y comunicar los resultados de forma rigurosa encierra complejidades metodológicas y psicológicas críticas.

---

### 1. Variables y datos categóricos

Una ==**variable**== es cualquier característica, número o cantidad que puede adoptar distintos valores en diferentes circunstancias o individuos. En función de su escala de medida, se clasifican en:

*   ==**Variables binarias (dicotómicas):**== Solo pueden tomar dos estados mutuamente excluyentes (ej. *vivo / muerto*, *sí / no*, *positivo / negativo*). Suelen codificarse numéricamente como $1$ y $0$, y se resumen de forma directa como recuentos y porcentajes.
*   ==**Variables categóricas nominales (no ordinales):**== Registran categorías cualitativas que carecen de un orden jerárquico natural (ej. hospital donde se interviene un paciente, país de residencia, grupo sanguíneo).
*   ==**Variables categóricas ordinales:**== Presentan un orden intrínseco o escala graduada, aunque la distancia matemática entre categorías consecutivas no sea homogénea ni medible (ej. rangos militares, niveles de dolor de 1 a 5, gravedad de una patología).
*   ==**Variables cuantitativas continuas agrupadas:**== Magnitudes numéricas continuas que se dividen artificialmente en intervalos o clases para facilitar su análisis como categorías (ej. clasificar el ==Índice de Masa Corporal== [$\text{IMC} = \frac{\text{peso (kg)}}{\text{estatura (m)}^2}$] en: *bajo peso* [$< 18{,}5$], *peso normal* [$18{,}5 - 25$], *sobrepeso* [$25 - 30$], *obesidad* [$30 - 35$] y *obesidad mórbida* [$> 35$]).

---

### 2. El reto operativo de la recogida y recuento de datos: El caso Bristol

A través de la investigación sobre la mortalidad en cirugía cardiaca pediátrica en el *Bristol Royal Infirmary* (1984–1995), Spiegelhalter ilustra que «contar» requiere ==**definiciones operativas estrictas**==:

1.  ==**Definición de entidades y eventos:**== Es indispensable acotar operativamente los términos de estudio. En Bristol se definió «niño» como menor de 16 años, «cirugía» como intervención a corazón abierto con baipás cardiopulmonar (agrupando múltiples cirugías por admisión en un único evento), y «fallecimiento» como muerte dentro de los 30 días posteriores a la operación (independientemente del lugar o causa aparente).
2.  ==**Calidad y discrepancia de las fuentes:**== Los registros administrativos rutinarios (como las Estadísticas de Episodios Hospitalarios, codificadas por personal administrativo) arrojaron un $14\,\%$ de mortalidad (62 muertes en 505 cirugías), mientras que el registro profesional de cirujanos consignó un $13\,\%$ (71 muertes en 653 operaciones). Ninguna fuente representa una «verdad objetiva» absoluta; todas reflejan decisiones de diseño, criterios de codificación y posibles errores de captura.
3.  ==**Muertes esperadas vs. «Muertes de más»:**== Si los pacientes de Bristol hubiesen experimentado el riesgo medio del Reino Unido ajustado a la época, se habrían esperado 32 muertes en lugar de las 62 registradas (30 muertes adicionales atribuibles al exceso de riesgo).

> [!CAUTION]
> **Peligro del concepto de «muertes evitables»:** En estadística comparativa, la mitad de los centros estarán por puro azar matemático por encima de la media de mortalidad esperada. Denominar ==«muertes evitables»== a las muertes por encima de la media es un error conceptual grave: la ==variación aleatoria natural== siempre produce desviaciones respecto al valor esperado.

---

### 3. La comunicación de proporciones y el sesgo de encuadre (*Framing*)

La percepción de los datos no depende únicamente de su valor numérico, sino del formato cognitivo en el que se presentan:

*   ==**Encuadre positivo vs. negativo (*Framing effect*):**==
    *   *Positivo:* «Tasa de supervivencia del 95 %».
    *   *Negativo:* «Tasa de mortalidad del 5 %».  
    Ambos expresan idéntica realidad matemática, pero el ==encuadre negativo== activa una percepción de riesgo y gravedad notablemente más elevada.
*   ==**Porcentajes vs. Frecuencias naturales absolutas:**== Convertir un porcentaje bajo en cifras absolutas sobre una población grande magnifica la alarma social. Por ejemplo, la afirmación «el 99 % de los jóvenes no son violentos» parece tranquilizadora; sin embargo, al reformularla en negativo y en números absolutos sobre una población de 1.000.000 de jóvenes, implica que «hay 10.000 jóvenes extremadamente violentos en las calles», produciendo el efecto opuesto.
*   ==**Tablas de clasificación (*League Tables*):**== Ordenar instituciones (colegios, hospitales) de forma jerárquica por sus tasas brutas sin ajustar por la gravedad del caso (==*case-mix*==) ni por variación aleatoria castiga injustamente a los centros de referencia que atienden los casos más complejos.

---

### 4. Métodos para comparar riesgos: Conceptos clave

Para ilustrar la comparación entre dos grupos respecto a una variable binaria (por ejemplo, el riesgo de desarrollar cáncer colorrectal según el consumo diario de beicon / carnes procesadas), se formalizan las siguientes métricas:

| Métrica | Definición conceptual | Ejemplo: Sin beicon vs. Con beicon |
| :--- | :--- | :--- |
| ==**Riesgo Absoluto ($RA$)**== | Probabilidad real directa de que el evento ocurra en un grupo específico ($p = \frac{\text{casos}}{\text{total}}$). | Grupo control: $RA_0 = 6\,\%$ ($0{,}06$).<br>Grupo expuesto: $RA_1 = 7\,\%$ ($0{,}07$). |
| ==**Diferencia de Riesgo Absoluto ($DRA$)**== | Reducción o incremento absoluto del riesgo: $DRA = |RA_1 - RA_0|$. | $7\,\% - 6\,\% = 1\,\%$ (o $0{,}01$). |
| ==**Riesgo Relativo ($RR$)**== | Razón entre el riesgo del grupo expuesto y el del grupo de referencia: $RR = \frac{RA_1}{RA_0}$. | $RR = \frac{0{,}07}{0{,}06} \approx 1{,}17$ o $1{,}18$ (incremento relativo del $18\,\% - 20\,\%$). |
| ==**Número Necesario a Tratar / Exponer ($NNT$ / $NNE$)**== | Número de individuos que deben someterse a la intervención o exposición para observar un evento adicional (o prevenirlo): $NNT = \frac{1}{DRA}$. | $NNT = \frac{1}{0{,}01} = 100$ personas comiendo beicon diariamente toda su vida para observar 1 caso adicional. |
| ==**Odds (Momio o Posibilidad)**== | Razón entre la probabilidad de que el evento ocurra ($p$) y la de que no ocurra ($1-p$): $\text{Odds} = \frac{p}{1-p}$. | Control: $\frac{6}{94} \approx 0{,}0638$.<br>Expuesto: $\frac{7}{93} \approx 0{,}0753$. |
| ==**Odds Ratio ($OR$)**== | Razón de ventajas o cociente entre dos odds: $OR = \frac{\text{Odds}_1}{\text{Odds}_0} = \frac{p_1/(1-p_1)}{p_0/(1-p_0)}$. | $OR = \frac{0{,}0753}{0{,}0638} \approx 1{,}18$. |

#### La trampa de la Odds Ratio y los riesgos comunes: Caso de las estatinas
*   Cuando un evento es ==**raro**== ($p$ es pequeño, ej. $< 5\,\%$), el denominador $1-p \approx 1$, por lo que $\text{Odds} \approx p$ y la ==$OR \approx RR$==.
*   Cuando un evento es ==**común**==, la $OR$ diverge drásticamente del $RR$. Por ejemplo, en un estudio sobre el dolor muscular con estatinas:
    *   Incidencia en control: $85\,\%$ ($p_0 = 0{,}85$).
    *   Incidencia con estatinas: $87\,\%$ ($p_1 = 0{,}87$).
    *   $RR = \frac{0{,}87}{0{,}85} = 1{,}02$ (aumento relativo del riesgo de solo el **$2\,\%$**).
    *   $\text{Odds}_0 = \frac{0{,}85}{0{,}15} = 5{,}67$; $\quad \text{Odds}_1 = \frac{0{,}87}{0{,}13} = 6{,}69$.
    *   $OR = \frac{6{,}69}{5{,}67} = 1{,}18$.
*   Un titular de prensa que confunda la $OR$ con el $RR$ anunciará un alarmante «incremento del riesgo del 20 %», cuando en la práctica se trata de un desplazamiento mínimo del $85\,\%$ al $87\,\%$.

> [!IMPORTANT]
> **Regla de comunicación científica:**
> 1. El ==**riesgo relativo ($RR$)**== por sí solo distorsiona la magnitud del efecto porque oculta el ==riesgo basal==.
> 2. Las ==**frecuencias esperadas**== sobre un denominador fijo (ej. «de cada 100 personas...») y las ==**matrices de iconos (*icon arrays*)**== son las herramientas cognitivas más eficaces e imparciales para comunicar riesgos a especialistas y público general.

---

### 5. Criterios de visualización para proporciones y variables categóricas

1.  ==**Línea base y ejes truncados:**==
    *   Truncar el eje vertical u horizontal (iniciando, por ejemplo, en $85\,\%$ o $95\,\%$ en lugar de $0\,\%$) magnifica desproporcionadamente pequeñas fluctuaciones aleatorias, engañando al lector.
    *   Sin embargo, situar la base en $0\,\%$ en fenómenos de alta supervivencia (como $98\,\%$ vs. $96\,\%$) puede aplanar las barras haciendo invisibles variaciones relevantes.
    *   *Solución metodológica (Alberto Cairo):* Emplear una ==**línea base lógica y justificada**== (por ejemplo, el estándar histórico o el umbral de alerta) y, cuando no se inicie en cero, sustituir las barras por ==**diagramas de puntos**==.
2.  ==**Inconveniencia de los gráficos circulares (*Pie Charts*):**==
    *   El sistema visual humano juzga con poca precisión las diferencias angulares y las áreas relativas en comparación con las longitudes en una escala común.
    *   El efecto 3D añade una distorsión sistemática adicional (hace que las porciones frontales parezcan artificialmente mayores).
    *   *Alternativa recomendada:* ==**Diagramas de barras horizontales**==, ordenados de forma descendente por valor o siguiendo un criterio conceptual lógico.

---
---

## Capítulo 2: Cómo resumir y comunicar los números. Muchísimos números

Cuando se trabaja con variables cuantitativas (continuas o discretas), la acumulación masiva de valores numéricos individuales desborda la cognición. El análisis exploratorio exige sintetizar esta información mediante medidas resumen (tendencia central y dispersión) y representaciones gráficas que revelen la estructura subyacente de la distribución sin desvirtuar la realidad.

---

### 1. Distribución de frecuencias y la «sabiduría de las multitudes»

*   ==**Concepto de Distribución:**== Es el patrón formal de variabilidad de los datos; describe los valores posibles que adopta una variable y la frecuencia con la que se observa cada uno.
*   ==**Experimento del peso del buey (Francis Galton, 1907) y tarro de gominolas (Spiegelhalter, 915 participantes):**==
    *   Galton observó que la mediana de las estimaciones del público en una feria rural sobre el peso de un buey faenado coincidió casi de forma exacta con el peso real (==sabiduría de las multitudes==).
    *   En el experimento de las gominolas (cifra real: 1.616), la estimación individual media estuvo muy alejada de la realidad, pero la ==mediana colectiva== (1.775) se situó más cerca del valor verdadero que el 90 % de los concursantes individuales.
*   ==**Asimetría (*Skewness*):**==
    *   La mayoría de las distribuciones empíricas en el mundo real presentan una ==**asimetría a la derecha (sesgo positivo)**==, caracterizada por una cola larga hacia los valores altos generada por unas pocas observaciones muy elevadas.
    *   En distribuciones asimétricas, las medidas de resumen estándar que asumen simetría fallan rotundamente.

---

### 2. Representación visual de distribuciones numéricas

*   ==**Diagrama de cinta (*Strip plot*) / Gráfico de puntos:**== Muestra cada observación como un punto individual a lo largo de un eje. Es ideal para conjuntos de datos pequeños, pues preserva todos los puntos y revela huecos o agrupaciones locales.
*   ==**Histograma:**== Divide el rango de la variable en intervalos contiguos (==*bins*==) y traza rectángulos cuya área/altura es proporcional a la frecuencia de observaciones.
    *   *Sensibilidad crítica:* Cambiar el ancho del intervalo altera profundamente la forma percibida de la distribución.
*   ==**Diagrama de caja y bigotes (*Box Plot*, John Tukey):**== Estructura gráfica sintética basada en los ==**estadísticos de orden**==:
    *   *Mediana ($Q_2$ / Percentil 50):* Línea divisoria central en la caja.
    *   *Caja central:* Abarca el rango entre el ==**primer cuartil ($Q_1$ / percentil 25)**== y el ==**tercer cuartil ($Q_3$ / percentil 75)**==. Representa el $50\,\%$ central de la población.
    *   *Rango Intercuartílico ($RIC$ o $IQR$):* Longitud de la caja: $RIC = Q_3 - Q_1$.
    *   *Bigotes:* Se extienden hasta el valor observado más lejano que caiga dentro de $1{,}5 \times RIC$ desde los bordes de la caja.
    *   ==*Valores atípicos (*Outliers*):*== Puntos que se encuentran a una distancia superior a $1{,}5 \times RIC$ de $Q_1$ o $Q_3$. Se dibujan de forma individualizada.

---

### 3. Transformaciones matemáticas de escala: La escala logarítmica

Cuando una variable presenta una distribución con un sesgo positivo extremo que abarca varios órdenes de magnitud (ej. ingresos económicos, poblaciones de países, estimaciones cuantitativas libres), la escala lineal compacta la inmensa mayoría de los datos en un bloque ilegible a la izquierda.

*   ==**Fundamento del logaritmo:**== El logaritmo en base 10 transforma multiplicaciones en sumas: $\log_{10}(10) = 1$, $\log_{10}(100) = 2$, $\log_{10}(1.000) = 3$.
*   ==**Justificación sustantiva:**== La ==**transformación logarítmica**== es la herramienta matemática adecuada cuando los errores o variaciones del proceso son de naturaleza ==**relativa o porcentual (multiplicativa)**== en lugar de absoluta (aditiva).
*   ==**Media geométrica:**== Corresponde a la media aritmética calculada sobre los valores logarítmicos transformada de nuevo a la escala original: $\text{MG} = \sqrt[n]{x_1 \cdot x_2 \cdots x_n} = 10^{\frac{1}{n} \sum \log_{10}(x_i)}$. Tiene la propiedad de no verse afectada por la dirección de una razón (ej. euros por naranja frente a naranjas por euro).

---

### 4. Medidas de tendencia central (Promedios)

Un «promedio» es una cifra única destinada a representar el centro de gravedad de una masa de datos. Existen tres definiciones fundamentales:

1.  ==**Media aritmética ($\bar{x}$):**==
    $$\bar{x} = \frac{1}{n}\sum_{i=1}^{n} x_i$$
    *   *Propiedades:* Utiliza algebraicamente todos los valores del conjunto de datos.
    *   *Inconveniente:* ==**No es robusta**==. Un único valor atípico o una distribución con cola larga arrastran la media hacia el extremo, dejando de reflejar la experiencia del individuo representativo.
2.  ==**Mediana:**== Valor que divide a la muestra ordenada exactamente en dos mitades (percentil 50).
    *   *Propiedades:* ==**Medida robusta**==. Su valor es insensible a magnitudes extremas o errores tipográficos en las colas.
3.  ==**Moda:**== El valor o intervalo que aparece con mayor frecuencia en la muestra.

> [!NOTE]
> **Distinción lingüística esencial:** Es fundamental distinguir entre el ==**«promedio de la variable»**== (media aritmética) y la ==**«variable de la persona promedio»**== (mediana). En variables económicas como la renta salarial o el valor inmobiliario, la media es siempre considerablemente superior a la mediana debido a la reducida élite de ingresos multimillonarios.

---

### 5. Medidas de dispersión y variabilidad

El valor central resulta insuficiente e incompleto si no se acompaña de una métrica que describa la amplitud o concentración de las observaciones:

1.  ==**Rango:**== Diferencia entre el valor máximo y el mínimo: $\text{Rango} = x_{\max} - x_{\min}$. Muy inestable y extremadamente vulnerable a valores atípicos.
2.  ==**Rango Intercuartílico ($RIC$):**== $RIC = Q_3 - Q_1$. Mide la dispersión del 50 % central de los datos. Es ==**robusto**== frente a outliers.
3.  ==**Varianza muestral ($s^2$):**== Promedio de las desviaciones al cuadrado de cada observación respecto a la media:
    $$s^2 = \frac{1}{n-1} \sum_{i=1}^{n} (x_i - \bar{x})^2$$
    *(Se divide por $n-1$ en inferencia muestral para obtener un estimador insesgado de la varianza poblacional).*
4.  ==**Desviación típica o estándar ($s$):**== Raíz cuadrada positiva de la varianza: $s = \sqrt{s^2}$. Devuelve la dispersión a las unidades originales de medida.
    *   ==*Regla empírica (solo válida para distribuciones unimodales y simétricas de tipo normal):*==
        *   Aproximadamente el $68\,\%$ de las observaciones se ubica en el intervalo $\bar{x} \pm 1s$.
        *   Aproximadamente el $95\,\%$ de las observaciones se ubica en $\bar{x} \pm 2s$.
        *   Aproximadamente el $99{,}7\,\%$ de las observaciones se ubica en $\bar{x} \pm 3s$.
    *   *Limitación:* Al basarse en diferencias al cuadrado respecto a la media, ==**no es robusta**==; un solo valor anómalo infla artificialmente la desviación típica.

---

### 6. Análisis bivariado: Relaciones entre variables continuas

Para analizar la asociación entre dos variables cuantitativas continuas ($X$ e $Y$), el análisis gráfico preliminar se realiza mediante un ==**diagrama de dispersión (*Scatter plot*)**==.

#### Coeficiente de correlación de Pearson ($r$)
Formalizado por Karl Pearson (1895), mide la fuerza y dirección de una ==**asociación estrictamente lineal**== entre dos variables continuas:

$$r = \frac{\sum_{i=1}^n (x_i - \bar{x})(y_i - \bar{y})}{\sqrt{\sum_{i=1}^n (x_i - \bar{x})^2}\sqrt{\sum_{i=1}^n (y_i - \bar{y})^2}}$$

*   Oscila en el intervalo cerrado $[-1, 1]$.
*   $r = 1$: Relación lineal perfecta positiva.
*   $r = -1$: Relación lineal perfecta negativa.
*   $r \approx 0$: Ausencia de tendencia lineal.

#### Coeficiente de correlación de rango de Spearman ($\rho$)
Aplica la fórmula de correlación sobre los ==**rangos u órdenes**== de las observaciones en lugar de sobre sus valores numéricos en bruto:
*   Mide ==**relaciones monótonas**== (si al crecer $X$, $Y$ crece de forma constante, aunque no siga una línea recta).
*   Es ==**robusto**== frente a valores atípicos severos.

> [!WARNING]
> **Peligros y axiomas de la correlación:**
> 1.  ==**$r = 0$ no implica ausencia de relación:**== Una relación funcional determinista curvilínea (parabólica o sinusoidal) o patrones geométricos complejos (como el ==*Datasaurus Dozen*== de Alberto Cairo o el cuarteto de Anscombe) pueden arrojar un coeficiente de Pearson idéntico a cero. **Siempre se debe graficar la nube de puntos**.
> 2.  ==**Sensibilidad a puntos de palanca:**== Un único dato atípico aislado puede transformar artificialmente una correlación nula en una correlación estadísticamente alta, o viceversa.
> 3.  ==**Correlación no implica causalidad:**== Una correlación entre el volumen quirúrgico de un hospital y su tasa de supervivencia no prueba que mayor práctica salve más vidas; la causalidad inversa (los centros con mejor reputación atraen a más pacientes) o variables de confusión no observadas pueden generar la asociación observada.

---

### 7. Principios avanzados de comunicación y narración visual de datos

1.  ==**Evolución temporal:**== Las series temporales deben representarse con gráficos de líneas. Cuando se comparan múltiples países o entidades a la vez, los *diagramas spaghetti* pueden volverse ininteligibles; se debe recurrir a gráficos de paneles múltiples (==*small multiples*==) o líneas de chispa (==*sparklines*==).
2.  ==**Visualización dinámica e interactiva:**== Popularizada por Hans Rosling (*Gapminder*), permite incorporar el tiempo y múltiples dimensiones (tamaño de burbuja, color, trayectoria) para desmantelar estereotipos erróneos y revelar patrones de evolución global convergentes.
3.  ==**Regla de oro de la comunicación estadística:**== Antes de seleccionar la codificación gráfica, es imprescindible identificar la audiencia, su nivel de competencia cuantitativa y el objetivo cognitivo del mensaje, evitando que la complejidad estética oculte la estructura de los datos.
