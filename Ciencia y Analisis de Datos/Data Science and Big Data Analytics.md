# Data Science and Big Data Analytics: Discovering, Analyzing, Visualizing and Presenting Data
**Editorial / Autores:** EMC Education Services (Wiley)  
**Nivel:** Máster Universitario en Ciencia y Análisis de Datos (UNED)

---

## Chapter 1: Introduction to Big Data Analytics

La explosión del volumen, velocidad y variedad de la información generada a nivel global ha redefinido los límites de la tecnología tradicional de gestión de datos. Este capítulo introduce las bases del fenómeno ==**Big Data**==, analiza la transición desde el ==*Business Intelligence* (BI)== tradicional hacia la ==Ciencia de Datos (*Data Science*)== y define los roles y arquitecturas que estructuran el nuevo ecosistema analítico.

---

### 1. Visión general de Big Data

Big Data no se define únicamente por el tamaño de los ficheros, sino por la confluencia de tres vectores dimensionales clásicos (a los que la práctica añade el valor y la veracidad):

1.  ==**Volumen:**== Cantidades masivas de datos generadas continuamente (pasando de terabytes a petabytes, exabytes y zettabytes).
2.  ==**Velocidad:**== Rapidez con la que se generan, transmiten y requieren procesarse los datos (a menudo en tiempo real o casi real, como flujos de sensores o transacciones financieras).
3.  ==**Variedad:**== Heterogeneidad de los formatos, procedentes de fuentes diversas que no encajan en tablas relacionales tradicionales.

---

### 2. Clasificación de las estructuras de datos

Desde la perspectiva del procesamiento computacional, los datos se categorizan en cuatro estructuras fundamentales:

```
[Datos Estructurados]     --> RDBMS, OLAP, CSV (Esquema estricto)
[Datos Semiestructurados] --> XML, JSON, logs de servidor (Auto-descriptivos, esquema flexible)
[Datos Cuasi-estructurados]-> Rutas de clickstream web (Formatos erráticos, marcas de tiempo inconsistentes)
[Datos No Estructurados]  --> Texto libre, PDFs, vídeo, audio, imágenes (Sin formato predefinido, ~80-90% del total)
```

1.  ==**Datos estructurados (*Structured data*):**== Poseen un tipo, longitud y formato claramente definidos según un modelo de datos o esquema rígido (==*schema-on-write*==). Constituyen únicamente el 10–20 % del universo digital.  
    *Ejemplos:* Tablas relacionales (==RDBMS==), cubos ==OLAP==, archivos CSV estándar y hojas de cálculo.
2.  ==**Datos semiestructurados (*Semi-structured data*):**== Ficheros de texto con patrones internos predecibles y marcas o etiquetas que permiten su segmentación e interpretación sintáctica, sin requerir un esquema relacional fijo.  
    *Ejemplos:* Archivos ==XML== (definidos mediante esquemas XSD), estructuras ==JSON== y registros de eventos (==*system logs*==).
3.  ==**Datos cuasi-estructurados (*Quasi-structured data*):**== Secuencias textuales con estructuras irregulares, formatos heterogéneos y marcadores discontinuos que exigen herramientas especializadas, reglas complejas de parseo y tiempo de procesamiento para dotarlos de estructura analítica.  
    *Ejemplo paradigmático:* Datos de flujo de clics web (==*clickstream*==), donde el orden temporal de URLs navegadas por un usuario refleja su intención pero presenta marcas de tiempo fluctuantes, parámetros no uniformes y valores perdidos.
4.  ==**Datos no estructurados (*Unstructured data*):**== Información que carece de una estructura interna formal o modelo relacional intrínseco. Representan entre el ==80 % y el 90 %== de los datos existentes en las organizaciones.  
    *Ejemplos:* Documentos de texto libre (informes médicos, contratos legales), correos electrónicos, imágenes satelitales, grabaciones de voz y flujos de vídeo.

---

### 3. Repositorios de datos: Perspectiva del analista

Los entornos corporativos albergan diferentes arquitecturas de almacenamiento con repercusiones críticas para el trabajo analítico:

| Tipo de repositorio | Características arquitectónicas | Impacto / Limitaciones para el Analista |
| :--- | :--- | :--- |
| ==**Hojas de cálculo y *Spreadmarts***== | Archivos locales aislados y bases de datos pequeñas de uso departamental. | Crean islas de información no gobernadas, redundancias y riesgo de inconsistencia de versiones. Dependen de extracciones manuales. |
| ==**Data Warehouse Corporativo (EDW)**== | Repositorios centralizados y estandarizados diseñados para transacciones seguras y reporting corporativo. | Excelente para BI e informes estables, pero rígido: agregar nuevas fuentes exige extensos ciclos de diseño de esquemas (ETL). Restringe la experimentación y el cómputo avanzado. |
| ==***Analytic Sandbox* (Workspace)**== | Espacios de trabajo analítico aislados del entorno productivo, gobernados pero flexibles. | Diseñados expresamente para ciencia de datos. Permiten cargar datos crudos sin esquema previo y combinarlos con datos productivos. |

> [!NOTE]
> **Dimensionamiento del *Analytic Sandbox*:**
> Como regla práctica de ingeniería de datos, el ==sandbox analítico== debe dimensionarse con una capacidad de almacenamiento de ==**5 a 10 veces**== el volumen de los conjuntos de datos originales. Esto permite almacenar datos en bruto, copias intermedias, tablas de ingeniería de atributos y datos derivados generados durante el modelado.

#### El paradigma del Cómputo dentro de la Base de Datos (*In-Database Analytics*)
En arquitecturas tradicionales, los datos se extraen de la base de datos y se transfieren a través de la red hacia la herramienta analítica local (ej. la memoria RAM de la máquina del analista). Con Big Data, mover petabytes por la red colapsa la infraestructura. El nuevo enfoque ==**lleva el algoritmo de cálculo a los datos**== (==*In-Database Analytics*== o procesamiento distribuido ==MapReduce/Hadoop==), ejecutando transformaciones y modelos estadísticos directamente en el clúster o motor de almacenamiento.

---

### 4. Estado de la práctica analítica: Business Intelligence (BI) vs. Data Science

Aunque ambas disciplinas extraen valor del conocimiento empresarial, difieren fundamentalmente en su horizonte temporal, naturaleza de preguntas, fuentes de datos y metodología:

| Dimensión | ==Business Intelligence (BI)== | ==Data Science== |
| :--- | :--- | :--- |
| **Preguntas clave** | *¿Qué ocurrió? ¿Dónde? ¿Con qué frecuencia?* | *¿Por qué ocurrió? ¿Qué sucederá después? ¿Cuál es la decisión óptima?* |
| **Horizonte temporal** | ==Retrospectivo== (pasado y presente reciente). | ==Prospectivo== (predictivo y prescriptivo). |
| **Enfoque analítico** | Descriptivo, orientado a métricas de rendimiento (KPIs) y cuadros de mando. | Exploratorio, predictivo, optimización de escenarios, descubrimiento inductivo. |
| **Tipo de datos** | Altamente estructurados, agregados, procedentes del EDW. | Estructurados, semiestructurados, cuasi-estructurados y no estructurados (datos crudos y masivos). |
| **Tolerancia al cambio** | Esquemas estables y gobernados; datos limpios y validados. | Esquemas dinámicos (==*schema-on-read*==); análisis tolerante a ruido y datos imperfectos. |
| **Herramientas típicas** | SQL, herramientas OLAP, cuadros de mando (Tableau, MicroStrategy, Cognos). | R, Python, SQL avanzado, Hadoop, Spark, librerías de *Machine Learning* (MADlib, Mahout). |

---

### 5. Arquitectura analítica tradicional frente a la arquitectura Big Data

*   ==**Cuellos de botella de la arquitectura tradicional:**==  
    En el flujo clásico: $\text{Fuentes operacionales} \xrightarrow{\text{ETL}} \text{EDW} \xrightarrow{} \text{Data Marts} \xrightarrow{} \text{Informes BI}$.  
    El proceso de ==Extracción, Transformación y Carga (ETL)== exige validar y estructurar exhaustivamente cada campo antes de almacenarlo. Esto genera una ==**alta latencia de incorporación de datos**==, fuerza a los analistas a trabajar con datos agregados o submuestras (lo que sesga o limita modelos predictivos avanzados) y excluye por completo los datos no estructurados.
*   ==**El ecosistema emergente de Big Data:**==  
    Permite absorber datos crudos a alta velocidad en un lago de datos (==*Data Lake*==) o sandbox, postergando la transformación analítica al momento del análisis (==*ELT/ETLT*==).

---

### 6. El ecosistema global de Big Data: Los 4 actores

1.  ==**Dispositivos de datos y la «Sensornet»:**== Sensores físicos, teléfonos inteligentes, dispositivos médicos, vehículos conectados y aplicaciones que registran y generan datos continuamente junto con sus metadatos asociados.
2.  ==**Recolectores de datos:**== Organizaciones (sitios web, grandes superficies, operadoras de telecomunicaciones, administraciones públicas) que capturan interacciones y transacciones en origen.
3.  ==**Agregadores de datos:**== Empresas intermediarias que integran, cruzan, estandarizan y enriquecen conjuntos de datos dispares para crear productos analíticos comerciales (ej. agencias de solvencia crediticia).
4.  ==**Usuarios / Compradores de datos:**== Científicos de datos y unidades de negocio que consumen estos conjuntos de datos masivos para optimizar operaciones, mitigar riesgos financieros o predecir patrones de consumo.

---

### 7. Roles clave en el equipo de Big Data y perfil del Científico de Datos

Para que una iniciativa de analítica avanzada prospere, se requiere un equipo interdisciplinar compuesto por:

*   ==**Talento analítico profundo (Data Scientists / Estadísticos):**== Diseñan algoritmos, modelan matemáticamente y descubren relaciones ocultas.
*   ==**Directivos y gestores informados en datos (*Data-savvy managers*):**== Plantean preguntas de negocio relevantes, comprenden las implicaciones estratégicas de los modelos y asignan presupuestos.
*   ==**Roles de habilitación técnica:**== Administradores de Bases de Datos (DBA), Arquitectos de Datos e Ingenieros de Datos.

#### El perfil competencial del Data Scientist
El científico de datos opera en la intersección de cuatro competencias críticas:
1.  ==**Habilidades cuantitativas:**== Sólida base en probabilidad, estadística inferencial, álgebra lineal y optimización matemática.
2.  ==**Habilidades técnicas y computacionales:**== Programación avanzada (R, Python), dominio de bases de datos relacionales y no relacionales (SQL, NoSQL) y entornos de procesamiento distribuido.
3.  ==**Pensamiento escéptico y método científico:**== Disposición a contrastar empíricamente cualquier asunción, cuestionar la procedencia de los datos y evitar sesgos de confirmación.
4.  ==**Visión de negocio y comunicación (*Storytelling*):**== Comprensión contextual del problema empresarial y capacidad de traducir hallazgos matemáticos abstractos en narrativas claras y decisiones operativas tangibles para audiencias no técnicas.

---
---

## Chapter 2: Data Analytics Lifecycle

A diferencia del desarrollo de software tradicional (que suele seguir modelos lineales tipo cascada o iteraciones cortas de desarrollo funcional), los proyectos de ciencia de datos son fundamentalmente exploratorios, hipotéticos y experimentales. El ==**Data Analytics Lifecycle**== proporciona un marco de gobernanza metódico, compuesto por seis fases cíclicas e interconectadas que aseguran el rigor científico sin estrangular la exploración creativa.

```
       +-------------------------------------------------+
       |                                                 |
       v                                                 |
[ Fase 1: Descubrimiento (Discovery) ]                   |
       |                                                 |
       v                                                 |
[ Fase 2: Preparación de Datos (Data Prep) ]             |
       |                                                 |
       v                                                 |
[ Fase 3: Planificación del Modelo (Model Planning) ]    |  (Flujo iterativo
       |                                                 |   con bucles
       v                                                 |   de retroalimentación)
[ Fase 4: Construcción del Modelo (Model Building) ]     |
       |                                                 |
       v                                                 |
[ Fase 5: Comunicación de Resultados ]                   |
       |                                                 |
       v                                                 |
[ Fase 6: Operacionalización (Operationalize) ]          |
       |                                                 |
       +-------------------------------------------------+
```

---

### 1. Marco conceptual y los 7 roles clave del proyecto

El ciclo sintetiza las mejores prácticas del ==**método científico**==, la metodología estándar de minería de datos ==**CRISP-DM**==, el marco ==**DELTA**== de Tom Davenport y el enfoque de información aplicada ==**AIE** (*Applied Information Economics*)== de Doug Hubbard.

Para coordinar la ejecución del ciclo de vida se articulan **siete roles funcionales**:

1.  ==**Usuario de Negocio (*Business User*):**== Experto en el dominio operativo que define la necesidad, valida la lógica sectorial del problema y consume las soluciones finales.
2.  ==**Patrocinador del Proyecto (*Project Sponsor*):**== Directivo de alto nivel que aporta la visión estratégica general, define los objetivos clave de negocio, financia la iniciativa y evalúa el retorno de la inversión (ROI).
3.  ==**Jefe de Proyecto (*Project Manager*):**== Coordina los recursos humanos y técnicos, gestiona el cronograma, supervisa el cumplimiento de hitos y mitiga retrasos.
4.  ==**Analista de BI (*Business Intelligence Analyst*):**== Aporta conocimiento sobre los datos existentes en los almacenes corporativos (EDW), los indicadores clave actuales (KPIs) y los reportes de rendimiento histórico.
5.  ==**Administrador de Bases de Datos (*DBA*):**== Aprovisiona, configura y asegura el entorno de bases de datos y sandboxes, gestionando permisos y políticas de acceso.
6.  ==**Ingeniero de Datos (*Data Engineer*):**== Especialista en arquitectura de datos y consultas de alto rendimiento; construye los conductos de extracción masiva (==*Big ETL*==), gestiona el flujo de ingesta y acondiciona los datos dentro del clúster o sandbox.
7.  ==**Científico de Datos (*Data Scientist*):**== Responsable metodológico de la formulación de hipótesis, exploración matemática, selección y calibración de algoritmos, evaluación de modelos y comunicación de hallazgos.

---

### 2. Las 6 fases del ciclo de vida analítico

---

#### Fase 1: Descubrimiento (*Discovery*)
En esta fase inicial, el equipo define el alcance del proyecto y evalúa su viabilidad antes de invertir recursos en análisis técnicos masivos.

*   ==**Comprensión del dominio del negocio (*Learning the business domain*):**== Estudio del contexto sectorial, terminología específica y antecedentes de iniciativas previas (identificando por qué funcionaron o fracasaron).
*   ==**Auditoría y evaluación de recursos:**== Inventario de personal disponible, infraestructura de cómputo, volumen temporal asignado y accesibilidad a los datos.
*   ==**Encuadre del problema (*Problem framing*):**== Traducir una necesidad de negocio abstracta (ej. «perdemos clientes») en una formulación formal de analítica avanzada (ej. «predecir la probabilidad de baja voluntaria en los próximos 60 días con un horizonte de intervención anticipada»).
*   ==**Entrevistas a los patrocinadores y definición de criterios de éxito:**== Fijar con precisión las métricas cuantitativas que determinarán si el proyecto es exitoso tanto a nivel analítico (precisión, sensibilidad) como operativo (ahorro económico, reducción de bajas).
*   ==**Desarrollo de Hipótesis Iniciales (*Initial Hypotheses - IH*):**== Formular conjeturas analíticas formales y contrastables que guiarán la investigación inicial.
*   ==**Identificación de fuentes de datos:**== Catálogo de fuentes internas y externas requeridas, considerando su formato, disponibilidad y volumen.

---

#### Fase 2: Preparación de los datos (*Data Preparation*)
Es la fase más intensiva en mano de obra y recursos del ciclo: habitualmente absorbe ==**el 50 % o más del tiempo total de un proyecto de ciencia de datos**==.

*   ==**Aprovisionamiento del *Analytic Sandbox*:**== Configuración del entorno aislado de alto rendimiento.
*   ==**Ejecución de ETLT (*Extract, Transform, Load, Transform*):**==  
    A diferencia del ETL corporativo clásico (donde los datos se limpian rigurosamente *antes* de entrar a la base de datos), en Big Data se aplica preferentemente ==**ELT o ETLT**==:
    *   Se extraen los datos crudos y se cargan directamente en el sandbox (==*Extract & Load*==).
    *   La transformación pesada se realiza mediante el motor distribuido dentro de la propia base de datos (==*Transform*==).
    *   *Razón metodológica:* Conservar una réplica exacta de los datos brutos no alterados permite detectar patrones sutiles o anomalías (ej. señales de fraude en transacciones atípicas) que un filtrado previo de limpieza eliminaría por error.
*   ==**Ingesta masiva (*Big ETL*):**== Uso de tecnologías distribuidas como ==Hadoop== y ==MapReduce== para paralelizar la ingesta masiva de archivos.
*   ==**Familiarización con los datos (*Data exploration & cataloging*):**== Detección de tipos de variables, coherencia de claves foráneas y solidez de metadatos.
*   ==**Acondicionamiento de datos (*Data conditioning*):**==
    *   Limpieza de ruido y resolución de inconsistencias léxicas.
    *   Tratamiento sistemático de ==valores ausentes (NAs / nulos)==: evaluar si se imputan (mediante media, mediana o modelos predictivos) o si se descartan registros.
    *   Normalización y escalado de formatos numéricos y fechas.
    *   Deduplicación de registros.
*   ==**Inspección y visualización preliminar:**== Empleo de histogramas, diagramas de caja y estadísticos descriptivos para auditar la distribución de variables y detectar valores atípicos (==*outliers*==).
*   *Herramientas habituales:* Hadoop, MapReduce, Alpine Miner, OpenRefine, Data Wrangler, PostgreSQL/MADlib, R, Python.

---

#### Fase 3: Planificación del modelo (*Model Planning*)
El equipo determina las técnicas analíticas idóneas, evalúa relaciones multidimensionales y selecciona los subconjuntos de datos sobre los que se construirán los algoritmos.

*   ==**Exploración de datos y selección de variables (*Feature Selection*):**==
    *   Identificación de la variable objetivo o dependiente ($Y$) y de las variables candidatas independientes o predictoras ($X_1, X_2, \dots, X_p$).
    *   Evaluación de la correlación entre variables para descartar redundancias y evitar la ==**multicolinealidad**== (que desestabiliza modelos de regresión lineal o logística).
*   ==**Selección de técnicas de modelado:**==  
    Elección del paradigma algorítmico en función de la naturaleza del problema:
    *   ==*Modelos descriptivos y no supervisados:*== Agrupamiento (==*Clustering*== como K-means), reglas de asociación (==*Market Basket Analysis / Apriori*==).
    *   ==*Modelos predictivos y supervisados:*== Clasificación (árboles de decisión, regresión logística, Naive Bayes), regresión continua o modelos de series temporales.
*   *Herramientas habituales:* R, Python, PL/R, SQL analítico, SAS/STAT, Matlab.

---

#### Fase 4: Construcción del modelo (*Model Building*)
Fase de computación intensiva en la que se ejecutan los algoritmos y se ajustan iterativamente sus parámetros.

*   ==**Partición de datos (*Train / Test / Validation Split*):**==  
    Para evitar el ==**sobreajuste (*overfitting*)**== (situación en la que el modelo memoriza el ruido del conjunto de entrenamiento pero fracasa al predecir casos nuevos), el conjunto de datos se divide aleatoriamente en:
    *   ==*Conjunto de entrenamiento (*Training set*):*== Se utiliza para calcular los parámetros del algoritmo.
    *   ==*Conjunto de prueba (*Testing/Validation set*):*== Datos no vistos por el modelo, reservados exclusivamente para medir su capacidad real de generalización.
*   ==**Ejecución y calibración de modelos:**== Entrenamiento de los algoritmos en el sandbox o infraestructura clúster.
*   ==**Evaluación del rendimiento y diagnóstico:**==
    *   En modelos de clasificación: Evaluación de precisión, sensibilidad (==*recall*==), especificidad, curvas ROC y área bajo la curva (==AUC==).
    *   En modelos de regresión: Coeficiente de determinación ($R^2$), error cuadrático medio (==RMSE==) y análisis de residuos.
*   *Herramientas habituales:* R, Python, Weka, Alpine Miner, MADlib, SAS Enterprise Miner, SPSS.

---

#### Fase 5: Comunicación de resultados (*Communicate Results*)
Tras la validación técnica, el equipo traduce las conclusiones a un lenguaje comprensible y útil para los responsables de la toma de decisiones.

*   ==**Evaluación de hipótesis y criterios de éxito:**== Determinar de forma explícita si las hipótesis iniciales (IH) formuladas en la Fase 1 fueron confirmadas o refutadas empíricamente.
*   ==**Cuantificación del impacto de negocio:**== Traducir las métricas estadísticas a métricas económicas operativas (ej. «la implantación de este modelo reducirá el fraude en 2,3 millones de euros anuales»).
*   ==**Documentación de hallazgos y limitaciones:**== Registrar con total transparencia las condiciones en las que el modelo es válido, las restricciones de los datos, los sesgos identificados y las asunciones no verificables.
*   ==**Doble nivel de entrega comunicativa:**==
    *   ==*Resumen ejecutivo:*== Centrado en valor, recomendaciones estratégicas y riesgos para el comité de dirección.
    *   ==*Informe técnico detallado:*== Código, matrices de confusión y diagnósticos para analistas y equipos de ingeniería.

---

#### Fase 6: Operacionalización (*Operationalize*)
El modelo abandona el entorno controlado del sandbox y se traslada a la operativa diaria de la organización.

*   ==**Proyecto Piloto acotado:**==  
    Antes de una implantación global, se recomienda desplegar el modelo en un entorno de producción controlado (o mediante pruebas A/B en un subconjunto de usuarios). Esto permite auditar la fiabilidad de los conductos de datos en vivo, medir la respuesta real de los usuarios y verificar si el rendimiento del modelo se mantiene estable en producción.
*   ==**Creación de conductos automatizados de puntuación (*Scoring Pipelines*):**== Integración del modelo en las aplicaciones operativas mediante APIs o trabajos batch programados.
*   ==**Monitorización continua y degradación del modelo (*Concept Drift*):**==  
    Los patrones del comportamiento humano y los mercados cambian con el tiempo. El rendimiento de un modelo en producción tiende a degradarse (*drift*). Es imprescindible implementar alertas de rendimiento y protocolos periódicos de reentrenamiento con nuevos datos.
*   ==**Entregables finales:**== Código fuente productivo, documentación de arquitectura, diccionarios de datos actualizados y plan de soporte operativo.

---

### 3. Caso de estudio práctico: GINA (*Global Innovation Network and Analysis*)

Para ilustrar de forma integral las seis fases en una gran corporación tecnológica, el libro documenta el proyecto ==**GINA**== de EMC:

1.  ==**Fase 1 (Descubrimiento):**== La dirección de innovación global de EMC requería un sistema riguroso para monitorizar, medir y catalizar la actividad investigadora e innovadora de sus centros de excelencia distribuidos por el mundo. El equipo definió como hipótesis inicial que el análisis semántico de datos formales (patentes registradas) e informales (conversaciones e ideas en plataformas internas) permitiría identificar a los «conectores clave» (==*key influencers*==) de la innovación que permanecían invisibles en el organigrama oficial.
2.  ==**Fase 2 (Preparación de datos):**== Se construyó un sandbox analítico centralizado integrando bases de datos de recursos humanos, repositorios de patentes y miles de propuestas textuales no estructuradas. Se ejecutaron pipelines ETLT para normalizar nombres y metadatos heterogéneos de diferentes países.
3.  ==**Fase 3 (Planificación del modelo):**== Se seleccionaron algoritmos de Procesamiento de Lenguaje Natural (PLN / Text Mining) para extraer palabras clave y conceptos tecnológicos, complementados con Análisis de Redes Sociales (==SNA - *Social Network Analysis*==) para modelar los grafos de colaboración.
4.  ==**Fase 4 (Construcción del modelo):**== Aplicación de agrupamiento (==*clustering*==) semántico con R y MADlib para descubrir comunidades emergentes de investigación técnica (ej. virtualización, almacenamiento flash, computación en la nube) y mapear los flujos de colaboración transfronteriza.
5.  ==**Fase 5 (Comunicación de resultados):**== Presentación interactiva ante la alta dirección evidenciando cómo determinados centros de excelencia operaban como nodos puente de transferencia tecnológica. Las hipótesis iniciales fueron validadas y se descubrieron clústeres interdisciplinares de innovación no planificados.
6.  ==**Fase 6 (Operacionalización):**== Se implementó una aplicación web interna permanente con cuadros de mando interactivos que permitía a los directores de I+D buscar expertos en tecnologías emergentes en tiempo real y asignar financiación a iniciativas de investigación basadas en datos objetivos.
