# Contexto del Proyecto: Resúmenes y Estudio de Ciencia y Análisis de Datos (UNED)

Este archivo centraliza el estado, metodología, correspondencia bibliográfica y reglas del proyecto para que cualquier sesión de trabajo pueda retomar las tareas de inmediato.

---

## 1. Información General del Proyecto y Entorno

* **Titulación / Ámbito:** Máster Universitario en Ciencia y Análisis de Datos (UNED).
* **Directorio de trabajo:** `/home/samuel/Work/Máster UNED/Ciencia y Analisis de Datos`
* **Documento de guía docente:** [`Indicaciones Asignatura.md`](file:///home/samuel/Work/M%C3%A1ster%20UNED/Ciencia%20y%20Analisis%20de%20Datos/Indicaciones%20Asignatura.md)
* **Objetivo:** Elaborar resúmenes universitarios de alto nivel, condensados pero rigurosos, de los libros de texto base, y apoyar el desarrollo de las prácticas y scripts asociados (MATLAB / R / Python).

---

## 2. Textos Base de la Asignatura y Archivos de Resumen

La asignatura articula su programa teórico en torno a dos manuales de referencia (presentes en formato PDF en la raíz del espacio de trabajo):

| Referencia en Asignatura | Archivo PDF original | Título de la obra y Autoría | Archivo de Resumen `.md` |
| :--- | :--- | :--- | :--- |
| **Texto Base 1** | [`data-science-and-big-data-analy-nieizv_book.pdf`](file:///home/samuel/Work/M%C3%A1ster%20UNED/Ciencia%20y%20Analisis%20de%20Datos/data-science-and-big-data-analy-nieizv_book.pdf) | *Data Science and Big Data Analytics: Discovering, Analyzing, Visualizing and Presenting Data*<br>EMC Education Services (Wiley) | [`Data Science and Big Data Analytics.md`](file:///home/samuel/Work/M%C3%A1ster%20UNED/Ciencia%20y%20Analisis%20de%20Datos/Data%20Science%20and%20Big%20Data%20Analytics.md) |
| **Texto Base 2** | [`TheArtOfStatistics.pdf`](file:///home/samuel/Work/M%C3%A1ster%20UNED/Ciencia%20y%20Analisis%20de%20Datos/TheArtOfStatistics.pdf) | *The Art of Statistics: How to Learn from Data* (*El arte de la estadística: Cómo aprender de los datos*)<br>David Spiegelhalter | [`The Art of Statistics.md`](file:///home/samuel/Work/M%C3%A1ster%20UNED/Ciencia%20y%20Analisis%20de%20Datos/The%20Art%20of%20Statistics.md) |

---

## 3. Estado Actual de los Resúmenes

* **Capítulos completados:**
  * **Texto Base 1:** Capítulos 1 y 2 incorporados íntegramente en [`Data Science and Big Data Analytics.md`](file:///home/samuel/Work/M%C3%A1ster%20UNED/Ciencia%20y%20Analisis%20de%20Datos/Data%20Science%20and%20Big%20Data%20Analytics.md).
  * **Texto Base 2:** Capítulos 1 y 2 incorporados íntegramente en [`The Art of Statistics.md`](file:///home/samuel/Work/M%C3%A1ster%20UNED/Ciencia%20y%20Analisis%20de%20Datos/The%20Art%20of%20Statistics.md).
* **Capítulos pendientes:** Pendientes de solicitud explícita del usuario.

---

## 4. Directrices Obligatorias de Elaboración de Resúmenes

1. **Control de alcance estricto:**  
   * **No resumir más capítulos de los que el usuario indique explícitamente.** Si el usuario solicita los capítulos 3 y 4, ceñirse exclusivamente a ellos, ni más ni menos.
2. **Síntesis sin pérdida de información (Alta densidad de conocimiento):**  
   * Los libros originales contienen gran cantidad de prosa narrativa accesoria. El resumen debe **condensar radicalmente la extensión**, pero sin omitir ninguna definición formal, procedimiento metodológico, ecuación, comparativa o lección de caso.
3. **Explicación conceptual completa:**  
   * Ningún concepto técnico se da por sabido. Todo concepto relevante debe explicarse de forma clara, rigurosa y comprensible para nivel de máster universitario.
4. **Estructura visual y anotaciones críticas:**  
   * Usar esquemas, tablas comparativas y fórmulas en LaTeX (`$...$` o `$$...$$`).
   * Incluir bloques destacados para aspectos cruciales:
     * `> [!NOTE]` para aclaraciones conceptuales y matices teóricos.
     * `> [!IMPORTANT]` para principios metodológicos clave.
     * `> [!WARNING]` o `> [!CAUTION]` para trampas habituales, falacias y riesgos de mala interpretación.
5. **Resaltado de conceptos clave:**
   * Utilizar la sintaxis de resaltado `==término clave==` (o `==**término clave**==`) para conceptos nucleares, definiciones y métricas fundamentales, asegurando una lectura rápida e integración óptima con Obsidian.
6. **Idioma:** Redactar en español formal y académico.
7. **Mantenimiento incremental:**  
   * Los nuevos capítulos se añaden al final del archivo `.md` respectivo del libro, manteniendo una estructura homogénea de cabeceras (`## Capítulo X: ...`, `### Subapartados`).

---

## 5. Hoja de Ruta de la Asignatura (Referencia de `Indicaciones Asignatura.md`)

* **Semana 1: Introducción a la Ciencia de Datos** *(Lecturas teóricas completadas)*
  * Capítulos 1 y 2 del Texto Base 1.
  * Capítulos 1 y 2 del Texto Base 2.
  * *Parte práctica:* Manejo básico de MATLAB, matrices, distribuciones de probabilidad y curso *MATLAB Onramp*. Archivo asociado en el workspace: [`Semana1_CyAD_Intro_Matlab.m`](file:///home/samuel/Work/M%C3%A1ster%20UNED/Ciencia%20y%20Analisis%20de%20Datos/Semana1_CyAD_Intro_Matlab.m).
* **Semana 2: Métodos Estadísticos** *(Próxima etapa teórica esperada)*
  * **Texto Base 1:** Capítulo 3 (*Review of Basic Data Analytic Methods Using R* / Métodos analíticos y estadísticos exploratorios).
  * **Texto Base 2:** Capítulos 3 y 4 (Cap. 3: *Poblaciones y mediciones*; Cap. 4: *¿Qué causa qué? / Causalidad*).
  * *Parte práctica:* Representación de datos, contraste de hipótesis, errores Tipo I y II, importación de datos en MATLAB y curso *Statistics Onramp*.
* **Semana 3: Regresión Lineal**
  * **Texto Base 1:** Capítulo 6.1 (*Linear Regression*).
  * **Texto Base 2:** Capítulos 5 y 6 (Cap. 5: *Modelizar relaciones mediante la regresión*; Cap. 6: *Algoritmos, analítica y predicción*).
  * *Parte práctica:* Ajuste de curvas con MATLAB, análisis de componentes principales (PCA), matrices de covarianza y coeficientes de correlación.

---

## 6. Instrucciones para Retomar la Sesión

Cuando el usuario inicie una nueva interacción o pida avanzar:
1. Consultar este archivo [`context.md`](file:///home/samuel/Work/M%C3%A1ster%20UNED/Ciencia%20y%20Analisis%20de%20Datos/context.md) para verificar el punto exacto de avance.
2. Identificar qué capítulos o tareas prácticas solicita el usuario.
3. Actualizar el archivo de resumen correspondiente ([`The Art of Statistics.md`](file:///home/samuel/Work/M%C3%A1ster%20UNED/Ciencia%20y%20Analisis%20de%20Datos/The%20Art%20of%20Statistics.md) o [`Data Science and Big Data Analytics.md`](file:///home/samuel/Work/M%C3%A1ster%20UNED/Ciencia%20y%20Analisis%20de%20Datos/Data%20Science%20and%20Big%20Data%20Analytics.md)) siguiendo estrictamente las directrices del apartado 4.
4. Actualizar el estado de este archivo [`context.md`](file:///home/samuel/Work/M%C3%A1ster%20UNED/Ciencia%20y%20Analisis%20de%20Datos/context.md) para reflejar los nuevos capítulos completados.
