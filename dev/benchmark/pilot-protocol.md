# Piloto del benchmark de consenso de dar

Versión del protocolo: 0.2 — enmienda prospectiva del 14 de octubre de 2026.

Estado: **enmienda preparada; nuevo preflight y piloto pendientes de ejecución**.
La versión 0.1 ya se ejecutó en `dev/benchmark/results/pilot-20261014-01` y sus
artefactos se conservan bajo su contrato original. La autorización actual abarca
los cambios de métodos y sus comprobaciones, no lanzar ajustes del preflight o
del piloto. Este documento no concede permiso para ejecutarlos, instalar
dependencias ausentes, descargar datos ni contratar cómputo.

La versión 0.1 se preparó en `devel`, inspeccionando dar 1.9.12, commit
`988d5431b706e53b0f6f3f2274bdfb5f79578f6c`. La inspección comenzó en
`7cb41e01ff3cf338798f4c39f6b8e4112940398d`; las dos actualizaciones recibidas
durante la preparación afectan a automatizaciones, no a estos métodos.
Cada nueva ejecución fijará el commit efectivo, los hashes del código y del
protocolo, el backend realmente utilizado y las versiones de dependencias.

### Enmienda 0.2 e historial

- La versión 0.1 evaluó ANCOM-BC2 básico con `pseudo_sens = FALSE` y LinDA con
  MicrobiomeStat 1.4 antes de la corrección local de ceros adaptativos. Sus
  resultados no se reinterpretan como llamadas robustas ni como LinDA corregido.
- La versión 0.2 activa `pseudo_sens = TRUE` y exige el voto nativo
  `diff_robust_<coeficiente>` de ANCOM-BC2. Para LinDA utiliza el backend local
  corregido que selecciona dar; el benchmark no aplica un parche propio.
- Se mantienen las semillas, los generadores y escenarios, la verdad conocida,
  las métricas, el panel, las reglas de consenso y los criterios A/B/C/D. La
  enmienda no aumenta efectos ni cambia umbrales a partir de los resultados.
- Antes de cualquier nuevo piloto se requiere un preflight nuevo autorizado,
  aprobado y guardado en otro directorio, con identidad coincidente del
  código/protocolo/runtime. No se reutilizan receipts ni ajustes anteriores bajo
  el contrato modificado y no se sobrescriben resultados históricos.

Repetir las mismas semillas permite comparar la corrección metodológica con la
ejecución anterior; no constituye una evaluación confirmatoria independiente.

## 1. Pregunta y alcance

¿Merece la pena evaluar a mayor escala el intercambio entre falsos
descubrimientos y sensibilidad de dos reglas **ya existentes** de consenso
direccional, frente a sus componentes individuales?

El piloto tiene dos objetivos distintos:

1. Comprobar que datos, contrastes, direcciones, llamadas significativas y
   artefactos se pueden evaluar correctamente con un coste asumible.
2. Detectar una ventaja práctica prometedora o un límite metodológico
   suficientemente informativo para justificar una evaluación independiente.

No demostrará control de FDR al 5 %, superioridad universal, nuevos biomarcadores
ni recuperación de cargas microbianas absolutas en pacientes. No pretende
optimizar motores, pesos o umbrales. Un resultado negativo no invalida dar como
software de integración y reproducibilidad.

Alcance: dos grupos independientes, una comparación tratamiento–control y
abundancia. El ajuste por batch se evalúa en dos condiciones. Quedan fuera del
piloto longitudinalidad, prevalencia como desenlace, ceros estructurales
dependientes del grupo, selección de pesos y benchmark completo de los ocho
motores del paquete.

## 2. Panel, referencia y reglas congeladas

Panel de cuatro motores: ANCOM-BC2, LinDA, ALDEx2 y MaAsLin3; una configuración
por motor. No son cuatro replicaciones independientes: utilizan los mismos datos.

Referencia individual fijada de antemano: **LinDA**. Se elige por su compatibilidad
con ambos diseños y porque constituye un comparador individual actual; no porque
se haya observado que gana. Se mostrarán también todos los otros componentes,
sin cambiar la referencia según la condición.

Comparador sencillo: **Wilcoxon sobre CLR**, solo en condiciones sin batch. No
forma parte del consenso. No se usa su normalización composicional por defecto:
su objetivo relativo sería diferente del objetivo principal del piloto.

Reglas principales, siempre sobre el mismo taxón, contraste y dirección:

| Regla | Votos necesarios | Pesos | Uso |
|---|---:|---|---|
| Mayoría estricta | 3 de 4 | Todos 1 | Regla principal para decidir ampliación |
| Unanimidad | 4 de 4 | Todos 1 | Regla conservadora secundaria |

Se ejecutarán los `bake()` y `cool()` reales, excluyendo expresamente el ID de
Wilcoxon si comparte el objeto preparado. No se sustituirán por un voto
reconstruido a partir de p-valores. No se añadirán unión, pesos aprendidos o
umbrales adicionales para buscar una regla ganadora. Un fallo no convierte el
panel en «3 de 3».

Usar `weights = NULL`, que asigna peso uno a cada motor no excluido. Si se
explicitan los pesos para el manifiesto, nombrar solo los cuatro IDs del panel,
no Wilcoxon: excluir un ID y conservarlo en el vector de pesos es inválido.

## 3. Modelo y configuración de los motores

En todas las condiciones: referencia `condition = control`, tratamiento como
numerador, `targets = "condition"`, `contrasts = "simple"`,
`na_action = "error"`, sin variable de tiempo. Fórmula `~ condition`, o
`~ condition + batch` en P5/P6. Condition y batch son factores de dos niveles.

Cada unidad sintética tendrá un ID y una anotación Species únicos. El rango
`Species` no producirá agregaciones adicionales ni mezclará nombres compartidos.
En G2, la agregación de los datos fuente se hará antes de generar muestras,
conservando la identidad de linaje definida por dar.

Perfil explícito, no «todos los valores por defecto»: BH y umbral 0.05, sin
rarefacción ni filtros adicionales de tamaño de efecto. Los controles nativos
de transformación, robustez y ceros especificados abajo son parte del método.

| Motor / ID | Parámetros explícitos |
|---|---|
| ANCOM-BC2 / `ancom_pilot` | `p_adj_method = "BH"`, `prv_cut = 0`, `lib_cut = 0`, `s0_perc = 0.05`, `struc_zero = FALSE`, `neg_lb = FALSE`, `alpha = 0.05`, `n_cl = 1`, `verbose = FALSE`, `rarefy = FALSE`, `engine_args = list(fit = list(pseudo_sens = TRUE))`; voto nativo `diff_robust_<coeficiente>` |
| LinDA / `linda_pilot` | `prev_filter = 0`, `mean_abund_filter = 0`, `max_abund_filter = 0`, `winsorize = TRUE`, `outlier_pct = 0.03`, `adaptive = TRUE`, `zero_handling = "pseudo-count"`, `pseudo_count = 0.5`, `corr_cut = 0.1`, `p_adj_method = "BH"`, `alpha = 0.05`, `n_cpus = 1`, `rarefy = FALSE` |
| ALDEx2 / `aldex_pilot` | `max_significance = 0.05`, `mc.samples = 128`, `denom = "all"`, `rarefy = FALSE`; etapa GLM con `fdr.method = "BH"` |
| MaAsLin3 / `maaslin_pilot` | `min_abundance = 0`, `min_prevalence = 0`, `min_variance = 0`, `normalization = "TSS"`, `transform = "LOG"`, `max_significance = 0.05`, `correction = "BH"`, `standardize = TRUE`, `median_comparison_abundance = TRUE`, `rarefy = FALSE`, `engine_args = list(fit = list(cores = 1), contrast = list(subtract_median = TRUE))` |
| Wilcoxon / `wilcox_pilot` | `norm_method = "clr"`, `max_significance = 0.05`, `p_adj_method = "BH"`, `rarefy = FALSE`, `engine_args = list(test = list(exact = FALSE))` |

Los parámetros restantes conservarán los valores de la versión congelada y se
registrarán en el manifiesto, junto con los argumentos efectivos de cada etapa.
No se copiarán sin revisión los controles relajados de la validación interna.
Las 128 instancias Monte Carlo de ALDEx2 no son 128 réplicas del benchmark.

### Precauciones verificadas en la implementación actual

- MaAsLin3: `median_comparison_abundance = TRUE` cambia el test, pero no centra
  automáticamente el coeficiente devuelto. Por eso se fija
  `contrast = list(subtract_median = TRUE)`, para que el signo usado por dar
  corresponda al efecto respecto a la mediana probada. La etapa de contraste de
  dar utiliza solo abundancia; no se mezclarán sus resultados con prevalencia.
  Fuente: [tutorial oficial, Median comparisons](https://www.bioconductor.org/packages/release/bioc/vignettes/maaslin3/inst/doc/maaslin3_tutorial.html).
- ANCOM-BC2: con `pseudo_sens = TRUE`, dar utiliza
  `diff_robust_<coeficiente>`, que incorpora la sensibilidad a pseudocuentas.
  El preflight exige tanto el argumento efectivo como la igualdad entre la
  clasificación de dar y esa columna nativa. Una columna robusta ausente es un
  fallo de contrato, no una invitación a usar `diff_<coeficiente>` o `q < 0.05`.
- LinDA: MicrobiomeStat 1.4 asigna `"Imputation"` en una rama adaptativa, pero
  comprueba después `"imputation"`. Dar corrige únicamente esa asignación en
  una copia local de la función, manteniendo su entorno nativo, winsorización,
  ajuste del modelo y corrección múltiple. No modifica la dependencia instalada.
  El preflight debe identificar el backend efectivo y verificar la corrección,
  además de conservar argumentos y salidas nativas. `adaptive = TRUE` por sí
  solo no demuestra que se haya ejecutado el backend correcto.

La clasificación `significant` obtenida por `tidy_results()` es la autoridad
para las llamadas individuales y el consenso. Se conservarán también p, q y
efecto, pero no se reconstruirá cada llamada como `q < 0.05`.

## 4. Generadores y verdad conocida

Se usarán **100 taxones** y perfiles independientes por muestra sintética.
Los taxones con señal estarán en la zona intermedia de abundancia, no se
seleccionarán por los resultados ni por supervivencia en los ajustes.

Ordenar los 100 taxones por abundancia basal creciente, desempatar por ID y fijar
las posiciones de señal de condición en
`21, 27, 33, 39, 45, 51, 57, 63, 69, 75`. Cuando sea equilibrada, alternar
signos positivo/negativo. Las señales de batch ocuparán
`22, 28, 34, 40, 46, 52, 58, 64, 70, 76`, también alternando signos.
Estas elecciones no se modificarán después de ver descubrimientos.

### G1: logística-normal y muestreo multinomial

Para muestra i y taxón j:

`log(A_ij) = mu_j + beta_j * condition_i + gamma_j * batch_i + epsilon_ij`.

- `mu_j`: secuencia equiespaciada entre `log(0.02)` y `log(20)`.
- Ruido logarítmico con SD 0.8 y bloques consecutivos de 10 taxones:
  `epsilon_ij = 0.8 * (sqrt(0.3) * z_i,block(j) + sqrt(0.7) * z_ij)`.
  Todos los z son normales estándar independientes; la correlación intrabloque
  es 0.3 antes del cierre composicional.
- `p_ij = A_ij / sum_j(A_ij)` y conteos multinomiales con profundidad `N_i`.
- `N_i = max(500, round(exp(log(D) + 0.5 * z_i)))`, independiente de condición,
  batch y abundancias. D se fija por celda.
- No se borrarán conteos posteriormente para crear dropout: los ceros de este
  brazo son ceros de muestreo. No se afirma evaluar ceros estructurales.

### G2: remuestreo de perfiles empíricos con perturbación sintética

Única fuente del piloto: `metaHIV_phy` ya distribuido en dar. No se descargarán
cohortes ni se usarán sus categorías clínicas para generar el desenlace.

Preparación reproducible de la reserva, antes de cualquier ajuste:

1. Agregar al rango Species mediante la identidad de taxones/linajes de dar.
2. Si hay un identificador inequívoco de donante, conservar una muestra por
   donante mediante orden de ID; si no lo hay, registrar «independencia de
   donantes no verificable». Esta es una reserva de perfiles, no validación
   clínica ni una nueva cohorte independiente.
3. Excluir taxones sin conteos en toda la reserva. Entre los detectados en al
   menos el 10 % de los perfiles, elegir los 100 de mayor prevalencia,
   desempatar por abundancia relativa media en la reserva completa y después por
   ID. Si hay menos de 100,
   detener la preparación y proponer una enmienda; no reducir J silenciosamente.
4. Guardar los IDs, conteos, muestras y checksum de esta reserva. Reordenar los
   taxones por abundancia basal para asignar las posiciones de señal anteriores.

Para cada nueva muestra, extraer con reemplazo un **perfil completo** de la misma
reserva, independientemente del grupo. Suavizar una sola vez:
`p0_ij = (count_ij + 0.5) / sum_j(count_ij + 0.5)`.
Implantar `A_ij = p0_ij * exp(beta_j * condition_i)`, cerrar a proporciones y
muestrear conteos multinomiales con la profundidad definida para G1.

La estructura multivariada basal es empírica, pero la verdad es la perturbación
matemática implantada, **no** una carga absoluta conocida de los pacientes.
Suavizado, reserva finita, elección de taxones y único contexto intestinal
limitan la generalización. Este brazo es una adaptación sencilla de implantación
de señal, no una reproducción exacta de SIMBA. La importancia de contrastar
realismo y generadores está documentada por
[Wirbel et al., 2024](https://link.springer.com/article/10.1186/s13059-024-03390-9).

### Tres tablas de verdad, sin mezclar objetivos

Guardar para cada taxón/contraste:

- `beta_latent`: efecto implantado beta, bajo el modelo sintético y su referencia
  de taxones mayoritariamente invariantes.
- `beta_clr_U = beta_j - mean(beta_U)`.
- `beta_median_U = beta_j - median(beta_U)`.

U es el universo común introducido en los ajustes. Los coeficientes son
logarítmicos naturales; para descubrir y evaluar dirección solo se necesita su
signo. Las escalas numéricas nativas no se promediarán ni compartirán un MAE.

Fijar nulidad como `abs(beta_objetivo) <= 1e-10` en cada tabla de verdad,
exclusivamente como tolerancia numérica de cancelación, no como umbral de efecto
biológico. Una señal es un taxón no nulo según la tabla que se está evaluando.

En las señales equilibradas estas verdades coinciden si el universo conserva el
equilibrio. En P4/P9, 10 efectos `log(4)` implican media `0.1 * log(4)` sobre los
100 taxones y mediana cero: los 90 taxones sin perturbación latente tienen efecto
CLR negativo, no efecto CLR nulo. Se mostrarán las puntuaciones por objetivo,
sin proclamar que todos esos descubrimientos sean errores de ALDEx2/Wilcoxon.
Esas celdas son diagnósticas y no entrarán en el criterio de ganancia principal.

El consenso no recibe un estimando común por el mero hecho de compartir fórmula.
Si el universo cambia o desaparecen señales, recalcular las verdades de referencia
correspondientes, registrar el desplazamiento y mostrar ambos objetivos. No
descartar réplicas ni redefinir retrospectivamente la verdad para favorecer una
regla. Si las condiciones principales dejan de compartir la misma hipótesis,
queda bloqueada la interpretación comparativa hasta una enmienda del protocolo.

## 5. Celdas y tamaño del piloto

20 réplicas independientes por celda. Todos los motores y ambos consensos se
evalúan sobre el mismo dataset de cada réplica.

| Celda | Generador | n por grupo | D | Efecto de condición | Batch / diseño | Finalidad |
|---|---|---:|---:|---|---|---|
| P1 | G1 | 30 | 15 000 | Ninguno | Sin batch | Nulo global |
| P2 | G1 | 30 | 15 000 | 5 positivos + 5 negativos, magnitud `log(2)` | Sin batch | Condición principal |
| P3 | G1 | 15 | 1 500 | Como P2 | Sin batch | Estrés de observabilidad |
| P4 | G1 | 30 | 15 000 | 10 positivos, magnitud `log(4)` | Sin batch | Discordancia de estimandos |
| P5 | G1 | 30 | 15 000 | Ninguno | Batch confundido; ajuste correcto | Nulo condicional |
| P6 | G1 | 30 | 15 000 | Como P2 | Batch confundido; ajuste correcto | Señal con ajuste |
| P7 | G2 | 30 | 15 000 | Ninguno | Sin batch | Nulo con ruido empírico |
| P8 | G2 | 30 | 15 000 | Como P2 | Sin batch | Condición principal empírica |
| P9 | G2 | 30 | 15 000 | Como P4 | Sin batch | Límite bajo otro generador |

En P5/P6, 6/30 controles y 24/30 tratados tendrán batch B; el resto batch A.
Batch causa cinco efectos positivos y cinco negativos de magnitud `log(2)` en
los taxones fijados para batch. No se confunde perfectamente con condición.
En P5 toda la verdad para el contraste de condición es nula **condicional a
batch**, aunque existan asociaciones marginales. Wilcoxon no se ejecutará en
P5/P6 ni se comparará allí como si ajustara por batch.

P3 altera tamaño y profundidad simultáneamente: diagnostica baja observabilidad,
no identifica por separado el efecto de cada factor. Quedan fuera por ahora la
profundidad desigual entre grupos y señales en taxones extremadamente raros.
Su incorporación se decidirá para la evaluación ampliada, no sobre la marcha.

Presupuesto por ejecución del piloto: **180 datasets, 720 ajustes del panel y 140 Wilcoxon = 860 ajustes**.
Los dos consensos reutilizan los resultados y no suman ajustes. El preflight
añadirá como máximo 20 ajustes; techo planificado total: **880**.

## 6. Universo, ejecución y fallos

No habrá filtrado por prevalencia dentro de cada réplica, salvo eliminar taxones
con conteo cero en **todas** las muestras mediante un único filtro común, ciego
al grupo. La preselección de G2 pertenece a la reserva, antes de la implantación.
No se aplicarán filtros distintos de prevalencia por motor. Las transformaciones
y la falta de estimaciones nativas se registrarán por separado.

Mantener verdad para los 100 taxones originales y para U. Para potencia
end-to-end, una señal eliminada o sin resultado es una señal no recuperada;
para diagnósticos condicionados, documentar explícitamente el denominador U.
Una fila ausente o con significancia NA no emite voto, pero tampoco se describe
como un test válido no significativo. Si `significant = TRUE` pero el efecto no
es finito o su dirección es neutral, marcar un fallo de contrato para este
protocolo y el consenso no disponible. El consenso actual puede agrupar la
dirección neutral: no borrar ni reemplazar esa llamada en el evaluador para
hacerla pasar como un voto direccional válido.

El runner aísla fallos y mide cada motor por separado. `prep()` propaga errores:
envolver solo la receta completa no permitiría atribuir todos los fallos. Se
conservan las salidas nativas de ajustes aislados y se construye el objeto
preparado de la receta común usando las validaciones del paquete. Esa
instrumentación es de desarrollo, no una ruta que altere resultados ni un
mecanismo para reutilizar ajustes de la versión anterior del protocolo.

En datasets completos, el resultado del `bake()/cool()` real deberá coincidir
exactamente con el conteo direccional explícito calculado solo como comprobación.
El control no sustituye la implementación que se evalúa.

Si un motor falla o alcanza un timeout, el consenso de cuatro motores queda
**no disponible**, no como una lista vacía con FDP cero. Se conserva el intento,
su dataset, duración, advertencias y causa. No cambiar semilla, bajar el umbral
de acuerdo ni sustituir el motor. Un reintento por un fallo externo documentado
usa los mismos datos/semillas y queda vinculado al intento original.

## 7. Métricas y resumen

Unidad: taxón × contraste; una sola dirección coherente por descubrimiento de
las reglas 3/4 y 4/4. Un caso con dirección ambigua o duplicados de clave es un
fallo de contrato, no dos descubrimientos.

Por dataset y procedimiento, para cada verdad aplicable:

| Métrica | Definición / interpretación |
|---|---|
| Descubrimientos R | Número de taxones seleccionados |
| FDP | `FP / max(R, 1)`; cero cuando un procedimiento válido no descubre nada |
| FDR empírica | Media de FDP entre réplicas válidas; no cociente de FP totales y R totales |
| Potencia | Taxones no nulos detectados / taxones no nulos según la verdad evaluada |
| Potencia direccional | Taxones no nulos detectados con signo correcto / taxones no nulos según la verdad evaluada |
| FDP direccional | `(FP + señales detectadas con signo incorrecto) / max(R, 1)` |
| Error de dirección condicionado | Signos incorrectos / señales verdaderas detectadas; NA si no hay ninguna |
| Disponibilidad | Señales retenidas, filas devueltas, q/efectos finitos y ausentes, completitud del panel |
| Coste | Tiempo por motor/dataset, CPU si está disponible y pico de memoria por proceso |

«Señales implantadas» es el denominador end-to-end para `beta_latent`, no para
todas las verdades: en P4/P9 son 10 señales latentes, pero 100 señales CLR sobre
el universo completo. Si un taxón queda fuera de U, conservar su verdad original
para la métrica end-to-end y declarar por separado cualquier puntuación
condicionada a U. No dividir descubrimientos CLR por las 10 perturbaciones latentes.

En nulos globales/condicionales, potencia y error de dirección son NA. La FDR
empírica es la proporción de datasets válidos con **al menos un descubrimiento**;
la tasa por taxón se mostrará solo como métrica secundaria.

Comparar mayoría y unanimidad con LinDA mediante diferencias **emparejadas por
dataset** de FDP y potencia direccional. Mostrar todos los componentes y sus
distribuciones; no construir una referencia «mejor motor de cada réplica».

Para nulos: intervalos binomiales exactos al 95 %. Para diferencias de FDP y
potencia: bootstrap emparejado de datasets, 2 000 remuestreos, intervalos
percentiles al 95 %, etiquetados exploratorios. No se usarán para seleccionar un
umbral ni para declaraciones de significación. No combinar celdas como si fueran
réplicas intercambiables de una misma población.

Informar N planificado, N válido y N emparejado. Ante fallos, además de la FDR
condicionada a casos válidos, mostrar el intervalo de identificación
`[sum(FDP_validos)/N, (sum(FDP_validos) + N_no_disponibles)/N]`.
Es un rango por resultados desconocidos, no un intervalo de confianza.
Para potencia end-to-end, mostrar también la versión conservadora con cero
recuperación en datasets no disponibles. No imputar FDP cero como resultado
principal de un fallo completo.

Con 20 réplicas, incluso cero datasets con hallazgos en un nulo da un límite
superior unilateral exacto al 95 % de aproximadamente **0.139**. El piloto no
puede certificar FDR ≤ 0.05. No se aplicarán los gates ni las tolerancias de la
baseline interna a estas decisiones.

## 8. Semillas, preflight y recursos

Semillas fijadas, conservadas de la versión 0.1 para esta enmienda:

- Datos: `20261014 + 10000 * cell_index + replicate`, con P1…P9 como índices 1…9.
- Motores: semilla de datos + `100000 * engine_index`, con orden ANCOM, LinDA,
  ALDEx2, MaAsLin3 y Wilcoxon. Registrar estado RNG y versión de R.
- Bootstrap del informe: `20271014 + cell_index`.
- Preflight separado: base `19800310`; no entra en estimaciones del piloto.

Una vez autorizada la ejecución, el preflight tendrá dos datasets pequeños de
señal fuerte equilibrada, uno por generador, cada uno ajustado dos veces con las
mismas semillas: 2 × 5 × 2 = 20 ajustes. Su objetivo es comprobar contratos,
orientación, argumentos efectivos, estabilidad reproducible y conservación de
IDs, no estimar rendimiento. Verificar que MaAsLin3 centra el coeficiente en la
etapa de contraste, que la salida de ANCOM coincide con su llamada robusta y
que LinDA utiliza el backend local verificado, no solo los argumentos previstos.

IDs exclusivos `PF_G1` y `PF_G2`, 100 taxones, 30 muestras por grupo, profundidad
mediana 15 000, diez señales alternadas ±`log(4)`, sin batch. Cada dataset usa
réplica 1 y su propio índice 1/2 en la fórmula de semilla del preflight; las dos
ejecuciones de cada motor reutilizan exactamente el mismo dataset y semilla.
La comprobación de centramiento observa el coeficiente y el RHS devueltos por el
runner nativo de MaAsLin3; no reconstruye una mediana con otra regla.

Congelar antes de este paso versiones, reserva G2, manifiesto de celdas y
configuraciones. Si falta una dependencia, detenerse y pedir autorización para
instalarla; no ejecutar un panel incompleto. Si un contrato falla, corregir solo
con autorización, versionar el cambio y repetir el preflight antes del piloto.

Plan de recursos inicial: ejecución local, máximo dos procesos concurrentes y
un hilo por proceso, sin paralelismo anidado. Timeout por ajuste: 10 minutos.
Techo del piloto: 24 horas de CPU estimadas/medidas, sin contratar cómputo ni
usar servicios externos. Registrar RAM disponible; no lanzar simultáneamente
procesos cuyo consumo previsto exceda el 75 % de ella. Tras los primeros cinco
datasets por celda, proyectar coste; si supera el techo, pausar y proponer una
revisión de recursos sin consultar ni cambiar resultados estadísticos.

La autorización actual no abarca lanzar el nuevo preflight ni el piloto. Ambos
requieren autorización de ejecución; el piloto exige además un receipt nuevo
aprobado que coincida con este contrato. La evaluación ampliada tampoco se lanza
automáticamente si el piloto resulta prometedor.

## 9. Criterios de decisión fijados de antemano

Los márgenes siguientes son **decisiones prácticas para un piloto**, no garantías
estadísticas ni umbrales calibrados de publicación.

### A. Viabilidad técnica: requisito previo

- Al menos 19/20 ajustes completos por motor/celda; también para Wilcoxon en
  sus celdas compatibles.
- Al menos 18/20 datasets con panel completo por celda.
- Al menos el 90 % de señales originales retenidas por el universo común en
  cada celda con señal; mostrar pérdidas de resultados nativos adicionalmente.
- Cero discrepancias inexplicadas de IDs, contraste, orientación o resultado de
  consenso entre `bake()/cool()` y su comprobación explícita.
- Hipótesis y universos interpretables, sin mezclar escalas de verdad.
- Coste proyectado dentro del techo y causas de fallo/ausencia documentadas.

Si falla A, decisión: **revisar la infraestructura/protocolo**, no declarar éxito
o fracaso científico ni eliminar el motor problemático para salvar el consenso.

### B. Señal de utilidad: ampliar para evaluar una posible ventaja

Usar exclusivamente mayoría 3/4 frente a la referencia LinDA, sin escoger entre
las dos reglas. Celdas principales: P2, P6 y P8; P3 es estrés secundario.

Debe haber resultados informativos: potencia direccional media de LinDA ≥ 0.20
en las celdas comparadas. Si no, el piloto es inconcluyente sobre ganancias;
solo se podrá proponer una enmienda prospectiva de observabilidad/intensidad.

Una vía de continuación exige el **mismo patrón** en P8 y en al menos una de
P2/P6, para no depender exclusivamente del generador paramétrico:

1. Caída media de FDP de al menos 0.05 puntos absolutos, conservando al menos
   el 80 % de la potencia direccional de LinDA; **o**
2. Aumento medio de potencia direccional de al menos 0.10 puntos absolutos,
   sin aumentar FDP más de 0.02 y con FDP media del consenso ≤ 0.10.

Además, no debe aparecer una pérdida de potencia direccional > 0.20 puntos
absolutos o un aumento de FDP > 0.05 en otra celda principal interpretable.
En P1/P5/P7, ≥ 4/20 datasets válidos con descubrimientos de la mayoría es una
alerta de inflación grave que impide esta vía hasta investigar la causa;
no significa que un resultado inferior a 4/20 certifique control del 5 %.

Mostrar incertidumbre, sensibilidad al dato individual y rendimiento frente a
todos los componentes. Si otro componente parece dominar al consenso, la
evaluación ampliada debe confrontarlo, no ignorarlo porque la referencia fija
era LinDA. Los márgenes de B solo justifican invertir en comprobar la señal.

### C. Límite científico útil: ampliar con otra pregunta, no con una promesa

Aunque no se cumpla B, puede merecer ampliación para estudiar límites si:

- Para la **mayoría 3/4**, P4 y P9 muestran un cambio ≥ 0.20 al pasar de verdad
  latente a CLR en la **misma métrica** (FDP o potencia direccional), y además
  hay desacuerdo observable entre los motores sobre al menos cinco taxones de
  media por dataset, explicable por su referencia de efecto y no errores de código.
  No basta el cambio de denominador: la diferencia algebraica de verdades ya se
  conoce; una ampliación investigaría sus consecuencias prácticas para la
  selección e interpretación de hallazgos;
  **o**
- La unanimidad pierde ≥ 0.30 puntos de potencia direccional frente a LinDA en
  P8 y al menos una de P2/P6, sin una reducción de FDP de al menos 0.05;
  **o**
- La alerta de inflación bajo nulo aparece en ambos brazos P1/P7 y sigue siendo
  explicable después de superar A.

Decisión: **proponer un estudio acotado de condiciones de fallo**, sujeto a
aprobación. No presentar esto como éxito de precisión del consenso.

### D. No ampliar el benchmark científico por ahora

Si A pasa pero no hay B ni C, y las celdas tienen señal informativa, no ampliar
por inercia. Si el consenso se limita a devolver menos taxones con pérdidas de
sensibilidad sin una mejora relevante, priorizar el artículo de software.
No cambiar pesos, panel o cutoff para fabricar una ventaja. Si hay un intervalo
demasiado amplio, clasificar la conclusión como inconcluyente, no equivalencia.

## 10. Entregables y eventual ampliación

Entregables futuros del piloto:

- Manifiestos de versiones, argumentos, semillas, taxones y reserva fuente.
- Conteos/metadata sintéticos y verdades por taxón/contraste, con hashes.
- Objetos preparados, recetas exportadas y resultados nativos/normalizados;
  sin alterar q, efectos o clasificaciones para conseguir coincidencia.
- Tabla de ejecución/ausencias y métricas por réplica/celda/procedimiento.
- Informe con plano FDP–potencia, coste, completitud, diagnóstico de estimandos
  y decisión A/B/C/D con sus denominadores.

Si se aprueba ampliar: congelar una nueva versión del protocolo antes de generar
datos finales; usar nuevas semillas y, cuando corresponda, reservas/cohortes
independientes. Las 20 réplicas del piloto no serán la evaluación final. Una
modificación motivada por el piloto es desarrollo, no confirmación independiente.

La ampliación deberá añadir profundidad desigual entre grupos y decidir si
requiere perfiles externos o un comparador
externo como ZicoSeq. No se exige incorporarlos al paquete. Longitudinalidad solo
se añadirá si forma parte de la pregunta o de las afirmaciones del manuscrito.

No se modifica ni se promueve `dev/scientific-validation/baselines/accepted-full.csv`.
La baseline de mantenimiento y la evaluación de rendimiento tienen fines distintos.

## 11. Fundamento y trazabilidad

- [Yang y Chen, 2022](https://link.springer.com/article/10.1186/s40168-022-01320-0):
  los consensos por votación no tienen garantizada una mejora de FDR bajo
  composicionalidad fuerte; justifica estudiar el intercambio error–potencia,
  no dar por demostrada una ventaja.
- [Wirbel et al., 2024](https://link.springer.com/article/10.1186/s13059-024-03390-9):
  el generador y su realismo influyen en los resultados; justifica un segundo
  brazo, sin atribuir a esta adaptación todas las propiedades de SIMBA.
- [MaAsLin3, tutorial oficial](https://www.bioconductor.org/packages/release/bioc/vignettes/maaslin3/inst/doc/maaslin3_tutorial.html):
  interpretación de mediana, coeficiente, abundancia y prevalencia.
- Código inspeccionado: `R/bake.R`, `R/plot_methods.R`, `R/model-adapters.R`,
  `R/engine-args.R`, funciones `step_*` y `dev/scientific-validation/R/`.

La nueva carpeta está bajo `dev/`, ya excluida de los builds del paquete. Este
protocolo no añade una API, una dependencia o un resultado científico al paquete.
