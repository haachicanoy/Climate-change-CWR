CWR Climate Change Analysis
H. Achicanoy & N. Castañeda
CIAT, 2015

RUTAS IMPORTANTES

# Proyecciones climáticas futuras
\\dapadfs\data_cluster_2\gcm\cmip5\downscaled

# Carpeta análisis 2014 de cambio climático
/curie_data/ncastaneda/threats

DESCRIPCIÓN DE LOS DATOS

Datos de ocurrencia: 26 cultivos con diferente número de parientes silvestres asociados. Cada archivo se encuentra escrito en formato .csv con las variables:

- Taxon. Nombre del taxón
- lon. Longitud medida en grados decimales
- lat. Latitud medida en grados decimales

Datos climáticos: 30 GCMs para el RCP 4.5 año 2050 que contienen las 19 variables climáticas de worldclim y la altitud además de las condiciones actuales con las mismas variables. Dichos datos se encuentran escritos en formato raster a una resolución de 5 arc min a nivel global.

DESCRIPCIÓN DE LAS ETAPAS A EJECUTAR

1. Por cada taxón cortar los rasters climáticos (por GCM y condiciones actuales) generando un buffer de 50 km alrededor de las coordenadas más extremas.
2. Generar los datos de background en los lugares definidos en la etapa anterior por taxón.
3. Por cada combinación taxón-GCM correr MaxEnt realizando una validación cruzada para 10 subparticiones de los datos (70% training, 30% testing).
4. Utilizar varios estadísticos de validación para evaluar los resultados (AUC, TSS, Kappa index, ...)
5. Guardar todos los resultados en formato NCDF
6. Realizar un análisis de varianza a partir de las estadísticas de evaluación por todas las combinaciones generadas
7. Ejecutar el ensamble de modelos para cada taxón
8. Evaluar el porcentaje de cambio o migración en la distribución de las especies bajo el escenario de cambio climático estudiado. Fijando un límite de migráción de ...

Problemas con: Phaseolus_persistentus (Crop: bean) [SOLUCIONADO]
...
Processing: bambara
Processing: bean
Error in (function (classes, fdef, mtable)  :
  unable to find an inherited method for function ‘raster’ for signature ‘"try-error"’
...

Continua el procesamiento a partir de Cajanus

Problemas con: Vigna_unguiculata_letouzeyi (Crop: cowpea) [SOLUCIONADO]
...
Processing: cowpea
Error in (function (classes, fdef, mtable)  :
  unable to find an inherited method for function ‘raster’ for signature ‘"try-error"’
In addition: Warning message:
In mclapply(1:length(fcRastersProc), cropRasters, mc.cores = length(fcRastersProc)) :
  all scheduled cores encountered errors in user code
...

Continua el procesamiento a partir de Daucus

Problemas con: Daucus_carota_drepanensis (Crop: daucus) [SOLUCIONADO]

Processing: daucus
Error in (function (classes, fdef, mtable)  :
  unable to find an inherited method for function ‘raster’ for signature ‘"try-error"’
In addition: Warning message:
In mclapply(1:length(fcRastersProc), cropRasters, mc.cores = length(fcRastersProc)) :
  all scheduled cores encountered errors in user code

Corrección de problemas identificados:
1. Especies con solo 1 coordenada geográfica
2. Especies donde la coordenadas extremas tienen una distancia cercana a 0 grados

# Identificación de un nuevo problema [Issue 1]

Processing: eggplantNHM
...
Processing: Solanum_nigriviolaceum [REVISAR!!!] [SOLUCIONADO] [Número reducido de coordenadas]
Error in .rasterObjectFromCDF(x, type = objecttype, band = band, ...) :
  cells are not equally spaced; you should extract values as points
In addition: Warning messages:
1: In min(rs) : no non-missing arguments to min; returning Inf
2: In max(rs) : no non-missing arguments to max; returning -Inf
3: In min(rs) : no non-missing arguments to min; returning Inf

Se deja este caso aparte y se continua procesando el resto de la información

Continua el procesamiento a partir de Eleusine

Problemas con: Lens_culinaris_tomentosus (Crop: lens) [SOLUCIONADO] [Número reducido de coordenadas]

Processing: Lens_culinaris_tomentosus
Error in .rasterObjectFromCDF(x, type = objecttype, band = band, ...) :
  cells are not equally spaced; you should extract values as points
In addition: Warning messages:
1: In min(rs) : no non-missing arguments to min; returning Inf
2: In max(rs) : no non-missing arguments to max; returning -Inf
3: In min(rs) : no non-missing arguments to min; returning Inf

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# Etapa de modelación
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

Comparación de tiempos de procesamiento

# Modelación
# Crop = 'avena', Taxon = 'Avena_abyssinica', GCM=1; Apróximadamente 40 segundos

# Proyección con predict function
      #    user  system elapsed
      # 425.976   2.470 398.787

# Proyección con make.projections function
      #   user  system elapsed
      # 83.118   3.869  87.053

DESCRIPCIÓN DE LA FUNCIÓN PARA MODELAR

Input: crop

Internamente se modela cada uno de los taxones asociados al crop. Para hacer esto necesitamos un loop o
paralelizar la función que realiza este paso. (PARA HACER)

Por taxón se realizarán 5 corridas mediante el procedimiento de validación cruzada. Las estadísticas de
cada corrida se deben almacenar en una tabla que debería contener:

- Regularized.training.gain: ganancia
- Training.AUC: AUC para los datos con los que se entrenó el modelo
- Test.AUC: AUC para los datos con los que se evaluó el modelo
- AUC.Standard.Deviation: desviación estándar de los AUC de entrenamiento
- Threshold: Esquina superior izquierda de la curva ROC

Falta calcular:

- Threshold (PARA HACER)

Se debe almacenar las estadísticas de todos los taxones en una tabla para obtener un resumen del cultivo.
(PARA HACER)
Adicional a esto, se podría crear un archivo de texto donde contemos con la información de: (PARA HACER)

- Número de presencias por taxón
- Verificar las features con los cuales se corrió el modelo
- ... (no se me ocurre nada más por el momento)

A partir de las 5 replicas realizadas por taxón, se calcula el promedio, el cual se umbraliza a partir
del threshold seleccionado (debe almacenarse), esto por cada uno de los 30 GCM's. Con estos 30 rasters
umbralizados calculamos el modelo ensamble a partir de la función de quantiles.

Archivos a almacenar en la carpeta: /curie_data/storage/climate_change/[crop]/maxent_modeling

- Tabla con estadísticas de resumen para el cultivo
- Archivo de texto con información de interés sobre las corridas del modelo

Archivos a almacenar en la carpeta: /curie_data/storage/climate_change/[crop]/maxent_modeling/[specie_name]

- [specie_name]_current.nc : contiene la distribución potencial umbralizada de la especie considerada en condiciones
actuales
- [specie_name]_future.nc  : contiene la distribución potencial umbralizada de la especie considerada en condiciones
futuras para cada uno de los 30 GCM's
- [specie_name]_ensemble.nc: contiene la distribución potencial umbralizada de la especie considerada en condiciones
futuras calculada mediante un procedimiento ensamble

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# Etapa de modelación incorporando nuevas ideas
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
@ Optimización de background
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

@@ Tener en consideración el artículo de Syfert et al (2013) e implementar la corrección por sampling bias
@@ en el código por cada especie

Tres escenarios en consideración
1. Polinomial: en base a ocurrencias extremas
2. Rectangular: en base a ocurrencias extremas a partir de los vertices
3. Modelación previa: utilizar un modelo de envelopment, BIOCLIM ...

Por cada escenario se tienen dos buffer para probar la bondad de ajuste del modelo por corrección de muestreo