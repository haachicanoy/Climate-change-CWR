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
8. Evaluar el porcentaje de cambio en la distribución de las especies bajo el escenario de cambio climático estudiado

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

# Identificación de un nuevo problema

Processing: eggplantNHM
...
Processing: Solanum_nigriviolaceum [REVISAR!!!] [PENDIENTE]
Error in .rasterObjectFromCDF(x, type = objecttype, band = band, ...) :
  cells are not equally spaced; you should extract values as points
In addition: Warning messages:
1: In min(rs) : no non-missing arguments to min; returning Inf
2: In max(rs) : no non-missing arguments to max; returning -Inf
3: In min(rs) : no non-missing arguments to min; returning Inf

Se deja este caso aparte y se continua procesando el resto de la información

Continua el procesamiento a partir de Eleusine