CWR Climate Change Analysis
H. Achicanoy & N. Castañeda
CIAT, 2015

RUTAS IMPORTANTES

# Proyecciones climáticas futuras
\\dapadfs\data_cluster_2\gcm\cmip5\downscaled

# Carpeta análisis 2014 de cambio climático
/curie_data/ncastaneda/threats

DESCRIPCIÓN DE LOS DATOS

Datos de ocurrencia: para 26 cultivos con diferente número de parientes silvestres asociados. Cada archivo se encuentra escrito en formato .csv con las variables:

- Taxon. Que corresponde al nombre del taxón
- lon. Longitud medida en grados decimales
- lat. Latitud medida en grados decimales

Datos climáticos: 30 GCMs para el RCP 4.5 año 2050 que contienen las 19 variables climáticas de worldclim y la altitud además de las condiciones actuales con las mismas variables. Dichos datos se encuentran escritos en formato raster a una resolución de 5 arc min a nivel global.

DESCRIPCIÓN DE LAS ETAPAS A EJECUTAR

1. Por cada taxón cortar los rasters climáticos (por GCM) generando un buffer de 50 km alrededor de las coordenadas más extremas.
2. 