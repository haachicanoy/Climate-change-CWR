# Arrange current climate
# H. Achicanoy
# CIAT, 2015

# In local PC

path <- "C:/Users/haachicanoy/Documents/current_clim"
library(raster)

var <- list.files(path,pattern=".bil$",full.names=T)
nam <- list.files(path,pattern=".bil$",full.names=F)
nam <- gsub(pattern=".bil",replacement="",nam)
var <- lapply(1:length(var),function(i){raster(var[[i]])})
if(!file.exists(paste0(path,"/ascii"))){dir.create(paste0(path,"/ascii"))}
lapply(1:length(var),function(i){writeRaster(var[[i]], filename=paste0(path,"/ascii/",nam[[i]],".asc"),format="ascii",overwrite=T);return("Done!")})

# Curie

path <- "/curie_data/storage/future_clm_cimp5/current"
library(raster)

var <- list.files(path,pattern=".asc$",full.names=T)
var <- lapply(1:length(var),function(i){raster(var[[i]])})
var <- stack(var)