# Climate change study
# H. Achicanoy & N. Castañeda
# CIAT, 2015

# Linux commands
# Verify space of storage
# df -h /curie_data/storage/

# Extraer información climatica de cada uno de los GCM para el escenario y periodo
# planteado

options(warn=-1)
gcm <- list.files("//dapadfs/data_cluster_2/gcm/cmip5/downscaled/rcp45/global_5min/", full.names=T)
gcm <- paste0(gcm,"/r1i1p1/2040_2069/")
gcm <- expand.grid(gcm,paste0("bio_",1:19))
gcm <- gcm[order(gcm[,1]),]
rownames(gcm) <- 1:nrow(gcm)
colnames(gcm) <- c("dir","biovariable")
gcm$GCM <- gcm$dir
gcm$GCM <- gsub(pattern="//dapadfs/data_cluster_2/gcm/cmip5/downscaled/rcp45/global_5min/",replacement="",gcm$GCM)
gcm$GCM <- gsub(pattern="/r1i1p1/2040_2069/",replacement="",gcm$GCM)
gcm$dir <- paste0(gcm$dir,gcm$biovariable)

gcmList <- unique(as.character(gcm$GCM))
env_dir <- "//dapadfs/workspace_cluster_6/CWR/CWR_PROJECT_CC_BD/future_clm_cimp5"

# Creating directories
lapply(gcmList,function(x){
  cat("\n\n\n Processing biovariables from GCM:",x,"\n\n")
  gcm_dir <- paste0(env_dir,"/",x)
  if(!file.exists(gcm_dir)){cat("Creating GCM directory\n"); dir.create(gcm_dir)} else{cat("GCM directory exists\n")}
}) # In this part is necessary to give permissions manually in local PC

# Creating ASCII files for maxent
# gcmList <- gcmList[-1]
library(raster)
lapply(gcmList,function(x){
  cat("\n\n\n Processing biovariables from GCM:",x,"\n\n")
  gcm_dir <- paste0(env_dir,"/",x)
  gcm_ras <- lapply(gcm$dir[which(gcm$GCM==x)], FUN=function(z){raster(paste0(z,"/w001001.adf"))})
  bios    <- as.character(gcm$biovariable[which(gcm$GCM==x)])
  cat("Printing results\n")
  mapply(x=gcm_ras,y=bios,function(x,y){writeRaster(x,paste0(gcm_dir,"/",y,".asc"),overwrite=T)})
})

rm(gcm); rm(bios); rm(env_dir); rm(gcmList)

#   for(i in 1:length(gcm_ras))
#   {
#     writeRaster(gcm_ras[[i]],paste0(gcm_dir,"/",bios[[i]],".asc"),overwrite=T)
#   }
# writeRaster(gcm_ras[[i]],paste0("\\\\dapadfs\\workspace_cluster_6\\CWR\\CWR_PROJECT_CC_BD\\future_clm_cimp5\\bcc_csm1_1\\",bios[[i]],".grid"),overwrite=F)
# writeRaster(gcm_ras[[i]],paste0("//dapadfs/workspace_cluster_6/CWR/CWR_PROJECT_CC_BD/future_clm_cimp5/",bios[[i]],".asc"),overwrite=F)
# writeRaster(gcm_ras[[i]],paste0("\\\\dapadfs\\workspace_cluster_6\\CWR\\CWR_PROJECT_CC_BD\\future_clm_cimp5\\bcc_csm1_1\\",bios[[i]],".tif"),format="GTiff",overwrite=T)
# writeRaster(gcm_ras[[i]],paste0("F:/",bios[[i]],".asc"),overwrite=T)

# Copiar información al directorio de destino
# Linux commands
cp -r /mnt/workspace_cluster_6/future_clm_cimp5* /curie_data/storage
# Eliminar información del directorio inicial
rm -r /mnt/workspace_cluster_6/future_clm_cimp5

# Copiar occurrence files al directorio de destino

main_dir <- "/curie_data/storage"
cc_dir <- paste0(main_dir,"/climate_change")
if(!file.exists(cc_dir)){dir.create(cc_dir)} else{"Climate change directory exists"}

cropList <- list.files("/curie_data/ncastaneda/threats")
cropList <- cropList[2:29]
cropList <- cropList[-12]

# Crear directorios de destino y organizar datos de ocurrencias

lapply(cropList, function(crop){
  cat("\n %%%%%%%%%%%% PROCESSING:",crop,"CROP %%%%%%%%%%%%\n\n")
  
  cat(" %%%%%%%% Creating directories\n")
  crop_dir <- paste0(cc_dir,"/",crop)
  if(!file.exists(crop_dir)){dir.create(crop_dir)} else{cat("Occurrence directory exists\n")}
  
  cat(" %%%%%%%% Arranging occurrence data\n")
  occTaxas <- list.files(path=paste0("/curie_data/ncastaneda/threats/",crop,"/occurrence_files"),pattern=".csv$")
  occTaxas <- gsub(pattern=".csv",replacement="",occTaxas)
  occFiles <- list.files(path=paste0("/curie_data/ncastaneda/threats/",crop,"/occurrence_files"),pattern=".csv$",full.names=T)
  mapply(x=occFiles, y=occTaxas, FUN=function(x,y){
    cat(" ... Processing:",y,"\n")
    data <- read.csv(x)
    data <- data[,c("Taxon","lon","lat")]
    data <- unique(data)
    rownames(data) <- 1:nrow(data)
    occ_dir <- paste0(crop_dir,"/occurrence_files")
    if(!file.exists(occ_dir)){dir.create(occ_dir)} else{cat("Occurrence directory exists\n")}
    write.csv(data, paste0(occ_dir,"/",y,".csv"),row.names=F)
  })
  
})

# Include altitude layer with the same resolution than others rasters
alt <- "/curie_data/ncastaneda/geodata/alt_2-5m/aligned-bio/alt.asc"
alt <- raster(alt)
bio_example <- raster("/curie_data/storage/future_clm_cimp5/bcc_csm1_1/bio_1.asc") # Whatever layer is good
alt2 <- raster::resample(alt,bio_example,method="bilinear")
rm(bio_example,alt)
gcm_dirs <- list.dirs("/curie_data/storage/future_clm_cimp5",full.names=T,recursive=F)
library(parallel)
writeRAS <- function(i){writeRaster(alt2,paste0(gcm_dirs[[i]],"/bio_0.asc"),overwrite=T)}
mclapply(1:length(gcm_dirs), writeRAS, mc.cores=10)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

## Realizar cortes por coordenadas de ocurrencias extremas

crops <- list.dirs("/curie_data/storage/climate_change",full.names=F,recursive=F) # Crops to analyse

## Function to crop the future and current climate

# Read future and current climate information for each model
library(parallel)
fcList <- list.dirs(paste0("/curie_data/storage/future_clm_cimp5"),full.names=T,recursive=F)
readRasters <- function(i)
{
  bios <- list.files(fcList[[i]],pattern=".asc$",full.names=T)
  rasters <- lapply(1:length(bios), function(j){library(raster); r <- raster(bios[[j]]); return(r)})
  rasters <- stack(rasters)
  return(rasters)
}
fcRasters <- mclapply(1:length(fcList), readRasters, mc.cores=20)
models <- list.dirs(paste0("/curie_data/storage/future_clm_cimp5"),full.names=F,recursive=F)
names(fcRasters) <- models
g <- gc()
rm(g)

# Function to crop raster files
createBackFile <- function(crop) # Run in line
{
  library(raster)
  library(ncdf)
  cat("Processing:",crop,"\n")
  spList <- list.files(paste0("/curie_data/storage/climate_change/",crop,"/occurrence_files"),pattern=".csv$",full.names=T)
  spName <- list.files(paste0("/curie_data/storage/climate_change/",crop,"/occurrence_files"),pattern=".csv$",full.names=F)
  spName <- gsub(pattern=".csv",replacement="",spName)
  lapply(1:length(spList), function(i)
  {
    cat("Processing:",spName[[i]],"\n")
    sp_occ <- read.csv(spList[[i]]) # Reading occurrence data
    sp_occ <- unique(sp_occ)
    
    if(nrow(sp_occ)>15) # Only more than 15 coordinates
    {
      extreme_coord <- c(min(sp_occ$lon),max(sp_occ$lon),min(sp_occ$lat),max(sp_occ$lat)) # Identify extreme coordinates
      extreme_coord <- extreme_coord+0.4495392 # Buffer of 50 km
      extreme_coord <- data.frame(lon=c(extreme_coord[1],extreme_coord[2]),lat=c(extreme_coord[4],extreme_coord[3]))
      
      eucDist <- sqrt(sum((extreme_coord[1,]-extreme_coord[2,])^2)); eucDist <- round(eucDist) # Calculate distance between extrem coordinates
      if(eucDist>0)
      {
        
        # Create directories
        fcDir <- paste0("/curie_data/storage/climate_change/",crop,"/future_climate")
        if(!file.exists(fcDir)){dir.create(fcDir)}
        fcDirSp <- paste0(fcDir,"/",spName[[i]])
        if(!file.exists(fcDirSp)){dir.create(fcDirSp)}
        
        lapply(1:length(fcRasters), function(k) # Run in line
        {
          fcRastersProc <- unstack(fcRasters[[k]])
          cells         <- cellsFromExtent(fcRastersProc[[1]],extent=extent(fcRastersProc[[1]])) # Cropping rasters
          infoCells     <- data.frame(cell=cells,xyFromCell(fcRastersProc[[1]],cells,spatial=F)); cells <- NULL
          infoCells     <- base::subset(infoCells,subset= x>=extreme_coord[1,1] & x<=extreme_coord[2,1] & y<=extreme_coord[1,2] & y>=extreme_coord[2,2], select=c(cell,x,y))
          library(parallel) # Run in parallel
          cropRasters <- function(j){z <- rasterFromCells(fcRastersProc[[j]], cells=infoCells$cell); z[] <- fcRastersProc[[j]][][infoCells$cell]; return(z)}
          fcRastersProc <- stack(mclapply(1:length(fcRastersProc),cropRasters,mc.cores=length(fcRastersProc)))
          fcRastersCrop <- fcRastersProc; rm(fcRastersProc)
          writeRaster(fcRastersCrop, filename=paste0(fcDirSp,"/",models[[k]],".nc"), format="CDF", overwrite=T) # Save results
          g <- gc(); rm(g)
        }
        )
        
      } else {cat("Distance between extreme coordinates is 0!\n")}
    } else {cat("Number of coordinates is limited for the analysis!\n")}
  }
  )
}

# Apply to all crops
lapply(1:length(crops),function(i){createBackFile(crops[[i]]); return("Done")})

# Function to create background files for each specie
# Step 1. Leer datos de ocurrencia para cada especie
# Step 2. Leer solo 1 raster de clima generado por especie
# Step 3. A partir de este raster generar 10000 puntos aleatorios como background
# Step 4. Guardar dichas coordenadas

crops <- list.dirs("/curie_data/storage/climate_change",full.names=F,recursive=F) # Crops to analyse

library(parallel)
backGenFunc <- function(i)
{
  library(dismo)
  cat("Processing:",crops[[i]],"\n")
  spList <- list.files(paste0("/curie_data/storage/climate_change/",crops[[i]],"/future_climate"),full.names=T,recursive=F)
  spList <- paste0(spList,".csv")
  spList <- gsub(pattern="future_climate",replacement="occurrence_files",spList)
  spName <- list.files(paste0("/curie_data/storage/climate_change/",crops[[i]],"/future_climate"),full.names=F)
  lapply(1:length(spList),function(j)
  {
    cat("Processing:",spName[[j]],"\n")
    sp_occ <- read.csv(spList[[j]])
    sp_occ <- unique(sp_occ)
    
    if(nrow(sp_occ)>15){
      # Read current climate
      r_file <- raster(paste0("/curie_data/storage/climate_change/",crops[[i]],"/future_climate/",spName[[j]],"/current.nc"))
      backFile <- as.data.frame(dismo::randomPoints(mask=r_file,n=10000,p=sp_occ[,c("lon","lat")])) # Generate background points
      backFile$Taxon <- spName[[j]]
      names(backFile)[1:2] <- c("lon","lat")
      backFile <- backFile[,c("Taxon","lon","lat")]
      backDir <- paste0("/curie_data/storage/climate_change/",crops[[i]],"/backgroundFiles")
      if(!file.exists(backDir)){dir.create(backDir)}
      write.csv(backFile,paste0(backDir,"/",spName[[j]],".csv"),row.names=F)
    } else {
      cat("Number of coordinates is limited for the analysis!\n")
    }
    
  })
  return("Done")
}
mclapply(1:length(crops),backGenFunc,mc.cores=20)

# Estrategia de modelación
# Step 1. Leer datos de ocurrencia por especie
# Step 2. Leer datos de background por especie
# Step 3. Leer información climática, presente y futuro (31 escenarios)
# Step 4. Correr algoritmo MaxEnt para los 31 diferentes escenarios por especie
# replicando 5 veces cada modelo mediante validación cruzada
# Step 5. Almacenar los resultados por cultivo, especie, modelo en formato .nc

crops <- list.dirs("/curie_data/storage/climate_change",full.names=F,recursive=F) # Crops to analyse

options(warn=-1)
if(!require(ff)){install.packages("ff"); library(ff)} else {library(ff)}
if(!require(ncdf)){install.packages("ncdf"); library(ncdf)} else {library(ncdf)}
if(!require(dismo)){install.packages("dismo"); library(dismo)} else {library(dismo)}
if(!require(parallel)){install.packages("parallel"); library(parallel)} else {library(parallel)}
if(!require(data.table)){install.packages("data.table"); library(data.table)} else {library(data.table)}
if(!require(PresenceAbsence)){install.packages("PresenceAbsence"); library(PresenceAbsence)} else {library(PresenceAbsence)}

stg.dir <- '/curie_data/storage'
src.dir <- paste(stg.dir,'/_scripts',sep='')

modelingStep <- function(crop)
{
  cat("Processing:",crop,"\n")
  # Species name
  spName <- list.files(paste0("/curie_data/storage/climate_change/",crop,"/backgroundFiles"),pattern=".csv",full.names=F)
  if(length(spName)!=0){
    spName <- gsub(pattern=".csv",replacement="",spName)
    
    # Occurrence data x specie
    taxList <- list.files(paste0("/curie_data/storage/climate_change/",crop,"/backgroundFiles"),pattern=".csv",full.names=T)
    taxList <- gsub(pattern="backgroundFiles",replacement="occurrence_files",taxList)
    taxList <- lapply(1:length(taxList),function(i){z <- read.csv(taxList[[i]]); return(z)})
    
    # Background data x specie
    bckList <- list.files(paste0("/curie_data/storage/climate_change/",crop,"/backgroundFiles"),pattern=".csv",full.names=T)
    bckList <- lapply(1:length(bckList),function(i){z <- read.csv(bckList[[i]]); return(z)})
    
    # Climate data x specie x GCM
    climDirSp <- paste0("/curie_data/storage/climate_change/",crop,"/future_climate/",spName)
    climData <- lapply(1:length(climDirSp),function(i){list.files(climDirSp[[i]],pattern=".nc$",full.names=T)})
    loadRasters <- function(i){lapply(1:length(climData[[i]]),function(j){stack(lapply(1:20,function(k){raster(paste(climData[[i]][[j]]),band=k)}))})}
    climData <- mclapply(1:length(climData),loadRasters,mc.cores=5)
    
    ## 1. Cross-validation process
    
    # Run MaxEnt with complete dataset of occurrences cross-validating by 5 folds
    # Using features: linear, quadratic and product
    # First index corresponds to taxon information, Second index corresponds to climatic information model
    fit <- dismo::maxent(x=climData[[7]][[4]],p=taxList[[7]][,c("lon","lat")],a=bckList[[7]][,c("lon","lat")],removeDuplicates=T,args=c("nowarnings","replicates=5","linear=true","quadratic=true","product=true","threshold=false","hinge=false"))
    source(paste(src.dir,'/do_projections.R',sep=''))
    cross.val.prj <- lapply(1:5,make.projections)
    cross.val.prj <- stack(cross.val.prj)
    
    ## 2. Testing
    metrics <- as.data.frame(fit@results)
    metrics$Metric <- rownames(metrics); rownames(metrics) <- 1:nrow(metrics)
    colnames(metrics) <- c(paste('Fold_',1:5,sep=''),'Fold_mean','Metric')
    metrics <- metrics[,c('Metric',paste('Fold_',1:5,sep=''),'Fold_mean')]
    # FALTA: 1. Calcular e incluir -> UpperLeftROC threshold
    metrics <- metrics[complete.cases(match(metrics$Metric,c('Regularized.training.gain','Training.AUC','Test.AUC','AUC.Standard.Deviation'))),]
    # FALTA: 2. Guardar resultados, mapas y metricos
    
    ## 3. Ensemble process
    
  } else {
    cat(crop,"doesn't have sufficient information for analyze\n")
  }
  return(cat("Done!\n"))
}

compareCurFut <- function(crop)
{
  # Compare current and future conditions
  # with percent change in species distribution
}
