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

# Seleccionar puntos de background para cada especie (o cultivo ... preguntar)

sp_occ <- read.csv(...)
sp_occ <- unique(sp_occ)
extreme_coord <- c(min(sp_occ$lon),max(sp_occ$lon),min(sp_occ$lat),max(sp_occ$lat))
extreme_coord <- extreme_coord+0.4495392 # Buffer of 50 km
ext <- extent(extreme_coord); rm(extreme_coord)
mask.test <- crop(mask,ext); rm(ext)

sampleRandom(...)
