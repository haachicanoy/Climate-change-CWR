# The process

crop <- "rice"

require(raster)
require(dismo)
require(rgdal)

# PREPARING CLIMATE LAYERS
# ==================================
cur.clim.dir <- paste("/curie_data/ncastaneda/gap-analysis/gap_",crop,"/biomod_modeling/current-clim",sep="")
# cur.clim.dir <- paste("/curie_data/cip/gap_potatoSpooner/biomod_modeling/current-clim/") # POTATO CIP

fut.clim.dir <- "/curie_data/ncastaneda/geodata/future_clm/cmip5/rcp4_5_2040_2069_ensemble"

mask <- paste("/curie_data/ncastaneda/gap-analysis/gap_",crop,"/masks/mask.asc",sep="")
# mask <- paste("/curie_data/cip/gap_ipomoeaCIP/masks/mask.asc") # POTATO CIP
mask <- raster(mask)
mask.e <- extent(mask)

o.dir.f <- paste("/curie_data/ncastaneda/threats/",crop,"/clim-chg/future-clim", sep="")
if(!file.exists(o.dir.f)){dir.create(o.dir.f)}

# Extracting future bioclim layers
asc.list <- c(2,3,8,9,13,14,15,18,19)
for(i in asc.list){
  cat("Preparing the layer bio_",i," for ",crop, "\n")
  rf <- raster(paste(fut.clim.dir,"/bio_",i,".asc",sep=""))
  rf <- crop(rf,mask.e)
  
  writeRaster(rf,paste(o.dir.f,"/bio_",i,".asc",sep=""), overwrite=T)
}

# Extracting elevation raster for future layers
alt <- "/curie_data/ncastaneda/geodata/alt_2-5m/alt.asc"
alt <- raster(alt)
alt.e <- crop(alt,mask.e)
writeRaster(alt.e,paste(o.dir.f,"/alt.asc",sep=""), overwrite=T)
rm(alt.e,alt)

# Extracting elevation raster for current layers
alt <- "/curie_data/ncastaneda/geodata/alt_2-5m/aligned-bio/alt.asc"
alt <- raster(alt)
alt.e <- crop(alt,mask.e)
writeRaster(alt.e,paste(cur.clim.dir,"/alt.asc",sep=""), overwrite=T)


# PREPARING SWD FILES
# ==================================

# Preparing occurrences swd (unique data is better)
main.dir <- paste("/curie_data/ncastaneda/threats/",crop,sep="")
occ_dir <- paste(main.dir,"/occurrence_files",sep="")
occList <- list.files(occ_dir, pattern=".csv", full.names=T)
for(occFile in occList){
  # Load data
  occ.data <- read.csv(occFile)
  occ.data <- occ.data[c("Taxon","lon","lat","bio_2","bio_3","bio_8","bio_9","bio_13","bio_14","bio_15","bio_18","bio_19")]
  xy <- occ.data[c("lon","lat")]
  occ.data$alt <- extract(alt,xy)
  occ.data <- occ.data[complete.cases(occ.data),]
  write.csv(occ.data,occFile, row.names=FALSE)
}

# Preparing background swd
eco <- paste("/curie_data/ncastaneda/gap-analysis/gap_",crop,"/masks/mask.asc",sep="")
# eco <- paste("/curie_data/cip/gap_potatoSpooner/masks/mask.asc") # POTATO CIP

eco <- raster(eco)

sppList <- list.files(occ_dir, pattern=".csv", full.names=T)
for(spp in sppList){
  #confirming if the species has enough data for modelling
  occ <- read.csv(spp)
  nocc <- nrow(occ)
  if(nocc > 10){
    spname <- tail(unlist(strsplit(spp,split="/")), n=1)
    spname <- substr(spname,1,nchar(spname)-4)
    size  = 10000
    cat("Creating", size, "background points for",spname, "\n")
    bck.points <- as.data.frame(sampleRandom(eco, size=size, xy=T, na.rm=T))
    
    asc.list <- c("bio_2", "bio_3", "bio_8", "bio_9", "bio_13", "bio_14", "bio_15", "bio_18", "bio_19", "alt")
    
    xy <- data.frame(X=bck.points[,"x"],Y=bck.points[,"y"])
    swd <- xy
    env_dir <- cur.clim.dir
    
    for(i in asc.list){
      cat("Reading environmental layer",i,"\n")
      rs <- raster(paste(env_dir,"/",i,".asc",sep=""))
      swd$NEW <- extract(rs,xy)
      k <- which(asc.list == i)
      names(swd)[k+2] <- i
    }
    swd <- swd[complete.cases(swd),]
#     outdir <- paste(main.dir,"/clim-chg/",spname,sep="")
#     if(!file.exists(outdir)){dir.create(outdir);dir.create(paste(outdir,"/background",sep=""))}
    write.csv(swd, paste(main.dir,"/clim-chg/",spname,"/background/bck.points.swd.csv", sep=""), row.names=T)
  }else{
    cat("Species with 10 or less coordinates \n")
  }
}

# MODELLING IN MAXENT
# ==================================

inputDir <- paste("/curie_data/ncastaneda/gap-analysis/gap_",crop,sep="")
cur.clim.dir <- paste(inputDir,"/biomod_modeling/current-clim",sep="")
maxentApp <- paste(inputDir, "/maxent_modeling/lib/maxent.jar", sep="")
NADir <- paste(inputDir, "/biomod_modeling/native-areas/asciigrids", sep="")
src.dir <- paste("/curie_data/ncastaneda/gap-analysis/gap_",crop,"/_scripts",sep="")
source(paste(src.dir,"/000.zipRead.R",sep=""))
OSys <- "linux"

main.dir <- paste("/curie_data/ncastaneda/threats/",crop,sep="")
occ_dir <- paste(main.dir,"/occurrence_files",sep="")

o.dir.f <- paste(main.dir,"/clim-chg/future-clim", sep="")

j.size <- "-mx2048m"
# j.size <- "-mx8192m"

climateChangeDist <- function(crop,spID,main.dir, occFile, j.size, maxentApp, OSys, cur.clim.dir, o.dir.f, NADir){
  
#   src.dir <- paste("/curie_data/ncastaneda/gap-analysis/gap_",crop,"/_scripts",sep="")
  src.dir <- "/curie_data/cip/gap_potatoSpooner/_scripts"
  source(paste(src.dir,"/000.zipRead.R",sep=""))
  
  cat("Processing",spID," \n")
  verFile <- paste(main.dir,"/clim-chg/",spID,"/maxent-projections", "/ps-", spID, ".run", sep="")
  if(!file.exists(verFile)){
    cat("Loading data... \n")
    inData <- read.csv(occFile)
    nOcc <- nrow(inData)
    if(nOcc > 10){
      outFileName <- occFile
      
      backoutdir <- paste(main.dir,"/clim-chg/",spID,"/background/", sep="")
      backFileSwd <- paste(backoutdir, "bck.points.swd.csv", sep="")
      
      spDir <- paste(main.dir,"/clim-chg/",spID,sep="")
      outFolder <- paste(spDir, "/maxent-projections", sep="")
      if (!file.exists(outFolder)) {
        dir.create(outFolder)
        dir.create(paste(spDir, "/maxent-metrics", sep=""))
      }
      
      cat("Fitting the model... \n")
      if (!file.exists(paste(outFolder,"/",spID,".html",sep=""))) {
        system(paste("java", j.size, "-jar", maxentApp, "-s", outFileName, "-e", backFileSwd, "-o", outFolder, "-P", "randomtestpoints=25", "nowarnings", "-a", "-z"), wait=TRUE)        
      }
      
      if (file.exists(paste(outFolder, "/", spID,".html", sep=""))) {
        cat("Model fitted successfully!", "\n")
        procSwitch <- T
      } else {
        cat("Error in computing... erasing the folder \n")
        if (OSys=="linux") {
          system(paste("rm", "-rv", outFolder))
        } else {
          unlink(outFolder)
        }
        procSwitch <- F
      }
      
      if(procSwitch){
        
        projectionList <- c("current","rcp4_5")
        cat("Projecting the model...", "\n")
        
        for (prj in projectionList) {
          
          if(prj == "current"){
            cat("Performing ", prj, "\n")
            projLayers <- cur.clim.dir
            suffix <- gsub("/", "_", prj)
            outGrid <- paste(outFolder, "/", spID, "_", suffix, sep="")
            lambdaFile <- paste(outFolder,"/",spID, ".lambdas",sep="")
            
            if (!file.exists(paste(outGrid,".asc",sep=""))) {
              bioexcluded <- paste("-N bio_1", "-N bio_4", "-N bio_5", "-N bio_6", "-N bio_7", "-N bio_10", "-N bio_11",
                                   "-N bio_12", "-N bio_16")
              
              system(paste("java", j.size, "-cp", maxentApp, "density.Project", lambdaFile, projLayers, outGrid, "nowarnings", "-r", "-a", "-z", bioexcluded), wait=TRUE)
            }
            if (file.exists(paste(outGrid, ".asc", sep=""))) {
              cat("Projection is OK!", "\n")
            } else {
              cat("Error in projecting", "\n")
            }
          }else{
            cat("Performing ", prj, "\n")
            projLayers <- o.dir.f
            suffix <- gsub("/", "_", prj)
            outGrid <- paste(outFolder, "/", spID, "_", suffix, sep="")
            lambdaFile <- paste(outFolder,"/",spID, ".lambdas",sep="")
            
            if (!file.exists(paste(outGrid,".asc",sep=""))) {
              system(paste("java", j.size, "-cp", maxentApp, "density.Project", lambdaFile, projLayers, outGrid, "nowarnings", "-r", "-a", "-z"), wait=TRUE)
            }
            if (file.exists(paste(outGrid, ".asc", sep=""))) {
              cat("Projection is OK!", "\n")
            } else {
              cat("Error in projecting", "\n")
            }
          }
        }
        
        cat("Thresholding the distribution models \n")
        threshData <- read.csv(paste(outFolder,"/maxentResults.csv",sep=""))
        thslds <- c("Maximum.training.sensitivity.plus.specificity.logistic.threshold")
        
        thrNames <- names(threshData)
        thePos <- which(thrNames == thslds)
        theVal <- threshData[1,thePos]
        
        currDist <- raster(paste(outFolder,"/",spID,"_current.asc",sep=""))
        currDistPR <- currDist
        currDistPR[which(currDistPR[] < theVal)] <- NA
        
        currDistPA <- currDist
        currDistPA[which(currDistPA[] < theVal)] <- 0
        currDistPA[which(currDistPA[] != 0)] <- 1
        
        futDist <- raster(paste(outFolder,"/",spID,"_rcp4_5.asc",sep=""))
        futDistPR <- futDist
        futDistPR[which(futDistPR[] < theVal)] <- NA
        
        futDistPA <- futDist
        futDistPA[which(futDistPA[] < theVal)] <- 0
        futDistPA[which(futDistPA[] != 0)] <- 2
        
        # Restricting distributions to known native areas
        NAGridName <- paste(NADir, "/", spID, "/narea.asc.gz", sep="")
        if (!file.exists(NAGridName)) {
          cat("The native area does not exist, generating one \n")
          NAGrid <- chullBuffer(inputDir, occFile, paste(NADir, "/", spID, sep=""), 500000)
          zipWrite(NAGrid, paste(NADir, "/", spID, sep=""), "narea.asc.gz")
        } else {
          cat("The native area exists, using it \n")
          NAGrid <- zipRead(paste(NADir, "/", spID, sep=""), "narea.asc.gz")
        }
        NAGrid[which(NAGrid[] != 1)] <- NA
        
        currDistPA <- currDistPA * NAGrid
        futDistPA <- futDistPA * NAGrid
        
        writeRaster(currDistPA,paste(outFolder,"/",spID,"_current-PA.tif",sep=""),overwrite=T)
        writeRaster(futDistPA,paste(outFolder,"/",spID,"_rcp4_5-PA.tif",sep=""),overwrite=T)
        
        # Analyzing change
        chng.rs <-  sum(currDistPA, futDistPA)
        chng.rs[which(chng.rs[] == 0)] <- NA
        
        writeRaster(chng.rs,paste(outFolder,"/",spID,"_clim-chng-PA.tif",sep=""))
        
        ftoTIF <- list.files(outFolder, pattern=".asc")
        cat("Converting to TIF ... \n")
        for (fz in ftoTIF) {
          fName <- paste(outFolder, "/", fz, sep="")
          if (OSys == "linux") {
            fRaster <- raster(fName)
            fz <- strsplit(fz,split = ".asc")[[1]]
            writeRaster(fRaster,paste(outFolder,"/",fz,".tif",sep=""), overwrite=T)
            file.remove(fName)
          } else {
            setwd("C:/Program Files/7-Zip")
            system(paste('7z.exe a -tgzip "',gsub(paste(fName,".gz",sep=""),pattern="/",replacement="\\\\"),'" "',gsub(fName,pattern="/",replacement="\\\\"),'"',sep=''),wait=T)
            file.remove(fName)
          }
        }
        
        #Run verification file
        verFile <- paste(outFolder, "/ps-", spID, ".run", sep="")
        opnFile <- file(verFile, open="w")
        cat("Modelled on", date(), file=opnFile)
        close.connection(opnFile)             
        return("Done")
        
      }else{
        cat("Species with invalid maxent model \n")
      }
    }else{
      cat("This species has less than 10 coords \n")
    }
  }else{
    cat("This species has been modelled already \n")
  }
}

main.dir <- paste("/curie_data/ncastaneda/threats/",crop,sep="")
occ_dir <- paste(main.dir,"/occurrence_files",sep="")
occList <- list.files(occ_dir, pattern=".csv", full.names=T)
for(occFile in occList){
  spID <- tail(unlist(strsplit(occFile,split="/")), n=1)
  spID <- substr(spID,1,nchar(spID)-4)
  
  main.dir <- paste("/curie_data/ncastaneda/threats/",crop,sep="")
  j.size <- "-mx1024m"
  inputDir <- paste("/curie_data/ncastaneda/gap-analysis/gap_",crop,sep="")
#   inputDir <- "/curie_data/cip/gap_potatoSpooner"
  maxentApp <- paste(inputDir, "/maxent_modeling/lib/maxent.jar", sep="")
  OSys <- "linux"
  cur.clim.dir <- paste(inputDir,"/biomod_modeling/current-clim",sep="")
  o.dir.f <- paste(main.dir,"/clim-chg/future-clim", sep="")
  NADir <- paste(inputDir, "/biomod_modeling/native-areas/asciigrids", sep="")
  
  x <- climateChangeDist(crop=crop,spID=spID, main.dir=main.dir, occFile=occFile, j.size=j.size, maxentApp=maxentApp, OSys=OSys,
                         cur.clim.dir=cur.clim.dir, o.dir.f=o.dir.f, NADir=NADir)
}

