#==========================
# SUMMARISING RESULTS
#==========================
cropList= c("avena", "bambara", "bean", "cajanusPaper", "cicer", "cowpea", "daucus", "eggplantNHM", "eleusine", "faba_bean", "helianthusSoils", "hordeum", 
                 "ipomoea", "lathyrus", "lens", "lima_bean", "malus","medicago","musa","pennisetum", "pisum", "potatoCIP", "rice", "secale", "sorghum", 
                 "triticum", "vetch")

for(crop in cropList){
  main.dir <- paste("/curie_data/ncastaneda/threats/",crop,sep="")
  occ_dir <- paste(main.dir,"/occurrence_files",sep="")
  occList <- list.files(occ_dir, pattern=".csv", full.names=T)
  
  i <- 1
  for(occFile in occList){
    spID <- tail(unlist(strsplit(occFile,split="/")), n=1)
    spID <- substr(spID,1,nchar(spID)-4)
    
    mxnt.dir <- paste(main.dir,"/clim-chg/",spID,"/maxent-projections", sep="")
    verFile <- paste(mxnt.dir,"/ps-",spID,".run",sep="")
    
    if(file.exists(verFile)){
      
      cat("Reading data for",spID,"\n")
      mxnt.file <- paste(mxnt.dir,"/maxentResults.csv",sep="")
      mxnt.file <- read.csv(mxnt.file)
      mxnt.file$CROP_CODE <- crop
      
      if(i==1){
        allMxnt <- mxnt.file
      }else{
        allMxnt <- rbind(allMxnt,mxnt.file)
      }
      i <- i+ 1
      
      outDir <- paste(main.dir,"/clim-chg/_summaries/",sep="");if(!file.exists(outDir)){dir.create(outDir)}
      outFile <- paste(outDir,"allMaxentResults.csv",sep="")
      write.csv(allMxnt,outFile, row.names=F)
      
    }else{
      cat("Species not modelled \n")
    }
  }
}

# All crops summarized - producing file to guide species richness

cropList= c("avena", "bambara", "bean", "cajanusPaper", "cowpea", "daucus", "eggplantNHM", "eleusine", "helianthusSoils", "hordeum", 
            "ipomoea", "lathyrus", "lens", "lima_bean", "malus","medicago","musa","pennisetum", "pisum", "potatoCIP", "rice", "secale", "sorghum", 
            "triticum", "vetch") # cicer and faba_bean were removed from the list as no species were modelled

for(crop in cropList){
  main.dir <- paste("/curie_data/ncastaneda/threats/",crop,sep="")
  outFile <- paste(main.dir,"/clim-chg/_summaries/allMaxentResults.csv",sep="")
  outFile <- read.csv(outFile)
  
  forRich <- outFile[,c("Species","Test.AUC","CROP_CODE")]
  forRich$IS_VALID <- NA
  
  species <- forRich$Species
  for(sp in species){
    
    iv <- forRich$Test.AUC[which(forRich$Species==paste(sp))]
    cat(sp,iv,"\n")
    if(iv < 0.7){
      hs <- 0
      forRich$IS_VALID[which(forRich$Species==paste(sp))] <- iv
    }else{
      iv <- 1
      forRich$IS_VALID[which(forRich$Species==paste(sp))] <- iv
    }
  }
  outDir <- paste(main.dir,"/clim-chg/_summaries/",sep="")
  write.csv(forRich,paste(outDir,"speciesRichness.csv",sep=""),row.names=F) 
}

# Preparing richness maps

speciesRichness_alt <- function(crop,bdir,prj) {
  idir <- paste(bdir, "/clim-chg/_summaries", sep="")
  #   mask <- raster(paste(bdir, "/masks/mask.asc", sep="")) # New line
  
  spList <- read.csv(paste(idir, "/speciesRichness.csv", sep=""))
  names(spList)[1] <- "TAXON"
  
  spList_buffer=spList[spList$IS_VALID==1,]
  results_1=list()
  if(length(spList_buffer$TAXON)==0){
    cat("None species have valid niche models \n")
  }else{
    for(spp in spList_buffer$TAXON){
      pos=which(spList_buffer$TAXON==spp)
      sppFolder <- paste(bdir, "/clim-chg/", spp, sep="")
      projFolder <- paste(sppFolder, "/maxent-projections/", sep="")
      pagrid <- paste(projFolder,"/",spp,"_",prj,"-PA.tif",sep="")
      results_1[[pos]] <- raster(pagrid)
    }
    
    cat("Writing richness raster \n")
    if(length(results_1)==0){
      cat("No valid distribution models to map \n")
    }else if(length(results_1)==1){
      results_sum_1=results_1[[1]]
    }else{
      results_sum_1=results_1[[1]]
      for(i in 2:length(results_1)){
        #       results_sum_1=extend(results_sum_1, mask) # New line
        #       results_1[[i]] = extend(results_1[[i]], mask) # New line
        results_sum_1=sum(results_sum_1,results_1[[i]], na.rm=T)
      }
    }
    results_sum <- results_sum_1
    writeRaster(results_sum,paste(idir,"/",crop,"-",prj,".tif",sep=""))
    cat("Done! \n")
  }
}

crop <- "bean"
main.dir <- paste("/curie_data/ncastaneda/threats/",crop,sep="")
bdir <- main.dir

for(crop in cropList){
  main.dir <- paste("/curie_data/ncastaneda/threats/",crop,sep="")
  bdir <- main.dir
  prjList <- c("current","rcp4_5")
  for(prj in prjList){
    cat("Doing richness map for ",crop,prj,"\n")
    x <- speciesRichness_alt(crop,bdir,prj)
  }
}

# Analyzing change

require(raster)

for(crop in cropList){
  sum.dir <- paste("/curie_data/ncastaneda/threats/",crop,"/clim-chg/_summaries", sep="")
  sppList <- read.csv(paste(sum.dir,"/speciesRichness.csv",sep=""))
  names(sppList)[1] <- "TAXON"
  
  sppList_buffer=sppList[sppList$IS_VALID==1,]
  
  for(spp in sppList_buffer$TAXON){
    cat("Analyzing change for",spp,"\n")
    in.dir <- paste("/curie_data/ncastaneda/threats/",crop,"/clim-chg/",spp, sep="")
    rs <- paste(in.dir,"/maxent-projections/",spp,"_clim-chng-PA.tif",sep="")
    rs <- raster(rs)
    rs <- as.data.frame(rs)
    rs <- table(rs)
    rs <- as.data.frame(rs)
    rs$TAXON <- spp
    rs$CROP <- crop
    names(rs)[1] <- "CAT"
    names(rs)[2] <- "FREQ"
    rs$SHARE <- NA
    
    one <- rs$FREQ[which(rs$CAT ==1)]
    two <- rs$FREQ[which(rs$CAT ==2)]
    thr <- rs$FREQ[which(rs$CAT ==3)]
    
    rs$SHARE[which(rs$CAT == 1)] <- -1*(one/sum(one,thr))
    rs$SHARE[which(rs$CAT == 2)] <- two/sum(one,thr)
    rs$SHARE[which(rs$CAT == 3)] <- thr/sum(one,thr)
    
    out.file <- write.csv(rs,paste(in.dir,"/maxent-metrics/clmChange.csv",sep=""),row.names=F)
    rm(rs)
  }
}

# Adding-up all analyses of change per genepool

for(crop in cropList){
  sum.dir <- paste("/curie_data/ncastaneda/threats/",crop,"/clim-chg/_summaries", sep="")
  sppList <- read.csv(paste(sum.dir,"/speciesRichness.csv",sep=""))
  names(sppList)[1] <- "TAXON"
  
  sppList_buffer=sppList[sppList$IS_VALID==1,]
  
  i <- 1
  for(spp in sppList_buffer$TAXON){
    cat("Summing  clmChange.csv for",crop,"\n")
    in.dir <- paste("/curie_data/ncastaneda/threats/",crop,"/clim-chg/",spp, sep="")
    in.file <- read.csv(paste(in.dir,"/maxent-metrics/clmChange.csv",sep=""))
    
    if(i==1){
      allChng <- in.file
    }else{
      allChng <- rbind(allChng,in.file)
    }
    i <- i+ 1
  }
  outFile <- paste(sum.dir,"/allclmChange.csv",sep="")
  write.csv(allChng,outFile, row.names=F)
}


# all genepools in a single file
sum.dir <- "/curie_data/ncastaneda/threats/_allCrops"
i <- 1
for(crop in cropList){
  cat("Summing  clmChange.csv for",crop,"\n")
  in.dir <- paste("/curie_data/ncastaneda/threats/",crop,"/clim-chg/_summaries", sep="")
  in.file <- read.csv(paste(in.dir,"/allclmChange.csv",sep=""))
  
  if(i==1){
    allChng <- in.file
  }else{
    allChng <- rbind(allChng,in.file)
  }
  i <- i+ 1
}
outFile <- paste(sum.dir,"/allclmChange.csv",sep="")
write.csv(allChng,outFile, row.names=F)

# Preparing graphs - MODERATE MIGRATION

allcrops <- read.csv("clipboard",header=T,sep="\t") # D:\PhD\Chapter6\allclmChange.xlsx filtered by CAT = 1 and 2

levels(allcrops$CROP) <- c("Oat","Bambara groundnut", "Bean", "Pigeonpea", "Cowpea", "Carrot",
                           "Eggplant", "Finger millet", "All crops","Sunflower", "Barley","Sweet potato",
                           "Grasspea", "Lentil","Lima bean","Apple","Alfalfa", "Banana and plantain",
                           "Pearl millet", "Pea","Potato","Rice","Rye","Sorghum", "Wheat","Vetch") # changing names

allcrops$CROP <- factor(allcrops$CROP, levels = rev(c("Finger millet","Cowpea","Potato","Carrot","Pearl millet","Oat","Lentil","Vetch", "Eggplant","Rye",
                                                  "Pea","All crops", "Grasspea","Banana and plantain","Sweet potato", "Bambara groundnut","Barley",
                                                  "Wheat","Sunflower","Alfalfa","Apple","Lima bean","Bean","Sorghum", "Rice", "Pigeonpea")), ordered=T)

require(ggplot2)
# p <- ggplot(allcrops, aes(x=reorder(factor(CROP),-SHARE.mod, FUN=median, na.rm=T), y=SHARE.mod)) #ordering per median (inverse order)
p <- ggplot(allcrops, aes(x=factor(CROP),y=SHARE.mod))
p <- p + geom_boxplot(fill = "grey80")
p <- p + theme_bw()
p <- p + coord_flip()
p <- p + labs(y="Overall suitability change (%)",x="")
p <- p + theme(text = element_text(size=20))
p <- p + geom_hline(linetype="longdash", colour="red",size=1)

p

outdir <- "D:/PhD/Chapter6"
ggsave(p,filename="allcropsclimchng-modmigr.tiff",width=10,height=8, path=outdir,dpi=300, compression="lzw") # final dir

# Preparing graphs - NO MIGRATION

allcrops.nomigr <- allcrops[which(allcrops$CAT==1),]

p <- ggplot(allcrops.nomigr, aes(x=reorder(factor(CROP),-SHARE.mod, FUN=median, na.rm=T), y=SHARE.mod)) #ordering per median (inverse order)
p <- p + geom_boxplot(fill = "grey80")
p <- p + theme_bw()
p <- p + coord_flip()
p <- p + labs(y="Overall suitability change (%)",x="")
p <- p + theme(text = element_text(size=20))
p <- p + geom_hline(linetype="longdash", colour="red",size=1)
p <- p + ylim(-100,0)
p <- p + scale_y_continuous(breaks=seq(-100,0,-50))
p

ggsave(p,filename="allcropsclimchng-nomigr.tiff",width=6,height=8, path=outdir,dpi=300, compression="lzw") # final dir
