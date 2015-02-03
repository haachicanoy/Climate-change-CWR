
cajanus_mollis <- list.files("C:/Users/haachicanoy/Desktop", pattern=".asc", full.names=T)
library(raster)
models <- lapply(cajanus_mollis, raster)
models <- stack(models)
plot(models)

quantile_ensemble <- function(x, quantile)
{
  tot_ensemble <- sum(x)
  tot_ensemble[which(tot_ensemble[] < nlayers(x)*quantile)] <- NA
  tot_ensemble[which(!is.na(tot_ensemble[]))] <- 1
  return(tot_ensemble)
}

test <- quantile_ensemble(x=models, quantile=0.7)
plot(test, xlim=c(60, 110), ylim=c(0, 50))
