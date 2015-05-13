# Do projections for MaxEnt runs
# H. Achicanoy
# CIAT, 2015

make.projections <- function(k)
{
  cat('\n\nProcessing fold:',k,'\n\n')
  lambdas.file <- strsplit(x=fit@models[[k]]@lambdas,split=',',fixed=TRUE)
  lambdas.file <- lapply(1:length(lambdas.file),function(i){z <- data.frame(t(lambdas.file[[i]])); return(z)})
  identify.lmd <- unlist(lapply(1:length(lambdas.file),function(i){z <- ncol(lambdas.file[[i]])==4; return(z)}))
  paramet.file <- lambdas.file[!identify.lmd]
  lambdas.file <- lambdas.file[identify.lmd]
  paramet.file <- Reduce(function(...) rbind(..., deparse.level=1), paramet.file)
  lambdas.file <- Reduce(function(...) rbind(..., deparse.level=1), lambdas.file)
  names(paramet.file) <- c("variable","value")
  names(lambdas.file) <- c("feature","lambda","min","max")
  
  lambdas.file$feature <- as.character(lambdas.file$feature)
  lambdas.file$lambda <- as.numeric(as.character(lambdas.file$lambda))
  lambdas.file$min <- as.numeric(as.character(lambdas.file$min))
  lambdas.file$max <- as.numeric(as.character(lambdas.file$max))
  
  paramet.file$variable <- as.character(paramet.file$variable)
  paramet.file$value <- as.numeric(as.character(paramet.file$value))
  
  q.feat <- grep(pattern="^",x=lambdas.file$feature,fixed=T) # Quadratic features
  p.feat <- grep(pattern="*",x=lambdas.file$feature,fixed=T) # Product features
  l.feat <- setdiff(1:nrow(lambdas.file),c(q.feat,p.feat)) # Linear features
  
  temp.dt <- climData[[1]][[1]][[1]]
  temp.dt <- ff(1,dim=c(ncell(temp.dt),20),vmode="double")
  lapply(1:20,function(i){z <- climData[[1]][[1]][[i]]; t <- getValues(z); cat('Processing: biovariable',i,'\n'); temp.dt[,i] <- t[]; return(cat("Done\n"))})
  temp.dt <- as.ffdf(temp.dt)
  names(temp.dt) <- paste0("variable.",1:20)
  temp.dt <- as.data.frame(temp.dt)
  temp.dt <- data.table(temp.dt)
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
  # Linear features
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
  
  linear.calcs <- lapply(lambdas.file$feature[l.feat],function(var)
  {
    lambdas.file.l <- lambdas.file[l.feat,]
    eval(parse(text=paste('result.by.var.l <- lambdas.file.l$lambda[which(lambdas.file.l$feature==var)]*((temp.dt[,',var,']-lambdas.file.l$min[which(lambdas.file.l$feature==var)])/(lambdas.file.l$max[which(lambdas.file.l$feature==var)]-lambdas.file.l$min[which(lambdas.file.l$feature==var)]))',sep='')))
    result.by.var.l <- data.frame(result.by.var.l)
    return(result.by.var.l)
  })
  linear.calcs <- Reduce(function(...) cbind(..., deparse.level=1), linear.calcs)
  names(linear.calcs) <- lambdas.file$feature[l.feat]
  linear.calcs <- as.data.table(linear.calcs)
  linear.calcs <- linear.calcs[,rowSums(.SD,na.rm=FALSE)]
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
  # Quadratic features
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
  
  temp.name <- strsplit(x=lambdas.file$feature[q.feat],split="^",fixed=TRUE)
  temp.name <- unlist(lapply(temp.name,function(set){return(set[1])}))
  lambdas.file$feature[q.feat] <- temp.name; rm(temp.name)
  quadratic.calcs <- lapply(lambdas.file$feature[q.feat],function(var)
  {
    lambdas.file.q <- lambdas.file[q.feat,]
    eval(parse(text=paste('result.by.var.q <- lambdas.file.q$lambda[which(lambdas.file.q$feature==var)]*((temp.dt[,',var,']^2-lambdas.file.q$min[which(lambdas.file.q$feature==var)])/(lambdas.file.q$max[which(lambdas.file.q$feature==var)]-lambdas.file.q$min[which(lambdas.file.q$feature==var)]))',sep='')))
    result.by.var.q <- data.frame(result.by.var.q)
    return(result.by.var.q)
  })
  quadratic.calcs <- Reduce(function(...) cbind(..., deparse.level=1), quadratic.calcs)
  names(quadratic.calcs) <- lambdas.file$feature[q.feat]
  quadratic.calcs <- as.data.table(quadratic.calcs)
  quadratic.calcs <- quadratic.calcs[,rowSums(.SD,na.rm=FALSE)]
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
  # Product features
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
  
  product.calcs <- lapply(lambdas.file$feature[p.feat],function(comb.var)
  {
    lambdas.file.p <- lambdas.file[p.feat,]
    eval(parse(text=paste('result.by.var.p <- lambdas.file.p$lambda[which(lambdas.file.p$feature==comb.var)]*((temp.dt[,',comb.var,']-lambdas.file.p$min[which(lambdas.file.p$feature==comb.var)])/(lambdas.file.p$max[which(lambdas.file.p$feature==comb.var)]-lambdas.file.p$min[which(lambdas.file.p$feature==comb.var)]))',sep='')))
    result.by.var.p <- data.frame(result.by.var.p)
    return(result.by.var.p)
  })
  product.calcs <- Reduce(function(...) cbind(..., deparse.level=1), product.calcs)
  names(product.calcs) <- lambdas.file$feature[p.feat]
  product.calcs <- as.data.table(product.calcs)
  product.calcs <- product.calcs[,rowSums(.SD,na.rm=FALSE)]
  
  fx.calcs <- cbind(linear.calcs,quadratic.calcs,product.calcs); rm(linear.calcs,quadratic.calcs,product.calcs)
  fx.calcs <- as.data.table(fx.calcs)
  fx.calcs <- fx.calcs[,rowSums(.SD,na.rm=FALSE)]
  
  S.x <- fx.calcs - paramet.file$value[which(paramet.file$variable=="linearPredictorNormalizer")]
  Q.x <- exp(S.x)/paramet.file$value[which(paramet.file$variable=="densityNormalizer")]
  L.x <- (Q.x*exp(paramet.file$value[which(paramet.file$variable=="entropy")]))/(1+Q.x*exp(paramet.file$value[which(paramet.file$variable=="entropy")]))
  
  prj_fn <- climData[[1]][[1]][[1]]
  cell <- 1:ncell(prj_fn)
  prj_fn[cell] <- L.x
  
  return(prj_fn)
  
}