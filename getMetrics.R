#To get and store in a csv (comma separated values) file the set of evaluation metrics and the thresholds that are likely to be used for model development. This is specific for the case when you have presence-only data, and have used a cross-validation approach, plus a model building using all the data.

cat("Metrics... \n")

crossValMxrFile <- paste(crossValDir, "/maxentResults.csv", sep="")
crossValMxrData <- read.csv(crossValMxrFile)
rm(crossValMxrFile)

nFolds <- 5
foldSuffix <- 'species'
FoldsinFile <- nrow(crossValMxrData) - 1

if (nFolds != FoldsinFile) {
  nFolds <- FoldsinFile
}
rm(FoldsinFile)

#Number of samples

totSamples <- crossValMxrData$X.Training.samples[1] + crossValMxrData$X.Test.samples[1]
trainSamples <- crossValMxrData$X.Training.samples[nFolds+1]
testSamples <- crossValMxrData$X.Test.samples[nFolds+1]

#AUC

trainAUCAvg <- crossValMxrData$Training.AUC[nFolds+1]
trainAUCStd <- sd(crossValMxrData$Training.AUC[1:nFolds])
testAUCAvg <- crossValMxrData$Test.AUC[nFolds+1]
testAUCStd <- sd(crossValMxrData$Test.AUC[1:nFolds])

someMets <- matrix(ncol=11, nrow=(nFolds))
nFolds <- nFolds - 1

for (fold in 0:nFolds) {
  
  #Loading basic files
  
  backPredFile <- paste(crossValDir, "/", foldSuffix, "_", fold, ".csv", sep="")
  backPredData <- read.csv(backPredFile)
  
  samPredFile <- paste(crossValDir, "/", foldSuffix, "_", fold, "_samplePredictions.csv", sep="")
  samPredData <- read.csv(samPredFile)
  
  omRatesFile <- paste(crossValDir, "/", foldSuffix, "_", fold, "_omission.csv", sep="")
  omRatesData <- read.csv(omRatesFile)
  
  trainPreds <- as.numeric(as.character(samPredData$Logistic.prediction[which(samPredData$Test.or.train == "train")]))
  matchTrain <- rep(1, times=length(trainPreds))
  
  testPreds <- as.numeric(as.character(samPredData$Logistic.prediction[which(samPredData$Test.or.train == "test")]))
  matchTest <- rep(1, times=length(testPreds))
  
  allPreds <- c(trainPreds, testPreds)
  matchAll <- rep(1, times=length(allPreds))
  
  backPreds <- backPredData[,3]
  matchBack <- rep(0.5, times=length(backPreds))
  
  #Correlation coefficients
  
  someMets[fold+1,1] <- cor(x=c(trainPreds, backPreds), y=c(matchTrain, matchBack))
  someMets[fold+1,2] <- cor(x=c(testPreds, backPreds), y=c(matchTest, matchBack))
  someMets[fold+1,3] <- cor(x=c(allPreds, backPreds), y=c(matchAll, matchBack))
  
  #LogDev
  
  someMets[fold+1,4] <- mean(-2 * log(trainPreds))
  someMets[fold+1,5] <- mean(-2 * log(testPreds))
  someMets[fold+1,6] <- mean(-2 * log(allPreds))
  
  #RMSQD
  
  someMets[fold+1,7] <- (mean((trainPreds - matchTrain)^2))^0.5
  someMets[fold+1,8] <- (mean((testPreds - matchTest)^2))^0.5
  someMets[fold+1,9] <- (mean((allPreds - matchAll)^2))^0.5
  
  #Prevalence (avg. over presence sites) threshold
  someMets[fold+1,10] <- mean(trainPreds); rm(trainPreds)
  
  #Upper-left corner ROC curve threshold
  
  Spec <- omRatesData$Fractional.area
  Sens <- 1-omRatesData$Training.omission
  logVals <- omRatesData$Corresponding.logistic.value
  
  absDif <- abs(1-(Sens+Spec)); rm(Spec); rm(Sens)
  
  if (length(logVals[which(absDif[] == min(absDif))]) > 1) {
    someMets[fold+1,11] <- mean(logVals[which(absDif[] == min(absDif))])
  } else {
    someMets[fold+1,11] <- logVals[which(absDif[] == min(absDif))]; rm(logVals); rm(absDif)
  }
}
rm(fold); rm(backPredFile); rm(samPredFile); rm(omRatesFile); rm(backPreds); rm(matchBack)
rm(matchTrain); rm(testPreds); rm(matchTest); rm(allPreds); rm(matchAll)
nFolds <- nFolds + 1

#Averaging the metrics

trainRAvg <- mean(someMets[,1])
trainRStd <- sd(someMets[,1])
testRAvg <- mean(someMets[,2])
testRStd <- sd(someMets[,2])
allRAvg <- mean(someMets[,3])
allRStd <- sd(someMets[,3])

trainLDAvg <- mean(someMets[,4])
trainLDStd <- sd(someMets[,4])
testLDAvg <- mean(someMets[,5])
testLDStd <- sd(someMets[,5])
allLDAvg <- mean(someMets[,6])
allLDStd <- sd(someMets[,6])

trainRMSQDAvg <- mean(someMets[,7])
trainRMSQDStd <- sd(someMets[,7])
testRMSQDAvg <- mean(someMets[,8])
testRMSQDStd <- sd(someMets[,8])
allRMSQDAvg <- mean(someMets[,9])
allRMSQDStd <- sd(someMets[,9])

#The metrics matrix and output file

metMatrix <- as.data.frame(matrix(ncol=25, nrow=1))
names(metMatrix) <- c("NSamples", "TrainSamples", "TestSamples", "TrainAUC", "TrainAUCSD", "TestAUC", "TestAUCSD", "TrainR", "TrainRSD", "TestR", "TestRSD", "AllR", "AllRSD", "TrainLogD", "TrainLogDSD", "TestLogD", "TestLogDSD", "AllLogD", "AllLogDSD", "TrainRMSQD", "TrainRMSQDSD", "TestRMSQD", "TestRMSQDSD", "AllRMSQD", "AllRMSQDSD")

metMatrix[1,1] <- totSamples; rm(totSamples)
metMatrix[1,2] <- trainSamples; rm(trainSamples)
metMatrix[1,3] <- testSamples; rm(testSamples)
metMatrix[1,4] <- trainAUCAvg; rm(trainAUCAvg)
metMatrix[1,5] <- trainAUCStd; rm(trainAUCStd)
metMatrix[1,6] <- testAUCAvg; rm(testAUCAvg)
metMatrix[1,7] <- testAUCStd; rm(testAUCStd)
metMatrix[1,8] <- trainRAvg; rm(trainRAvg)
metMatrix[1,9] <- trainRStd; rm(trainRStd)
metMatrix[1,10] <- testRAvg; rm(testRAvg)
metMatrix[1,11] <- testRStd; rm(testRStd)
metMatrix[1,12] <- allRAvg; rm(allRAvg)
metMatrix[1,13] <- allRStd; rm(allRStd)
metMatrix[1,14] <- trainLDAvg; rm(trainLDAvg)
metMatrix[1,15] <- trainLDStd; rm(trainLDStd)
metMatrix[1,16] <- testLDAvg; rm(testLDAvg)
metMatrix[1,17] <- testLDStd; rm(testLDStd)
metMatrix[1,18] <- allLDAvg; rm(allLDAvg)
metMatrix[1,19] <- allLDStd; rm(allLDStd)
metMatrix[1,20] <- trainRMSQDAvg; rm(trainRMSQDAvg)
metMatrix[1,21] <- trainRMSQDStd; rm(trainRMSQDStd)
metMatrix[1,22] <- testRMSQDAvg; rm(testRMSQDAvg)
metMatrix[1,23] <- testRMSQDStd; rm(testRMSQDStd)
metMatrix[1,24] <- allRMSQDAvg; rm(allRMSQDAvg)
metMatrix[1,25] <- allRMSQDStd; rm(allRMSQDStd)

outMetsFile <- paste(metricsDir, "/metrics.csv", sep="")
out <- write.csv(metMatrix, outMetsFile, quote=F, row.names=F); rm(metMatrix); rm(outMetsFile)

#Ten percentile threshold
tpThreshAvg <- crossValMxrData$X10.percentile.training.presence.logistic.threshold[nFolds+1]
tpThreshStd <- sd(crossValMxrData$X10.percentile.training.presence.logistic.threshold[1:nFolds])

#Prevalence threshold
prThreshAvg <- mean(someMets[,10])
prThreshStd <- sd(someMets[,10])

#Fixed threshold
fxThresh <- 0.5

#Maximum training sensitivity plus specificity
mspsThreshAvg <- mean(crossValMxrData$Maximum.training.sensitivity.plus.specificity.logistic.threshold[1:nFolds])
mspsThreshStd <- sd(crossValMxrData$Maximum.training.sensitivity.plus.specificity.logistic.threshold[1:nFolds])

#Equal training sensitivity and specificity
esasThreshAvg <- mean(crossValMxrData$Equal.training.sensitivity.and.specificity.logistic.threshold[1:nFolds])
esasThreshStd <- sd(crossValMxrData$Equal.training.sensitivity.and.specificity.logistic.threshold[1:nFolds])

#Balance training omission
btoThreshAvg <- mean(crossValMxrData$Balance.training.omission..predicted.area.and.threshold.value.logistic.threshold[1:nFolds])
btoThreshStd <- sd(crossValMxrData$Balance.training.omission..predicted.area.and.threshold.value.logistic.threshold[1:nFolds])

rm(crossValMxrData)
rm(nFolds)

#ROC Curve based threshold (upper left corner prob.)
rcbThreshAvg <- mean(someMets[,11])
rcbThreshStd <- sd(someMets[,11]); rm(someMets)

#The thresholds matrix and output file

threshMatrix <- as.data.frame(matrix(ncol=13, nrow=1))
names(threshMatrix) <- c("TenPercentile", "TenPercentileSD", "Prevalence", "PrevalenceSD", "FixedValue", "MaxTrainSensSpec", "MaxTrainSesnsSpecSD", "EqualTrainSensSpec", "EqualTrainSensSpecSD", "BalanceTrainOmission", "BalanceTrainOmissionSD", "UpperLeftROC", "UpperLeftROCSD")

threshMatrix[1,1] <- tpThreshAvg; rm(tpThreshAvg)
threshMatrix[1,2] <- tpThreshStd; rm(tpThreshStd)
threshMatrix[1,3] <- prThreshAvg; rm(prThreshAvg)
threshMatrix[1,4] <- prThreshStd; rm(prThreshStd)
threshMatrix[1,5] <- fxThresh; rm(fxThresh)
threshMatrix[1,6] <- mspsThreshAvg; rm(mspsThreshAvg)
threshMatrix[1,7] <- mspsThreshStd; rm(mspsThreshStd)
threshMatrix[1,8] <- esasThreshAvg; rm(esasThreshAvg)
threshMatrix[1,9] <- esasThreshStd; rm(esasThreshStd)
threshMatrix[1,10] <- btoThreshAvg; rm(btoThreshAvg)
threshMatrix[1,11] <- btoThreshStd; rm(btoThreshStd)
threshMatrix[1,12] <- rcbThreshAvg; rm(rcbThreshAvg)
threshMatrix[1,13] <- rcbThreshStd; rm(rcbThreshStd)

outThresholdFile <- paste(metricsDir, "/thresholds.csv", sep="")
write.csv(threshMatrix, outThresholdFile, quote=F, row.names=F); rm(threshMatrix); rm(outThresholdFile)
