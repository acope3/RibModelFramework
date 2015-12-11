library(ribModel)
rm(list=ls())
#TEST BLOCK
load("RFPObject.Rdat")
p <- new(RFPParameter)
p$setCurrentAlphaParameter(currentAlpha)
p$setProposedAlphaParameter(proposedAlpha)
p$setCurrentLambdaPrimeParameter(currentLambdaPrime)
p$setProposedLambdaPrimeParameter(proposedLambdaPrime)

t <- p$getTraceObject()
t$setSphiTraces(paramBase$sPhiTraces)
t$setSphiAcceptanceRatioTrace(paramBase$sphiAcceptRatTrace)

t$setSynthesisRateTrace(paramBase$synthRateTrace)
t$setSynthesisRateAcceptanceRatioTrace(paramBase$synthAcceptRatTrace)
t$setMixtureAssignmentTrace(paramBase$mixAssignTrace)
t$setMixtureProbabilitiesTrace(paramBase$mixProbTrace)
t$setCspAcceptanceRatioTrace(paramBase$cspAcceptRatTrace)
p$setRFPTrace(t)
p$numMixtures <- paramBase$numMix
p$numMutationCategories <- paramBase$numMut
p$numSelectionCategories <- paramBase$numSel
#END TEST BLOCK

trace <- p$getTraceObject()

pdf("test1.pdf")


plot(mcmc) #plots the whole logliklihood trace

#Here I take a subset of the trace values for the logliklihood trace and plot them.
#The primary reason for doing this is the "jump" that throws the scale of the graph
#at the beginning is removed by taking out the beginning values.
loglik.trace <- mcmc$getLogLikelihoodTrace()
start <- length(loglik.trace) * 0.5 #the multiplier determines how much of the beginning trace is 
#eliminated.

logL <- logL <- mean(loglik.trace[start:length(loglik.trace)]) #get the mean for the subset
plot(loglik.trace[start:length(loglik.trace)], type="l", main=paste("logL:", logL), xlab="Sample", ylab="log(Likelihood)")
grid (NULL,NULL, lty = 6, col = "cornsilk2")


plot(trace, what = "MixtureProbability") #right now, will be straight line (mix =1)
plot(trace, what = "Mphi")
plot(trace, what = "Sphi")
plot(trace, what = "ExpectedPhi")
loglik.trace <- mcmc$getLogLikelihoodTrace()
acf(loglik.trace)
acf(loglik.trace[start:length(loglik.trace)])
dev.off()



pdf("test2.pdf", width = 11, height = 20)
#plot(trace, what = "Expression", geneIndex = 905) #used to make sure gene stabalized, not really needed now
plot(trace, what = "Alpha", mixture = 1)
plot(trace, what = "LambdaPrime", mixture = 1)
dev.off()







pdf("test3.pdf")
cat <- 1
proposal <- FALSE
alphaList <- numeric (61)
lambdaPrimeList <- numeric (61)
waitingTimes <- numeric(61)
phiList <- numeric(genome$getGenomeSize())
ids <- numeric(genome$getGenomeSize())
codonList <- codons()
i <- 1
for (i in 1:61)
{
  codon <- codonList[i]
  alphaList[i] <- parameter$getAlphaPosteriorMeanForCodon(cat, samples * 0.5, codon)
  lambdaPrimeList[i] <- parameter$getLambdaPrimePosteriorMeanForCodon(cat, samples * 0.5, codon)
  waitingTimes[i] <- alphaList[i] * lambdaPrimeList[i]
}

for (geneIndex in 1:genome$getGenomeSize()) {
  phiList[geneIndex] <- parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples * 0.5, geneIndex, 1)
}

for (i in 1:genome$getGenomeSize())
{
  g <- genome$getGeneByIndex(i, FALSE)
  ids[i] <- g$id
}






plot(NULL, NULL, xlim=range(alphaList, na.rm = T), ylim=range(lambdaPrimeList), 
     main = "Correlation Between Alpha and Lambda Prime", xlab = "alpha", ylab = "lambdaPrime")
upper.panel.plot(alphaList, lambdaPrimeList)


#corrolation between RFPModel and Premal's data
#X <- read.table("../data/rfp/codon.specific.translation.rates.table.csv", header = TRUE, sep =",")
#X <- X[order(X[,1]) , ]

#XM <- matrix(c(X[,1], X[,2]), ncol = 2, byrow = FALSE)
#Y <- data.frame(codonList[-c(62,63,64)], waitingTimes)
#colnames(Y) <- c("Codon", "PausingTime")
#Y <- Y[order(Y[,1]) , ]

#plot(NULL, NULL, xlim=range(XM[,2], na.rm = T), ylim=range(Y[,2]), 
#     main = "Correlation Between Premal and RFP Model Pausing Times", xlab = "True Values", ylab = "Run Values")
#upper.panel.plot(XM[,2], Y[,2])
dev.off()