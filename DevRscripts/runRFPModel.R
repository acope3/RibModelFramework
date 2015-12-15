rm(list=ls())
library(ribModel)



#read genome
genome <- initializeGenomeObject(file = "../data/rfp/rfp.counts.by.codon.and.gene.GSE63789.wt.csv", FALSE)


#initialize parameter object
sphi_init <- 2
numMixtures <- 1
mixDef <- "allUnique"
geneAssignment <- c(rep(1, genome$getGenomeSize()))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, model= "RFP", split.serine = TRUE, mixture.definition = mixDef)
#parameter <- new(RFPParameter, "30restartFile.rst")

# initialize MCMC object
samples <- 100
thining <- 10
adaptiveWidth <- 10
mcmc <- initializeMCMCObject(samples=samples, thining=thining, adaptive.width=adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)


# get model object
model <- initializeModelObject(parameter, "RFP")
setRestartSettings(mcmc, "restartFile.rst", 10, TRUE)


#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model, 8)
)










# plots different aspects of trace
trace <- parameter$getTraceObject()
writeParameterObject(parameter, file="RFPObject.Rdat")
writeMCMCObject(mcmc, file="MCMCObject.Rdat")


pdf("RFP_Genome_allUnique_startCSP_True_startPhi_true_adaptSphi_True.pdf")


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



pdf("RFP_CSP_Values_Mixture1.pdf", width = 11, height = 20)
#plot(trace, what = "Expression", geneIndex = 905) #used to make sure gene stabalized, not really needed now
plot(trace, what = "Alpha", mixture = 1)
plot(trace, what = "LambdaPrime", mixture = 1)
dev.off()







pdf("correlationBetweenAlphaAndLambdaPrime.pdf")
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





#Will write csv files based off posterior for alpha, lambda prime, and phi
m <- matrix(c(codonList[-c(62,63,64)], alphaList), ncol = 2, byrow = FALSE)
colnames(m) <- c("Codon", "Alpha")
write.table(m, "RFPAlphaValues.csv", sep = ",", quote = F, row.names = F, col.names = T)

m <- matrix(c(codonList[-c(62,63,64)], lambdaPrimeList), ncol = 2, byrow = FALSE)
colnames(m) <- c("Codon", "LambdaPrime")
write.table(m, "RFPLambdaPrimeValues.csv", sep = ",", quote = F, row.names = F, col.names = T)


m <- matrix(c(ids, phiList, phiList), ncol = 3, byrow = FALSE)
colnames(m) <- c("Gene", "PhiValue", "PhiValue")
write.table(m, "RFPPhiValues.csv", sep = ",", quote = F, row.names = F, col.names = T)











load("RFPtraces1.Rdat")
curloglikeTrace <- mcmc$getLogLikelihoodTrace()
cursPhiTraces <- trace$getSphiTraces()
cursphiAcceptRatTrace <- trace$getSphiAcceptanceRatioTrace()
cursynthRateTrace <- trace$getSynthesisRateTrace()
cursynthAcceptRatTrace <- trace$getSynthesisRateAcceptanceRatioTrace()
curmixAssignTrace <- trace$getMixutreAssignmentTrace()
curmixProbTrace <- trace$getMixtureProbabilitiesTrace()
curcspAcceptRatTrace <- trace$getCspAcceptanceRatioTrace()
curalphaTrace <- trace$getAlphaParameterTrace()
curlambdaPrimeTrace <- trace$getLambdaPrimeParameterTrace()

max <- samples + 1
fullAlphaTrace <- alphaTrace
for (size in 1:length(curalphaTrace))
{
  for (csp in 1:length(curalphaTrace[[size]]))
  {
    fullAlphaTrace[[size]][[csp]] <- c(fullAlphaTrace[[size]][[csp]], curalphaTrace[[size]][[csp]][2:max])
  }
}

fullLambdaPrimeTrace <- lambdaPrimeTrace
for (size in 1:length(curlambdaPrimeTrace))
{
  for (csp in 1:length(curlambdaPrimeTrace[[size]]))
  {
    fullLambdaPrimeTrace[[size]][[csp]] <- c(fullLambdaPrimeTrace[[size]][[csp]], curlambdaPrimeTrace[[size]][[csp]][2:max])
  }
}

fullcspAcceptRatTrace <- cspAcceptRatTrace
for (size in 1:length(curcspAcceptRatTrace))
{
  fullcspAcceptRatTrace[[size]]<- c(fullcspAcceptRatTrace[[size]], curcspAcceptRatTrace[[size]])
}

fullmixProbTrace <- mixProbTrace
for (size in 1:length(curmixProbTrace))
{
  fullmixProbTrace[[size]]<- c(fullmixProbTrace[[size]], curmixProbTrace[[size]][2:max])
}


fullmixAssignTrace <- mixAssignTrace
for (size in 1:length(curmixAssignTrace))
{
  fullmixAssignTrace[[size]]<- c(fullmixAssignTrace[[size]], curmixAssignTrace[[size]][2:max])
}

fullsynthAcceptRatTrace <- synthAcceptRatTrace
for (size in 1:length(cursynthAcceptRatTrace))
{
  for (csp in 1:length(cursynthAcceptRatTrace[[size]]))
  {
    fullsynthAcceptRatTrace[[size]][[csp]] <- c(fullsynthAcceptRatTrace[[size]][[csp]], cursynthAcceptRatTrace[[size]][[csp]])
  }
}


fullsynthRateTrace <- synthRateTrace
for (size in 1:length(cursynthRateTrace))
{
  for (csp in 1:length(cursynthRateTrace[[size]]))
  {
    fullsynthRateTrace[[size]][[csp]] <- c(fullsynthRateTrace[[size]][[csp]], cursynthRateTrace[[size]][[csp]][2:max])
  }
}



fullsphiAcceptRatTrace <- sphiAcceptRatTrace
fullsphiAcceptRatTrace <- c(fullsphiAcceptRatTrace, cursphiAcceptRatTrace)




fullsPhiTraces <- sPhiTraces
for (size in 1:length(cursPhiTraces))
{
  fullsPhiTraces[[size]]<- c(fullsPhiTraces[[size]], cursPhiTraces[[size]][2:max])
}

max <- max - 1
fullloglikeTrace <- loglikeTrace
fullloglikeTrace <- c(fullloglikeTrace, curloglikeTrace[2:max])




# Write the traces for the run.


# -------Sphi Traces ----------#
sPhiTraces <- trace$getSphiTraces()
m <- matrix(nrow = samples + 1, ncol = length(sPhiTraces), byrow = FALSE)
colnames(m) <- c("Mixture 1") 
for (i in 1:length(sPhiTraces)){
  m[,i] <- sPhiTraces[[i]]
}
write.table(m, "SphiTraces.csv", sep = ",", quote = F, row.names = F, col.names = T)
# -------End Sphi Traces ------#


# -------Sphi Acceptance Ratio Trace -------#
sphiAcceptRatTrace <- trace$getSphiAcceptanceRatioTrace()
m <- matrix(nrow = 1, ncol = length(sphiAcceptRatTrace), byrow = FALSE)
names <- paste(c(rep("sample", length(sphiAcceptRatTrace))), 1:length(sphiAcceptRatTrace))
colnames(m) <- names
m[1,] <- sphiAcceptRatTrace
write.table(m, "SphiAcceptanceRatioTrace.csv", sep = ",", quote = F, row.names = F, col.names = T)
# -------End Sphi Acceptance Raiot Trace ------#


# -------Synthesis Rate Trace ----------#
synthRateTrace <- trace$getSynthesisRateTrace()
for (cat in 1:length(synthRateTrace)){
  m <- matrix(nrow = samples+1, ncol = length(synthRateTrace[[cat]]), byrow = FALSE)
  names <- paste(c(rep("gene", length(synthRateTrace[[cat]]))), 1:length(synthRateTrace[[cat]]))
  colnames(m) <- names
  for (i in 1:length(synthRateTrace[[cat]])){
    m[,i] <- synthRateTrace[[cat]][[i]]
  }
  fileName <- paste("SynthesisRateTrace_ExpressionCategory", cat, "csv", sep = ".")
  write.table(m, file = fileName, sep = ",", quote = F, row.names = F, col.names = T)
}
# -------End Synthesis Rate Trace ------#


# ---------Synthesis Acceptance Ratio Trace ---------#
synthAcceptRatTrace <- trace$getSynthesisRateAcceptanceRatioTrace()
for (cat in 1:length(synthAcceptRatTrace)){
  m <- matrix(nrow = length(synthAcceptRatTrace), ncol = length(synthAcceptRatTrace[[cat]]), byrow = FALSE)
  names <- paste(c(rep("gene", length(synthAcceptRatTrace[[cat]]))), 1:length(synthAcceptRatTrace[[cat]]))
  colnames(m) <- names
  for (i in 1:length(synthAcceptRatTrace[[cat]])){
    m[,i] <- synthAcceptRatTrace[[cat]][[i]]
  }
  fileName <- paste("SynthesisRateAcceptanceRatioTrace_ExpressionCategory", cat, "csv", sep = ".")
  write.table(m, file = fileName, sep = ",", quote = F, row.names = F, col.names = T)
}
# ---------End Synthesis Acceptance Ratio Trace --------#


# --------- Mixture Assignment Trace -------#

mixAssignTrace <- trace$getmixtureP

# --------- End Mixture Assignment Trace -------#
groupList <- parameter$getGroupList()
m <- matrix(nrow = samples + 1, ncol = length(groupList))

colnames(m) <- c(groupList)
for (codon in groupList){
  tmp <- trace$getAlphaParameterTraceByMixtureElementForCodon(1, codon)
  m[,codon] <- tmp
}

