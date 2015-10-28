library(ribModel)
rm(list=ls())


#read genome
genome <- initializeGenomeObject(file = "../data/rfp/simulatedRFPData.csv", FALSE)


#initialize parameter object
sphi_init <- 2
numMixtures <- 1
mixDef <- "allUnique"
geneAssignment <- c(rep(1, genome$getGenomeSize()))
#parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, model= "RFP", split.serine = TRUE, mixture.definition = mixDef)
parameter <- new(RFPParameter, "1000restartFile.rst")

#init from "true" values

phiVals <- parameter$readPhiValues( "../data/rfp/RFPPhiValues.csv" )
parameter$initializeSynthesisRateByList(phiVals)
parameter$initMutationSelectionCategories(c("../data/rfp/RFPAlphaValues.csv"), 1, "Alpha")
parameter$initMutationSelectionCategories(c("../data/rfp/RFPLambdaPrimeValues.csv"), 1, "LambdaPrime")


# initialize MCMC object
samples <- 100
thining <- 10
adaptiveWidth <- 10
mcmc <- initializeMCMCObject(samples=samples, thining=thining, adaptive.width=adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
# get model object
model <- initializeModelObject(parameter, "RFP")

setRestartSettings(mcmc, "restartFile.rst", adaptiveWidth, TRUE)

#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model, 8)
)


# plots different aspects of trace
trace <- parameter$getTraceObject()
pdf("RFP_Genome_allUnique_startCSP_True_startPhi_true_adaptSphi_True.pdf")
plot(mcmc)
plot(trace, what = "MixtureProbability")
plot(trace, what = "SPhi")
plot(trace, what = "ExpectedPhi")
loglik.trace <- mcmc$getLogLikelihoodTrace()
acf(loglik.trace)
dev.off()

plot(trace, what = "Expression", geneIndex = 905)
pdf("RFP_CSP_Values_Mixture1.pdf", width = 11, height = 12)
plot(trace, what = "Alpha", mixture = 1)
plot(trace, what = "LambdaPrime", mixture = 1)
dev.off()


pdf("correlationBetweenAlphaAndLambdaPrime.pdf")
cat <- 1
proposal <- FALSE
alphaList <- numeric (61)
alphaList2 <- numeric (61)
lambdaPrimeList <- numeric (61)
lambdaPrimeList2 <- numeric (61)
LambdaPrimeLis3t <- numeric (61)
lambdaPrime.ci <-numeric(61)
phiList <- numeric(genome$getGenomeSize())
codonList <- codons()
i <- 1
for (i in 1:61)
{
  codon <- codonList[i]
  alphaList[i] <- parameter$getAlphaPosteriorMeanForCodon(cat, 100, codon)
  lambdaPrimeList[i] <- parameter$getLambdaPrimePosteriorMeanForCodon(cat, 100, codon)
  lambdaPrime.ci[i] <- parameter$getLambdaPrimeVarianceForCodon(cat, 100, codon, TRUE)
}
for (geneIndex in 1:genome$getGenomeSize()) {
  phiList[geneIndex] <- parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(100, geneIndex, 1)
}
log10(phiList) -> phiList

A <- read.table("../data/rfp/RFPAlphaValues.csv", header =TRUE, sep = ",")
LP1 <- read.table("../data/rfp/RFPLambdaPrimeValues.csv", header =TRUE, sep = ",")
PHI <- log10(read.table("../data/rfp/RFPPhiValues.csv", header = TRUE, sep = ",")[,2])

for(i in length(LP)){
  LP[i] <- (1/LP[i])
}

plot(NULL, NULL, xlim=range(alphaList, na.rm = T), ylim=range(A[,2]), 
     main = "Correlation Between Initial and Simulated Alphas", xlab = "estimated alpha", ylab = "true alpha")
upper.panel.plot(alphaList, A[,2])


plot(NULL, NULL, xlim=range(lambdaPrimeList, na.rm = T), ylim=range(LP1[,2]), 
     main = "Correlation Between Initial and Simulated Lambda Primes", xlab = "estimated lambda prime", ylab = "true lambda prime")
upper.panel.plot(lambdaPrimeList, LP1[,2],sd.y = 1.96 * sqrt(lambdaPrime.ci))

plot(NULL, NULL, xlim=range(LP1[,2], na.rm = T), ylim=range(lambdaPrimeList), 
     main = "Correlation Between Initial and Simulated Lambda Primes", xlab = "true lambda prime", ylab = "estimated lambda prime")
upper.panel.plot(LP1[,2], lambdaPrimeList,sd.y = 1.96 * sqrt(lambdaPrime.ci))

plot(NULL, NULL, xlim=range(phiList, na.rm = T), ylim=range(PHI), 
     main = "Correlation Between Initial and Simulated Synthesis Rates", xlab = "true Synthesis Rate", ylab = "simulated synthesis rate")
upper.panel.plot(phiList, PHI)


dev.off()
