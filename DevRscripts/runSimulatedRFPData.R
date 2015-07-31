library(ribModel)
rm(list=ls())


#read genome
genome <- initializeGenomeObject(file = "../ribModel/data/simulatedRFPData.csv", FALSE)


#initialize parameter object
sphi_init <- 2
numMixtures <- 1
mixDef <- "allUnique"
geneAssignment <- c(rep(1, genome$getGenomeSize()))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, model= "RFP", split.serine = TRUE, mixture.definition = mixDef)


#init from "true" values

phiVals <- parameter$readPhiValues( "../ribModel/data/RFPPhiValues.csv" )
parameter$initializeSynthesisRateByRandom(phiVals)
parameter$initMutationSelectionCategories(c("../ribModel/data/RFPAlphaValues.csv"), 1, "Alpha")
parameter$initMutationSelectionCategories(c("../ribModel/data/RFPLambdaPrimeValues.csv"), 1, "LambdaPrime")


# initialize MCMC object
samples <- 1000
thining <- 10
adaptiveWidth <- 10
mcmc <- initializeMCMCObject(samples=samples, thining=thining, adaptive.width=adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
# get model object
model <- initializeModelObject(parameter, "RFP")

setRestartSettings(mcmc, "restartFile.rst", adaptiveWidth, TRUE)

#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model)
)


# plots different aspects of trace
trace <- parameter$getTraceObject()
pdf("RFP_Genome_allUnique_startCSP_True_startPhi_true_adaptSphi_True.pdf")
plot(mcmc)
plot(trace, what = "MixtureProbability")
plot(trace, what = "SPhi")
plot(trace, what = "ExpectedPhi")
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
lambdaPrimeList <- numeric (61)
codonList <- codons()
i <- 1
for (i in 1:61)
{
  codon <- codonList[i]
  alphaList[i] <- parameter$getParameterForCategory(cat, 0, codon, FALSE)
  lambdaPrimeList[i] <- parameter$getParameterForCategory(cat, 1, codon, FALSE)
}


plot(NULL, NULL, xlim=range(alphaList, na.rm = T), ylim=range(lambdaPrimeList), 
     main = "Correlation Between Alpha and Lambda Prime", xlab = "alpha", ylab = "lambdaPrime")
upper.panel.plot(alphaList, lambdaPrimeList)
dev.off()

A <- read.table("../ribModel/data/RFPAlphaValues.csv", header =TRUE, sep = ",")
LP <- read.table("../ribModel/data/RFPLambdaPrimeValues.csv", header =TRUE, sep = ",")

plot(NULL, NULL, xlim=range(alphaList, na.rm = T), ylim=range(A[,2]), 
     main = "Correlation Between Initial and Simulated Alphas", xlab = "true alpha", ylab = "simulated alpha")
upper.panel.plot(alphaList, A[,2])

plot(NULL, NULL, xlim=range(lambdaPrimeList, na.rm = T), ylim=range(LP[,2]), 
     main = "Correlation Between Initial and Simulated Lambda Primes", xlab = "true lambda prime", ylab = "simulated lambda prime")
upper.panel.plot(lambdaPrimeList, LP[,2])
dev.off()