library(ribModel)
rm(list=ls())


#read genome
genome <- initializeGenomeObject(file = "../data/rfp/rfp.counts.by.codon.and.gene.GSE63789.wt.csv", FALSE)


#initialize parameter object
sphi_init <- 2
numMixtures <- 1
mixDef <- "allUnique"
geneAssignment <- c(rep(1, genome$getGenomeSize()))
#parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, model= "RFP", split.serine = TRUE, mixture.definition = mixDef)
parameter <- new(RFPParameter, "500KrestartFile.rst")

# initialize MCMC object
samples <- 100
thining <- 10
adaptiveWidth <- 10
mcmc <- initializeMCMCObject(samples=samples, thining=thining, adaptive.width=adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)


# get model object
model <- initializeModelObject(parameter, "RFP")
setRestartSettings(mcmc, "restartFile.rst", 100, TRUE)


#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model, 8)
)



# plots different aspects of trace
trace <- parameter$getTraceObject()


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
  phiList[geneIndex] <- parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(100, geneIndex, 1)
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


