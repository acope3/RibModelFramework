library(ribModel)
rm(list=ls())


#read genome
genome <- initializeGenomeObject(file = "../ribModel/data/rfp.counts.by.codon.and.gene.GSE63789.wt.csv", FALSE)


#initialize parameter object
sphi_init <- 2
numMixtures <- 1
mixDef <- "allUnique"
geneAssignment <- c(rep(1, genome$getGenomeSize()))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, model= "RFP", split.serine = TRUE, mixture.definition = mixDef)


# initialize MCMC object
samples <- 10
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


pdf("RFP_CSP_Values_Mixture1.pdf", width = 11, height = 12)
plot(trace, what = "Expression", geneIndex = 905)
plot(trace, what = "Alpha", mixture = 1)
plot(trace, what = "LambdaPrime", mixture = 1)
dev.off()


pdf("correlationBetweenAlphaAndLambdaPrime.pdf")
cat <- 1
proposal <- FALSE
alphaList <- numeric (61)
lambdaPrimeList <- numeric (61)
waitingTimes <- numeric(61)
codonList <- codons()
i <- 1
for (i in 1:61)
{
  codon <- codonList[i]
  alphaList[i] <- parameter$getParameterForCategory(cat, 0, codon, FALSE)
  lambdaPrimeList[i] <- parameter$getParameterForCategory(cat, 1, codon, FALSE)
  waitingTimes[i] <- alphaList[i] * lambdaPrimeList[i]
}
plot(NULL, NULL, xlim=range(alphaList, na.rm = T), ylim=range(lambdaPrimeList), 
     main = "Correlation Between Alpha and Lambda Prime", xlab = "alpha", ylab = "lambdaPrime")
upper.panel.plot(alphaList, lambdaPrimeList)
dev.off()


#corrolation between RFPModel and Premal's data
X <- read.table("codon.specific.translation.rates.table.csv", header = TRUE, sep =",")
X <- X[order(X[,1]) , ]

XM <- matrix(c(X[,1], X[,2]), ncol = 2, byrow = FALSE)
Y <- matrix(c(codonList[-c(62,63,64)], waitingTimes), ncol = 2, byrow = FALSE)
colnames(Y) <- c("Codon", "PausingTime")
Y <- Y[order(Y[,1]) , ]

plot(NULL, NULL, xlim=range(XM[,2], na.rm = T), ylim=range(Y[,2]), 
     main = "Correlation Between True and RFP Model Pausing Times", xlab = "True Values", ylab = "Run Values")
upper.panel.plot(XM[,2], Y[,2])
dev.off()

#Will write csv files based off posterior for alpha, lambda prime, and phi
alphaValues <- numeric(61)
lambdaPrimeValues <- numeric(61)
size <- genome$getGenomeSize()
phis <- numeric(size)
ids <- numeric(size)
for (i in 1:61)
{
  codon <- codonList[i]
  alphaValues[i] <- parameter$getAlphaPosteriorMeanForCodon(1, 5, codon)
  lambdaPrimeValues[i] <- parameter$getLambdaPrimePosteriorMeanForCodon(1, 5, codon)
}

for (i in 1:size)
{
  phis[i] <- parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(5, i, 1)
}

for (i in 1:size)
{
  g <- genome$getGeneByIndex(i)
  ids[i] <- g$id
}

m <- matrix(c(codonList[-c(62,63,64)], alphaValues), ncol = 2, byrow = FALSE)
colnames(m) <- c("Codon", "Alpha")
write.table(m, "RFPAlphaValues.csv", sep = ",", quote = F, row.names = F, col.names = T)

m <- matrix(c(codonList[-c(62,63,64)], lambdaPrimeValues), ncol = 2, byrow = FALSE)
colnames(m) <- c("Codon", "LambdaPrime")
write.table(m, "RFPLambdaPrimeValues.csv", sep = ",", quote = F, row.names = F, col.names = T)


m <- matrix(c(ids, phis, phis), ncol = 3, byrow = FALSE)
colnames(m) <- c("Gene", "PhiValue", "PhiValue")
write.table(m, "RFPPhiValues.csv", sep = ",", quote = F, row.names = F, col.names = T)


