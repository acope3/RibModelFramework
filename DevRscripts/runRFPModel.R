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

alphaValues <- numeric(61)
lambdaPrimeValues <- numeric(61)
for (i in 1:61)
{
  codon <- codonList[i]
  alphaValues[i] <- parameter$getAlphaPosteriorMeanForCodon(1, 30, codon)
  lambdaPrimeValues[i] <- parameter$getLambdaPrimePosteriorMeanForCodon(1, 30, codon)
}

m <- matrix(c(codonList[-c(62,63,64)], alphaValues), ncol = 2, byrow = FALSE)
colnames(m) <- c("Codon", "Alpha")
write.table(m, "RFPAlphaValues.csv", sep = ",", quote = F, row.names = F, col.names = T)

m <- matrix(c(codonList[-c(62,63,64)], lambdaPrimeValues), ncol = 2, byrow = FALSE)
colnames(m) <- c("Codon", "LambdaPrime")
write.table(m, "RFPLambdaPrimeValues.csv", sep = ",", quote = F, row.names = F, col.names = T)

phis <- parameter$getCurrentSynthesisRateForMixture(1)
m <- matrix(c(rep("void", genome$getGenomeSize()), phis, phis), ncol = 3, byrow = FALSE)
colnames(m) <- c("Gene", "PhiValue", "PhiValue")
write.table(m, "RFPPhiValues.csv", sep = ",", quote = F, row.names = F, col.names = T)
