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

plot(mcmc)

# plots different aspects of trace
trace <- parameter$getTraceObject()
pdf("RFP_Genome_allUnique_startCSP_True_startPhi_true_adaptSphi_True.pdf")
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
alphaList <- list()
lambdaPrimeList <- list()
codonList <- codons()
i <- 1
for (i in 1:61)
{
  codon <- codonList[i]
  alphaList[i] <- parameter$getParameterForCategory(cat, 0, codon, FALSE)
  lambdaPrimeList[i] <- parameter$getParameterForCategory(cat, 1, codon, FALSE)
}

reg <- lm(as.matrix(alphaList)~as.matrix(lambdaPrimeList))
plot(alphaList, lambdaPrimeList,
     main="correlationBetweenAlphaAndLambdaPrime", 
     sub=paste("r^2 = ", format(summary(reg)$r.squared, digits=3), sep=""))
abline(reg, col="red", lwd=2)
dev.off()
