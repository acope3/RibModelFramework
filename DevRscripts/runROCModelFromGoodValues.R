rm(list=ls())
library(ribModel)
#read genome
genome <- initializeGenomeObject(file = "../ribModel/data/simulatedAllUniqueR.fasta")

#initialize parameter object
sphi_init <- 2
numMixtures <- 2
mixDef <- "allUnique"
#geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457), rep(1, 3903))
#geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457))
geneAssignment <- c(rep(1,500), rep(2,500))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)

#parameter <- initializeParameterObject(restart.file = "2000restartFile.rst")

phivals <- parameter$readPhiValues( "../ribModel/data/simulatedAllUniqueR_phi.csv")
parameter$initializeSynthesisRateByRandom(phivals)
parameter$initMutationCategories(c("../ribModel/data/simulated_mutation0.csv", "../ribModel/data/simulated_mutation1.csv") , 2)
parameter$initSelectionCategories(c("../ribModel/data/simulated_selection0.csv", "../ribModel/data/simulated_selection1.csv") , 2)
# initialize MCMC object
samples <- 1000
thining <- 10
adaptiveWidth <- 10
mcmc <- initializeMCMCObject(samples, thining, adaptive.width=adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
# get model object
model <- initializeModelObject(parameter, "ROC")

setRestartSettings(mcmc, "restartFile.rst", adaptiveWidth*20, TRUE)
#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model, 4)
)


#plots log likelihood trace, possibly other mcmc diagnostics in the future
pdf("simulated_Genome_allUnique_startCSP_VGAM_startPhi_SCUO_adaptSphi_True.pdf")
plot(mcmc)
loglik.trace <- mcmc$getLogLikelihoodTrace()
acf(loglik.trace)

# plots different aspects of trace
trace <- parameter$getTraceObject()
plot(trace, what = "MixtureProbability")
plot(trace, what = "SPhi")
plot(trace, what = "ExpectedPhi")

mixtureAssignment <- unlist(lapply(1:genome$getGenomeSize(),  function(geneIndex){parameter$getEstimatedMixtureAssignmentForGene(samples, geneIndex)}))
expressionValues <- unlist(lapply(1:genome$getGenomeSize(), function(geneIndex){
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixtureAssignment[geneIndex])
  parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples, geneIndex, expressionCategory)
}))
expressionValues <- log10(expressionValues)
obs.phi <- log10(read.table("../ribModel/data/simulatedAllUniqueR_phi.csv", sep=",", header=T)[, 2])
plot(NULL, NULL, xlim=range(expressionValues, na.rm = T) + c(-0.1, 0.1), ylim=range(obs.phi) + c(-0.1, 0.1), 
     main = "Synthesis Rate", xlab = "true values", ylab = "estimated values")
upper.panel.plot(obs.phi[mixtureAssignment == 1], expressionValues[mixtureAssignment == 1], col="black")
upper.panel.plot(obs.phi[mixtureAssignment == 2], expressionValues[mixtureAssignment == 2], col="red")
legend("topleft", legend = paste("Mixture Element", 1:numMixtures), 
       col = ribModel:::.mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")

plot(parameter, what = "Mutation")
plot(parameter, what = "Selection")
dev.off()


#plot(trace, what = "Expression", geneIndex = 999, mixture = 2)

pdf("simulated_Genome_allUnique_startCSP_VGAM_startPhi_SCUO_adaptSphi_True_mix1.pdf", width = 11, height = 12)
mixture <- 1
plot(trace, what = "Mutation", mixture = mixture)
plot(trace, what = "Selection", mixture = mixture)

# plots model fit (cub plot)
plot(model, genome, parameter, samples = samples*0.1, mixture = mixture, main = "Codon Usage Plot")


names.aa <- aminoAcids()
selection <- c()
mutation <- c()
codon.storage <- c()
csp.m <- read.table("../ribModel/data/simulated_mutation0.csv", sep=",", header=T)
csp.e <- read.table("../ribModel/data/simulated_selection0.csv", sep=",", header=T)
csp <- rbind(csp.m,csp.e)
idx.eta <- 41:80
idx.mu <- 1:40
for(aa in names.aa)
{
  if(aa == "M" || aa == "W" || aa == "X") next
  codons <- AAToCodon(aa, T)
  codon.storage <- c(codon.storage,codons)
  for(i in 1:length(codons))
  {
    selection <- c(selection, parameter$getSelectionPosteriorMeanForCodon(mixture, samples*0.1, codons[i]))
    mutation <- c(mutation, parameter$getMutationPosteriorMeanForCodon(mixture, samples*0.1, codons[i]))
  }
}
plot(NULL, NULL, xlim=range(csp[idx.mu, 3], na.rm = T), ylim=range(mutation), 
     main = "Mutation", xlab = "true values", ylab = "estimated values")
upper.panel.plot(csp[idx.mu, 3], mutation)
plot(NULL, NULL, xlim=range(csp[idx.eta, 3], na.rm = T), ylim=range(selection), 
     main = "Selection", xlab = "true values", ylab = "estimated values")
upper.panel.plot(csp[idx.eta, 3], selection)
dev.off()



pdf("simulated_Genome_allUnique_startCSP_VGAM_startPhi_SCUO_adaptSphi_True_mix2.pdf", width = 11, height = 12)
mixture <- 2
plot(trace, what = "Mutation", mixture = mixture)
plot(trace, what = "Selection", mixture = mixture)

# plots model fit (cub plot)
plot(model, genome, parameter, samples = samples*0.1, mixture = mixture, main = "Codon Usage Plot")


names.aa <- aminoAcids()
selection <- c()
mutation <- c()
csp.m <- read.table("../ribModel/data/simulated_mutation1.csv", sep=",", header=T)
csp.e <- read.table("../ribModel/data/simulated_selection1.csv", sep=",", header=T)
csp <- rbind(csp.m,csp.e)
idx.eta <- 41:80
idx.mu <- 1:40
for(aa in names.aa)
{
  if(aa == "M" || aa == "W" || aa == "X") next
  codons <- AAToCodon(aa, T)
  for(i in 1:length(codons))
  {
    selection <- c(selection, parameter$getSelectionPosteriorMeanForCodon(mixture, samples*0.1, codons[i]))
    mutation <- c(mutation, parameter$getMutationPosteriorMeanForCodon(mixture, samples*0.1, codons[i]))
  }
}
plot(NULL, NULL, xlim=range(csp[idx.mu, 3], na.rm = T), ylim=range(mutation), 
     main = "Mutation", xlab = "true values", ylab = "estimated values")
upper.panel.plot(csp[idx.mu, 3], mutation)
plot(NULL, NULL, xlim=range(csp[idx.eta, 3], na.rm = T), ylim=range(selection), 
     main = "Selection", xlab = "true values", ylab = "estimated values")
upper.panel.plot(csp[idx.eta, 3], selection)
dev.off()
