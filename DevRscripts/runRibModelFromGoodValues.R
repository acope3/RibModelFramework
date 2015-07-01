library(ribModel)
rm(list=ls())
#read genome
genome <- initializeGenomeObject(fasta.file = "../ribModel/data/simulatedAllUniqueR.fasta")

#initialize parameter object
sphi_init <- 2
numMixtures <- 2
mixDef <- "allUnique"
#geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457), rep(1, 3903))
#geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457))
geneAssignment <- c(rep(1,500), rep(2,500))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)



phivals <- parameter$readPhiValues( "../ribModel/data/simulatedAllUniqueR_phi.csv")
parameter$initializeExpressionByRandom(phivals)
parameter$initMutationSelectionCategories(c("../ribModel/data/simulated_CSP0.csv", "../ribModel/data/simulated_CSP1.csv") , 2, "Selection")
parameter$initMutationSelectionCategories(c("../ribModel/data/simulated_CSP0.csv", "../ribModel/data/simulated_CSP1.csv") , 2, "Mutation")
# initialize MCMC object
samples <- 100
thining <- 10
adaptiveWidth <- 10
mcmc <- initializeMCMCObject(samples, thining, adaptive.width=adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
# get model object
model <- initializeModelObject(parameter, "ROC")

setRestartSettings(mcmc, "restartFile.rst", adaptiveWidth, TRUE)
#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model)
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


plot(trace, what = "Expression", geneIndex = 999, mixture = 2)

pdf("simulated_Genome_allUnique_startCSP_VGAM_startPhi_SCUO_adaptSphi_True_mix2.pdf", width = 11, height = 12)
mixture <- 2
plot(trace, what = "Mutation", mixture = mixture)
plot(trace, what = "Selection", mixture = mixture)

# plots model fit (cub plot)
plot(model, genome, parameter, samples = samples*0.1, mixture = mixture, main = "Codon Usage Plot")


names.aa <- aminoAcids()
selection <- c()
mutation <- c()
csp <- read.table("../ribModel/data/simulated_CSP1.csv", sep=",", header=T)
idx.eta <- grepl(pattern = "[A-Z].[A-Z]{3}.Delta.eta", x = as.character(csp[,1]))
idx.mu <- grepl(pattern = "[A-Z].[A-Z]{3}.log.mu", x = as.character(csp[,1]))
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
plot(NULL, NULL, xlim=range(csp[idx.mu, 2], na.rm = T), ylim=range(mutation), 
     main = "Mutation", xlab = "true values", ylab = "estimated values")
upper.panel.plot(csp[idx.mu, 2], mutation)
plot(NULL, NULL, xlim=range(csp[idx.eta, 2], na.rm = T), ylim=range(selection), 
     main = "Selection", xlab = "true values", ylab = "estimated values")
upper.panel.plot(csp[idx.eta, 2], selection)
dev.off()
