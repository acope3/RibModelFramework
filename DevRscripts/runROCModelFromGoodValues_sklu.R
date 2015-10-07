rm(list=ls())
library(ribModel)
  
with.phi <- FALSE
  
if (with.phi) {
  genome <- initializeGenomeObject(file = "../ribModel/data/simulatedAllUniqueR.fasta", expression.file = "../ribModel/data/simulatedAllUniqueR_phi.csv")
  #genome <- initializeGenomeObject(file = "../ribModel/data/simulatedOneMix.fasta", expression.file = "../ribModel/data/simulatedOneMix_simphi.csv")
} else {
  genome <- initializeGenomeObject(file = "../ribModel/data/realGenomes/Skluyveri.fasta")
  #genome <- initializeGenomeObject(file = "../ribModel/data/twoMixtures/simulatedAllUniqueR.fasta")
  #genome <- initializeGenomeObject(file = "../ribModel/data/genome_2000.fasta")
  #genome <- initializeGenomeObject(file = "../ribModel/data/twoMixtures/simulatedKluyveri_full.fasta")
  #genome <- initializeGenomeObject(file = "../../../organisms/human/data/human_genome_cds_brain.fasta")
}
 
#initialize parameter object
sphi_init <- c(1, 1)
numMixtures <- 2
mixDef <- "allUnique"
#geneAssignment <- c(rep(1,4860), rep(2,457)) # simulated Skluyveri
geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457), rep(1, 3903)) # S.kluyveri full genome
#geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457))
#geneAssignment <- rep(1,328) # human brain
#geneAssignment <- c(rep(1,500), rep(2,500))
#geneAssignment <- rep(1,2000)
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)

#parameter <- initializeParameterObject(restart.file = "2000restartFile.rst")

phivals <- parameter$readPhiValues( "../ribModel/data/Skluyveri_GSM552569.csv")
parameter$initializeSynthesisRateByRandom(phivals)
parameter$initMutationCategories(c("../ribModel/data/Skluyveri_mutation_ChrA.csv", "../ribModel/data/Skluyveri_mutation_ChrCleft.csv") , 2)
parameter$initSelectionCategories(c("../ribModel/data/Skluyveri_selection_ChrA.csv", "../ribModel/data/Skluyveri_selection_ChrCleft.csv") , 2)
# initialize MCMC object
samples <- 1000
thining <- 10
adaptiveWidth <- 10
divergence.iteration <- 0
mcmc <- initializeMCMCObject(samples, thining, adaptive.width=adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
# get model object
model <- initializeModelObject(parameter, "ROC", with.phi = with.phi)
  
setRestartSettings(mcmc, "restartFile.rst", adaptiveWidth*200, TRUE)
#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model, 4, divergence.iteration)
)

#plots log likelihood trace, possibly other mcmc diagnostics in the future

pdf("Skluyveri_Kluyveri_allUnique_startCSP_VGAM_startPhi_SCUO_adaptSphi_True.pdf", width = 11, height = 12)
plot(mcmc)
loglik.trace <- mcmc$getLogLikelihoodTrace()
acf(loglik.trace)
convergence.test(mcmc, n.samples = 500, plot=T)

# plots different aspects of trace
trace <- parameter$getTraceObject()
plot(trace, what = "MixtureProbability")
plot(trace, what = "Sphi")
plot(trace, what = "Mphi")
if (with.phi) {
  plot(trace, what = "Aphi")
  plot(trace, what = "Sepsilon")
}
plot(trace, what = "ExpectedPhi")

mixtureAssignment <- unlist(lapply(1:genome$getGenomeSize(),  function(geneIndex){parameter$getEstimatedMixtureAssignmentForGene(samples, geneIndex)}))
expressionValues <- unlist(lapply(1:genome$getGenomeSize(), function(geneIndex){
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixtureAssignment[geneIndex])
  parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples, geneIndex, expressionCategory)
}))
expressionValues <- log10(expressionValues)
obs.phi <- log10(read.table("../ribModel/data/realGenomes/Skluyveri_phi.csv", sep=",", header=T)[, 2])
#obs.phi <- log10(read.table("../ribModel/data/twoMixtures/simulatedAllUniqueR_phi.csv", sep=",", header=T)[, 2])
#obs.phi <- log10(read.table("../ribModel/data/genome_2000.phi.csv", sep=",", header=T)[, 2])
plot(NULL, NULL, xlim=range(obs.phi) + c(-0.1, 0.1), ylim=range(expressionValues, na.rm = T) + c(-0.1, 0.1), 
     main = "Synthesis Rate", xlab = "true values", ylab = "estimated values")
for(k in 1:numMixtures){
  upper.panel.plot(obs.phi[mixtureAssignment == k], expressionValues[mixtureAssignment == k], col=ribModel:::.mixtureColors[k])
}
legend("topleft", legend = paste("Mixture Element", 1:numMixtures), 
       col = ribModel:::.mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")

plot(parameter, what = "Mutation")
plot(parameter, what = "Selection")
dev.off()


#plot(trace, what = "Expression", geneIndex = 999, mixture = 2)

pdf("Skluyveri_Kluyveri_allUnique_startCSP_VGAM_startPhi_SCUO_adaptSphi_True_mix1.pdf", width = 11, height = 12)
mixture <- 1 
plot(trace, what = "Mutation", mixture = mixture)
plot(trace, what = "Selection", mixture = mixture)

# plots model fit (cub plot)
plot(model, genome, parameter, samples = samples*0.1, mixture = mixture, main = "Codon Usage Plot")


names.aa <- aminoAcids()
selection <- c()
mutation <- c()
codon.storage <- c()
csp.m <- read.table("../ribModel/data/realGenomes/Skluyveri_mutation_ChrA.csv", sep=",", header=T)
csp.e <- read.table("../ribModel/data/realGenomes/Skluyveri_selection_ChrA.csv", sep=",", header=T)
#csp.m <- read.table("../ribModel/data/twoMixtures/simulated_mutation0.csv", sep=",", header=T)
#csp.e <- read.table("../ribModel/data/twoMixtures/simulated_selection0.csv", sep=",", header=T)
#csp.m <- read.table("../ribModel/data/simulatedOneMix_mutation.csv", sep=",", header=T)
#csp.e <- read.table("../ribModel/data/simulatedOneMix_selection.csv", sep=",", header=T)

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



pdf("Skluyveri_Kluyveri_allUnique_startCSP_VGAM_startPhi_SCUO_adaptSphi_True_mix2.pdf", width = 11, height = 12)
mixture <- 2
plot(trace, what = "Mutation", mixture = mixture)
plot(trace, what = "Selection", mixture = mixture)

# plots model fit (cub plot)
plot(model, genome, parameter, samples = samples*0.1, mixture = mixture, main = "Codon Usage Plot")


names.aa <- aminoAcids()
selection <- c()
mutation <- c()
csp.m <- read.table("../ribModel/data/realGenomes/Skluyveri_mutation_ChrCleft.csv", sep=",", header=T)
csp.e <- read.table("../ribModel/data/realGenomes/Skluyveri_selection_ChrCleft.csv", sep=",", header=T)
#csp.m <- read.table("../ribModel/data/twoMixtures/simulated_mutation1.csv", sep=",", header=T)
#csp.e <- read.table("../ribModel/data/twoMixtures/simulated_selection1.csv", sep=",", header=T)
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
