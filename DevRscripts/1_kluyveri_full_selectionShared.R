rm(list=ls()) 
library(ribModel) 
set.seed(2935642) 
with.phi <- FALSE

if (with.phi) { 
genome <- initializeGenomeObject(file ="../data/realGenomes/Skluyveri.fasta", expression.file = "../data/simulatedAllUniqueR_phi.csv") 
} else { 
genome <- initializeGenomeObject(file = "../data/realGenomes/Skluyveri.fasta") 
} 

sphi_init <- c(1,1)
numMixtures <- 2
#mixDef <- "allUnique"
mixDef <- "selectionShared"
geneAssignment <- c(rep(1,961), rep(2,457), rep(1, 3403), rep(1, 500))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef) 

phivals <- parameter$readPhiValues( "../data/realGenomes/Skluyveri_phi.csv")
parameter$initializeSynthesisRateByList(phivals)
parameter$initMutationCategories(c("../data/realGenomes/Skluyveri_mutation_ChrA.csv", "../data/realGenomes/Skluyveri_mutation_ChrCleft.csv") , 2)
parameter$initSelectionCategories(c("../data/realGenomes/Skluyveri_selection_ChrA.csv", "../data/realGenomes/Skluyveri_selection_ChrCleft.csv") , 1)


samples <- 100
thining <- 50
adaptiveWidth <- 1000000
divergence.iteration <- 0
mcmc <- initializeMCMCObject(samples, thining, adaptive.width=adaptiveWidth, est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE) 

setRestartSettings(mcmc, "1_kluyveri_full.rst", adaptiveWidth*50, TRUE) 

model <- initializeModelObject(parameter, "ROC", with.phi = with.phi) 

start <- Sys.time() 
system.time(runMCMC(mcmc, genome, model, 4, divergence.iteration)) 
end <- Sys.time() 
end - start 


pdf("1_kluyveri_full_global.pdf", width = 11, height = 12) 
plot(mcmc) 
loglik.trace <- mcmc$getLogLikelihoodTrace() 
acf(loglik.trace) 
convergence.test(mcmc, n.samples = 500, plot=T) 

trace <- parameter$getTraceObject() 
plot(trace, what = "MixtureProbability") 
plot(trace, what = "Sphi") 
plot(trace, what = "Mphi") 
if (with.phi) { 
plot(trace, what = "Aphi") 
plot(trace, what = "Sepsilon") 
} 
plot(trace, what = "ExpectedPhi") 

mixtureAssignment <- unlist(lapply(1:genome$getGenomeSize(),  function(geneIndex){parameter$getEstimatedMixtureAssignmentForGene(samples*0.1, geneIndex)})) 
expressionValues <- unlist(lapply(1:genome$getGenomeSize(), function(geneIndex){ 
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixtureAssignment[geneIndex]) 
  parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples*0.1, geneIndex, expressionCategory) 
})) 
expressionValues <- log10(expressionValues) 
obs.phi <- log10(read.table("../data/realGenomes/Skluyveri_phi.csv", sep=",", header=T)[, 2]) 
plot(NULL, NULL, xlim=range(obs.phi) + c(-0.1, 0.1), ylim=range(expressionValues, na.rm = T) + c(-0.1, 0.1), 
main = "Synthesis Rate", xlab = "True values", ylab = "Estimated values") 
for(k in 1:numMixtures){ 
  upper.panel.plot(obs.phi[mixtureAssignment == k], expressionValues[mixtureAssignment == k], col=ribModel:::.mixtureColors[k]) 
} 
legend("topleft", legend = paste("Mixture Element", 1:numMixtures), 
       col = ribModel:::.mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n") 

mixprob <- do.call("rbind", lapply(1:genome$getGenomeSize(), function(geneIndex){
  props <- parameter$getEstimatedMixtureAssignmentProbabilitiesForGene(samples*0.1, geneIndex)
  gene <- genome$getGeneByIndex(geneIndex, F)
  c(gene$id, props)
}))
plot(NULL, NULL, xlim = c(1-0.2, numMixtures+0.2), ylim = c(0, 1), xlab="Mixture Assignment", ylab = "Probability")
for(k in 1:numMixtures){
  points(mixtureAssignment[mixtureAssignment == k] + runif(length(mixtureAssignment[mixtureAssignment == k]), -0.2, 0.2), 
         as.numeric(mixprob[mixtureAssignment == k, k+1]), col = ribModel:::.mixtureColors[k])
}
legend("bottomleft", legend = paste("Mixture Element", 1:numMixtures),
       col = ribModel:::.mixtureColors[1:numMixtures], pch = rep(1, numMixtures), bty = "n")

write.table(file="1_kluyveri_full.csv", x = mixprob, sep = ",", row.names = F, quote = F, col.names = F)

plot(parameter, what = "Mutation") 
plot(parameter, what = "Selection") 
dev.off() 

observed.mutation <- c("../data/realGenomes/Skluyveri_mutation_ChrA.csv", "../data/realGenomes/Skluyveri_mutation_ChrCleft.csv")
observed.selection <- c("../data/realGenomes/Skluyveri_selection_ChrA.csv", "../data/realGenomes/Skluyveri_selection_ChrCleft.csv")
for(k in 1:numMixtures){ 
   pdf(paste("1_kluyveri_full_mixture_", k, ".pdf", sep=""), width = 11, height = 12) 
   mixture <- k 
   plot(trace, what = "Mutation", mixture = mixture) 
   plot(trace, what = "Selection", mixture = mixture) 
   plot(model, genome, parameter, samples = samples*0.1, mixture = mixture, main = "Codon Usage Plot") 
   names.aa <- aminoAcids() 
   selection.ci <- c() 
   mutation.ci <- c() 
   selection <- c() 
   mutation <- c() 
   codon.storage <- c() 
   csp.m <- read.table(observed.mutation[k], sep=",", header=T) 
   csp.e <- read.table(observed.selection[k], sep=",", header=T) 
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
       selection.ci <- c(selection.ci, parameter$getSelectionVarianceForCodon(mixture, samples*0.1, codons[i], TRUE)) 
       mutation <- c(mutation, parameter$getMutationPosteriorMeanForCodon(mixture, samples*0.1, codons[i])) 
       mutation.ci <- c(mutation.ci, parameter$getMutationVarianceForCodon(mixture, samples*0.1, codons[i], TRUE)) 
     } 
   } 
   plot(NULL, NULL, xlim=range(csp[idx.mu, 3], na.rm = T), ylim=range(mutation), 
        main = "Mutation", xlab = "True values", ylab = "Estimated values") 
   upper.panel.plot(x = csp[idx.mu, 3], y = mutation, sd.y = 1.96*sqrt(mutation.ci)) 
   plot(NULL, NULL, xlim=range(csp[idx.eta, 3], na.rm = T), ylim=range(selection), 
        main = "Selection", xlab = "True values", ylab = "Estimated values") 
   upper.panel.plot(x = csp[idx.eta, 3], y = selection, sd.y = 1.96*sqrt(selection.ci)) 
   dev.off() 
} 
