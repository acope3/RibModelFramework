rm(list=ls()) 
library(ribModel)
seeds <- read.table(file = "seed.txt")[,1]
task.id <- 1 #as.numeric(Sys.getenv("SGE_TASK_ID"))

set.seed(446141)
with.phi <- FALSE

if (with.phi) {
  genome <- initializeGenomeObject(file = "../data/twoMixtures/simulatedAllUniqueR.fasta", expression.file = "../data/simulatedAllUniqueR_phi.csv") 
} else {
  genome <- initializeGenomeObject(file = "../data/twoMixtures/simulatedAllUniqueR.fasta") 
}

sphi_init <- c(1,1)
numMixtures <- 2
mixDef <- "allUnique"
geneAssignment <- sample(c(1,2), size = length(genome), replace = TRUE, prob = c(0.3, 0.7)) #c(rep(1,500), rep(2,500))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)
#parameter <- initializeParameterObject(genome, sphi = sphi_init, geneAssignment = geneAssignment, numMixtures = numMixtures, restart.file = "5001_simulated_allUnique.rst")

samples <- 500
thining <- 10
adaptiveWidth <- 10000
divergence.iteration <- 0
mcmc <- initializeMCMCObject(samples, thining, adaptive.width=adaptiveWidth, est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE) 

setRestartSettings(mcmc, paste(task.id, "_simulated_allUnique.rst", sep=""), adaptiveWidth*50, TRUE) 

model <- initializeModelObject(parameter, "ROC", with.phi = with.phi) 

start <- Sys.time() 
system.time(runMCMC(mcmc, genome, model, 4, divergence.iteration)) 
end <- Sys.time() 
end - start 


pdf(paste(task.id, "_simulated_global_allUnique.pdf", sep=""), width = 11, height = 12) 
plot(mcmc)
loglik.trace <- mcmc$getLogLikelihoodTrace()[-1]
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

mixtureAssignment <- unlist(lapply(1:length(genome),  function(geneIndex){parameter$getEstimatedMixtureAssignmentForGene(samples*0.1, geneIndex)})) 
expressionValues <- unlist(lapply(1:length(genome), function(geneIndex){ 
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixtureAssignment[geneIndex]) 
  parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples*0.1, geneIndex, expressionCategory) 
})) 
expressionValues <- log10(expressionValues) 
obs.phi <- log10(read.table("../data/twoMixtures/simulatedAllUniqueR_phi.csv", sep=",", header=T)[, 2]) 
plot(NULL, NULL, xlim=range(obs.phi) + c(-0.1, 0.1), ylim=range(expressionValues, na.rm = T) + c(-0.1, 0.1), 
     main = "Synthesis Rate", xlab = "True values", ylab = "Estimated values") 
for(k in 1:numMixtures){
  if(sum(mixtureAssignment == k) == 0) next
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

write.table(file=paste(task.id, "_simulated_allUnique.csv", sep=""), x = mixprob, sep = ",", row.names = F, quote = F, col.names = F)

gene.ids <- unlist(lapply(1:genome$getGenomeSize(), function(geneIndex){ 
  gene <- genome$getGeneByIndex(geneIndex, F)
  return(gene$id)
})) 
gene.obs.ids <- as.character(read.table("../data/twoMixtures/simulatedAllUniqueR_phi.csv", sep=",", header=T)[, 1])
mixtureAssignment <- mixtureAssignment[gene.ids %in% gene.obs.ids]
expressionValues <- expressionValues[gene.ids %in% gene.obs.ids]
gene.ids <- gene.ids[gene.ids %in% gene.obs.ids]
mat <- cbind(gene.ids, expressionValues, mixtureAssignment)


write.table(file=paste(task.id, "_simulated_filtered_expression_allUnique.csv", sep=""), x = mat, sep=",", quote = F, row.names = F, col.names = F)

plot(parameter, what = "Mutation", samples = samples*0.1) 
plot(parameter, what = "Selection", samples = samples*0.1) 
dev.off() 

observed.mutation <- c("../data/twoMixtures/simulated_mutation0.csv", "../data/twoMixtures/simulated_mutation1.csv")
observed.selection <- c("../data/twoMixtures/simulated_selection0.csv", "../data/twoMixtures/simulated_selection1.csv")
for(k in 1:numMixtures){
  if(sum(mixtureAssignment == k) == 0) next
  pdf(paste(task.id, "_simulated_mixture_", k, "_allUnique.pdf", sep=""), width = 11, height = 12) 
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
