library(ribModel)
rm(list=ls())

parameter <- new(ROCParameter)
parameter <- loadParameterObject(parameter, "ROCParameter.Rdat", "ROC")
mcmc <- loadMCMCObject("MCMCObject.Rdat")
genome <- initializeGenomeObject(file = "../data/realGenomes/Skluyveri.fasta")



# plots different aspects of trace
trace <- parameter$getTraceObject()

pdf("test1.pdf")
plot(mcmc)
plot(trace, what = "MixtureProbability")
plot(trace, what = "SPhi")
plot(trace, what = "ExpectedPhi")
dev.off()
plot(trace, what = "Expression", geneIndex = 905)
pdf("test2.pdf", width = 11, height = 12)
plot(trace, what = "Mutation", mixture = 1)
plot(trace, what = "Selection", mixture = 1)
mixtureAssignment <- unlist(lapply(1:genome$getGenomeSize(),  function(geneIndex){parameter$getEstimatedMixtureAssignmentForGene(samples, geneIndex)}))
expressionValues <- unlist(lapply(1:genome$getGenomeSize(), function(geneIndex){
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixtureAssignment[geneIndex])
  parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples, geneIndex, expressionCategory)
}))
expressionValues <- log10(expressionValues)
obs.phi <- log10(read.table("../data//twoMixtures/simulatedAllUniqueR_phi.csv", sep=",", header=T)[, 2])
plot(NULL, NULL, xlim=range(expressionValues, na.rm = T) + c(-0.1, 0.1), ylim=range(obs.phi) + c(-0.1, 0.1), 
     main = "Synthesis Rate", xlab = "true values", ylab = "estimated values")
upper.panel.plot(obs.phi[mixtureAssignment == 1], expressionValues[mixtureAssignment == 1], col="black")
upper.panel.plot(obs.phi[mixtureAssignment == 2], expressionValues[mixtureAssignment == 2], col="red")
legend("topleft", legend = paste("Mixture Element", 1:numMixtures), 
       col = ribModel:::.mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")
# plots model fit (cub plot)
plot(model, genome, parameter, samples = samples*0.1, mixture = 1, main = "S. kluyveri Chr (A,B,Cleft) Codon Usage Plot")
dev.off()

plot(parameter, what = "Mutation", main = "Mutation Correlation, Not shared")
plot(parameter, what = "Selection", main = "Selecion Correlation, Not shared")

##-----------------------------------------##

#creates a plot Mixture Prob vs log10(phi) (Not part of the package)
num.genes <- genome$getGenomeSize()
mixtureAssignment <- do.call("rbind", lapply(1:num.genes,  function(geneIndex){parameter$getEstimatedMixtureAssignmentProbabilitiesForGene(samples*0.1, geneIndex)}))
plot(NULL, NULL, xlim = c(0,1), ylim = c(-2, 1), xlab = "Probability beeing in Mixture 2", ylab = "log10(phi)")
colors <- c("black", "red")
for(mixture in 1:2)
{
  genes.in.mixture <- which(round(mixtureAssignment[, mixture]+1) == mixture)
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixture)
  
  # need expression values to know range
  num.genes <- length(genes.in.mixture)
  expressionValues <- unlist(lapply(genes.in.mixture, function(geneIndex){
    parameter$getExpressionPosteriorMeanByMixtureElementForGene(samples*0.1, geneIndex, expressionCategory)
  }))
  points(mixtureAssignment[genes.in.mixture, mixture], log10(expressionValues), col=colors[mixture])
}