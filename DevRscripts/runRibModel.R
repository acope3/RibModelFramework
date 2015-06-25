library(ribModel)

#read genome
genome <- initializeGenomeObject(fasta.file = "../ribModel/data/simulatedAllUniqueR.fasta")

#initialize parameter object
sphi_init <- 2;
numMixtures <- 2;
mixDef <- "allUnique";
geneAssignment <- c(rep(1,500), rep(2,500))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)



 phivals <- parameter$readPhiValues( "../ribModel/data/simulatedSelectionSharedR_phi.csv" )
 parameter$initializeSynthesisRateByRandom(phivals)
parameter$initMutationSelectionCategories(c("../ribModel/data/simulated_CSP0.csv", "../ribModel/data/simulated_CSP1.csv") , 2, "Selection")
parameter$initMutationSelectionCategories(c("../ribModel/data/simulated_CSP0.csv", "../ribModel/data/simulated_CSP1.csv") , 2, "Mutation")
# initialize MCMC object
samples <- 200
thining <- 10
adaptiveWidth <- 100
mcmc <- initializeMCMCObject(samples, thining, adaptive.width=adaptiveWidth, 
                     est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
# get model object
model <- initializeModelObject("ROC")

model$setParameter(parameter)
mcmc$run(genome, model)

#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model, parameter)
)

#plots log likelihood trace, possibly other mcmc diagnostics in the future
plot(mcmc)

# plots different aspects of trace
trace <- parameter$getTraceObject()
plot(trace, what = "MixtureProbability")
plot(trace, what = "SPhi")
plot(trace, what = "ExpectedPhi")
plot(trace, what = "Expression", geneIndex = 905)
plot(trace, what = "Mutation", mixture = 1)
plot(trace, what = "Selection", mixture = 1)

# plots model fit (cub plot)
plot(model, genome, parameter, samples = samples*0.1, mixture = 2, main = "S. kluyveri Chr (A,B,Cleft) Codon Usage Plot")

##-----------------------------------------##

#creates a plot Mixture Prob vs log10(phi) (Not part of the package)
num.genes <- genome$getGenomeSize()
mixtureAssignment <- do.call("rbind", lapply(1:num.genes,  function(geneIndex){parameter$getEstimatedMixtureAssignmentProbabilitiesForGene(samples*0.1, geneIndex)}))
plot(NULL, NULL, xlim = c(0,1), ylim = c(-2, 1), xlab = "Probability beeing in Mixture 2", ylab = "log10(phi)")
colors <- c("black", "red")
for(mixture in 1:2)
{
  genes.in.mixture <- which(round(mixtureAssignment[, mixture]+1) == mixture)
  expressionCategory <- parameter$getExpressionCategoryForMixture(mixture)
  
  # need expression values to know range
  num.genes <- length(genes.in.mixture)
  expressionValues <- unlist(lapply(genes.in.mixture, function(geneIndex){
    parameter$getExpressionPosteriorMeanByMixtureElementForGene(samples*0.1, geneIndex, expressionCategory)
  }))
  points(mixtureAssignment[genes.in.mixture, mixture], log10(expressionValues), col=colors[mixture])
}


