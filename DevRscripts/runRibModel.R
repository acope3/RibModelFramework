library(ribModel)

#read genome
genome <- initializeGenomeObject(fasta.file = "../ribModel/data/Skluyveri_ChrA_ChrB_andCleft.fasta")

#initialize parameter object
sphi_init <- 2;
numMixtures <- 2;
mixDef <- "allUnique";
geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)

# initialize MCMC object
samples <- 100
thining <- 10
adaptiveWidth <- 10
mcmc <- initializeMCMCObject(samples, thining, adaptive.width=adaptiveWidth, 
                     est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
# get model object
model <- initializeModelObject("ROC")

#run mcmc on genome with parameter using model
runMCMC(mcmc, genome, model, parameter);

#plots log likelihood trace, possibly other mcmc diagnostics in the future
plot(mcmc)
# plots different aspects of parameter
plot(parameter, what = "MixtureProbability")
plot(parameter, what = "SPhi")
plot(parameter, what = "ExpectedPhi")
plot(parameter, what = "Expression", geneIndex = 905)
plot(parameter, what = "Mutation", mixture = 1)
plot(parameter, what = "Selection", mixture = 1)

# plots model fit (cub plot)
plot(model, genome, parameter, samples = samples*0.1, mixture = 1, main = "S. kluyveri Chr (A,B,Cleft) Codon Usage Plot")

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
    parameter$getExpressionPosteriorMeanByExpressionCategoryForGene(samples*0.1, geneIndex, expressionCategory)
  }))
  points(mixtureAssignment[genes.in.mixture, mixture], log10(expressionValues), col=colors[mixture])
}


