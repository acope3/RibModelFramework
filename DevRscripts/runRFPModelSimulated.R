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


#init from "true" values

phiVals <- parameter$readPhiValues( "../ribModel/data/RFPPhiValues.csv" )
parameter$initializeSynthesisRateByRandom(phiVals)
parameter$initMutationSelectionCategories(c("../ribModel/data/RFPAlphaValues.csv"), 1, "Alpha")
parameter$initMutationSelectionCategories(c("../ribModel/data/RFPLambdaPrimeValues.csv"), 1, "LambdaPrime")


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

