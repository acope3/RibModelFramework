library(ribModel)
rm(list=ls())
genome <- initializeGenomeObject(file = "../ribModel/data/simulatedAllUniqueR.fasta")

#initialize parameter object
sphi_init <- 1
numMixtures <- 2
mixDef <- "allUnique"
geneAssignment <- c(rep(1,500), rep(2,500))

parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)

# initialize MCMC object
samples <- 1000
thining <- 10
adaptiveWidth <- 10
divergence.iteration <- 50
mcmc <- initializeMCMCObject(samples, thining, adaptive.width=adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)

# get model object
model <- initializeModelObject(parameter, "ROC", with.phi = FALSE)

#run mcmc on genome with parameter using model
n.cores <- 4
start <- Sys.time()
system.time(
  runMCMC(mcmc, genome, model, n.cores, divergence.iteration)
)
end <- Sys.time()
end - start