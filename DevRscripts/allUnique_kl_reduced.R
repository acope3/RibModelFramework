rm(list=ls())
library(ribModel)
seeds <- read.table(file = "seed.txt")[,1]
task.id <- as.numeric(Sys.getenv("SGE_TASK_ID"))

set.seed(seeds[task.id])
#set.seed(446141)
with.phi <- FALSE

if (with.phi) {
  genome <- initializeGenomeObject(file = "../data/realGenomes/Skluyveri_reduced.fasta", expression.file = "../data/simulatedAllUniqueR_phi.csv") 
} else {
  genome <- initializeGenomeObject(file = "../data/realGenomes/Skluyveri_reduced.fasta") 
}

sphi_init <- c(1,1)
numMixtures <- 2
mixDef <- "allUnique"
geneAssignment <- c(rep(1,399), rep(2,457), rep(1, 644))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)

samples <- 1000
thining <- 50
adaptiveWidth <- 10
divergence.iteration <- 0
mcmc <- initializeMCMCObject(samples, thining, adaptive.width=adaptiveWidth, est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE, est.mix = TRUE) 

setRestartSettings(mcmc, paste("_", task.id, "_kluyveri_reduced_allUnique.rst", sep=""), adaptiveWidth*50, TRUE) 

model <- initializeModelObject(parameter, "ROC", with.phi = with.phi) 

start <- Sys.time() 
system.time(runMCMC(mcmc, genome, model, 8, divergence.iteration)) 
end <- Sys.time() 
end - start 


writeParameterObject(parameter, file=paste(task.id, "_kl_reduced_allUnique_ROCParameter.Rda", sep=""))
writeMCMCObject(mcmc, file=paste(task.id, "_kl_reduced_allUnique_MCMCObject.Rda", sep=""))
