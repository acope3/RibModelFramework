rm(list=ls()) 
library(ribModel) 

seeds <- read.table(file = "seed.txt")[,1]
task.id <- 1

set.seed(seeds[task.id]) 
with.phi <- FALSE

if (with.phi) { 
  genome <- initializeGenomeObject(file ="../data/realGenomes/Skluyveri.fasta", expression.file = "../data/simulatedAllUniqueR_phi.csv") 
} else { 
  genome <- initializeGenomeObject(file = "../data/realGenomes/Skluyveri.fasta") 
} 

sphi_init <- c(1,1)
numMixtures <- 2
mixDef <- "selectionShared"
geneAssignment <- c(rep(1,961), rep(2,457), rep(1, 3903))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef) 

samples <- 2000
thining <- 50
adaptiveWidth <- 10
divergence.iteration <- 20
mcmc <- initializeMCMCObject(samples, thining, adaptive.width=adaptiveWidth, est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE) 

setRestartSettings(mcmc, paste(task.id, "_kluyveri_full_selectionShared.rst", sep=""), adaptiveWidth*50, TRUE) 

model <- initializeModelObject(parameter, "ROC", with.phi = with.phi) 

start <- Sys.time() 
system.time(runMCMC(mcmc, genome, model, 4, divergence.iteration)) 
end <- Sys.time() 
end - start 


writeParameterObject(parameter, file=paste(task.id, "_ROCParameter.Rdat", sep=""))
writeMCMCObject(mcmc, file=paste(task.id, "_MCMCObject.Rdat", sep=""))

