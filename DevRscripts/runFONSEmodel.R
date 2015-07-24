library(ribModel)
rm(list=ls())
#read genome
genome <- initializeGenomeObject(file = "../ribModel/data/Skluyveri_ChrA_andCleft.fasta")

#initialize parameter object
sphi_init <- 2
numMixtures <- 2
mixDef <- "allUnique"
#geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457), rep(1, 3903))
#geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457))
geneAssignment <- c(rep(1,448), rep(2,457))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, model= "FONSE", split.serine = TRUE,
                                       mixture.definition = mixDef)

# initialize MCMC object
samples <- 100
thining <- 10
adaptiveWidth <- 10
mcmc <- initializeMCMCObject(samples=samples, thining=thining, adaptive.width=adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
# get model object
model <- initializeModelObject(parameter, "FONSE")

setRestartSettings(mcmc, "restartFile.rst", adaptiveWidth, TRUE)
#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model)
)