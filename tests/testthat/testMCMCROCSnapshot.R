library(AnaCoDa)

library(testthat)
local_edition(3)
Config/testthat/edition: 3

rm(list=ls(all.names=TRUE))
context("MCMC with ROC")

fileName = file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")
expressionFile = file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR_phi_withPhiSet.csv")
selectionMainFile = file.path("UnitTestingData", "testMCMCROCFiles", "selection_1.csv")
selectionHtFile = file.path("UnitTestingData", "testMCMCROCFiles", "selection_2.csv")
mutationMainFile = file.path("UnitTestingData", "testMCMCROCFiles", "mutation_1.csv")
mutationHtFile = file.path("UnitTestingData", "testMCMCROCFiles", "mutation_2.csv")
mcmcSaveFile <- file.path("UnitTestingOut", "testMCMCROCobject.Rda")

sphi_init <- c(1,1)
numMixtures <- 2
mixDef <- "allUnique"

samples <- 10
thinning <- 10
adaptiveWidth <- 10
divergence.iteration <- 0

#run MCMC - with phi
mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, adaptive.width = adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)

set.seed(446141)

genome <- initializeGenomeObject(file = fileName, observed.expression.file = expressionFile, match.expression.by.id=FALSE)

geneAssignment <- sample(c(1,2), size = length(genome), replace = TRUE, prob = c(0.3, 0.7)) #c(rep(1,500), rep(2,500))

parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)
parameter$initSelectionCategories(c(selectionMainFile, selectionHtFile), 2,F)
parameter$initMutationCategories(c(mutationMainFile, mutationHtFile), 2,F)

model <- initializeModelObject(parameter, "ROC", with.phi = TRUE) 

outFile <- file.path("UnitTestingOut", "testMCMCROCLogPhi.txt")

sink(outFile)
runMCMC(mcmc, genome, model, 1, divergence.iteration)
sink()

print(mcmc$getLogLikelihoodTrace()[10])
#known logLik for est_phi <- -2250948

#2022-09-09: begin snapshot testing by Elizabeth Barnes
T_or_F <- data.frame(est_phi = c("T", "T", "T", "F", "F", "F"),
                 est_dM = c("T", "F", "T", "T", "F", "F"),
                 est_dE = c("T", "F", "F", "T", "T","F"),
                 logLik = c(-2250948, -2250948, -2250948, "NULL", "NULL", "NULL"))


test_cond <- c("est_phi", "est_dM", "est_dE")
test_val <- c("logLik", "vardM", "vardE")
test_msg <- paste0("Test cond: est_phi", "Test value: logLik", print(mcmc$getLogLikelihoodTrace()[10]))

for(est_phi in T_or_F) {
  test_cond <- T_or_F$est_phi=="T"

test_that(test_msg, {
  local_edition(3)
  expect_snapshot(T_or_F$est_phi[4], mcmc$getLogLikelihoodTrace()[10])
})
      

}


             
