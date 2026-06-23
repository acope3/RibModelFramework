library(AnaCoDa)
library(testthat)

local_edition(3)

fileName         <- file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")
expressionFile   <- file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR_phi_withPhiSet.csv")
selectionMainFile <- file.path("UnitTestingData", "testMCMCROCFiles", "selection_1.csv")
selectionHtFile  <- file.path("UnitTestingData", "testMCMCROCFiles", "selection_2.csv")
mutationMainFile <- file.path("UnitTestingData", "testMCMCROCFiles", "mutation_1.csv")
mutationHtFile   <- file.path("UnitTestingData", "testMCMCROCFiles", "mutation_2.csv")

sphi_init          <- c(1, 1)
numMixtures        <- 2
mixDef             <- "allUnique"
samples            <- 10
thinning           <- 10
adaptiveWidth      <- 10
divergence.iteration <- 0

set.seed(446141)
genome <- initializeGenomeObject(file = fileName,
                                  observed.expression.file = expressionFile,
                                  match.expression.by.id = FALSE)
geneAssignment <- sample(c(1, 2), size = length(genome), replace = TRUE,
                          prob = c(0.3, 0.7))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment,
                                        split.serine = TRUE, mixture.definition = mixDef)
parameter$initSelectionCategories(c(selectionMainFile, selectionHtFile), 2, F)
parameter$initMutationCategories(c(mutationMainFile, mutationHtFile), 2, F)
model <- initializeModelObject(parameter, "ROC", with.phi = TRUE)

mcmc <- initializeMCMCObject(samples = samples, thinning = thinning,
                              adaptive.width = adaptiveWidth,
                              est.expression = TRUE, est.csp = TRUE, est.hyper = TRUE)

outFile <- file.path("UnitTestingOut", "testMCMCROCLogPhi.txt")
sink(outFile)
runMCMC(mcmc, genome, model, 1, divergence.iteration)
sink()

test_that("ROC MCMC logLikelihood trace final value (all estimation, seed 446141)", {
  expect_snapshot(mcmc$getLogLikelihoodTrace()[samples])
})
