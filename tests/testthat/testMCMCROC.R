library(testthat)
library(ribModel)

context("MCMC with ROC")

# This file currently checks the logLikelihood recorded at iteration 100, between an old, hard-coded test and a current Unit Test.

# Possible implementation change: take the logLikelihood value and hard code it here, and compare via
# mcmc$getLogLikelihoodTrace(), which returns a vector. Get the average of these values
# and compare it with the hard-coded average of logLikelihoodTrace.


# In R, file.path is faster than paste
fileName = file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")
expressionFile = file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR_phi_withPhiSet.csv")

# Ensure the input files exist.
test_that("file exists: simulatedAllUniqueR.fasta", {
  expect_equal(file.exists(fileName), T)
})
test_that("file exists: simulatedAllUniqueR_phi_withPhiSet.csv", {
  expect_equal(file.exists(expressionFile), T)
})

set.seed(446141)

# TODO: Try this with FALSE as well.
with.phi <- TRUE # FALSE

if (with.phi) {
  genome <- initializeGenomeObject(file = fileName, expression.file = expressionFile) 
} else {
  genome <- initializeGenomeObject(file = fileName) 
}

sphi_init <- c(1,1)
numMixtures <- 2
mixDef <- "allUnique"
geneAssignment <- sample(c(1,2), size = length(genome), replace = TRUE, prob = c(0.3, 0.7)) #c(rep(1,500), rep(2,500))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)

samples <- 10
thining <- 10
adaptiveWidth <- 10
divergence.iteration <- 0
mcmc <- initializeMCMCObject(samples = samples, thining = thining, adaptive.width = adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE) 

model <- initializeModelObject(parameter, "ROC", with.phi = with.phi) 

#outFile = file.path("UnitTestingOut", "testMCMCROCLog.txt")

#sink(outFile)
runMCMC(mcmc, genome, model, 1, divergence.iteration)
#sink()

test_that("identical MCMC-ROC input, same log likelihood", {
  knownLogLikelihood <- -825482
  testLogLikelihood <- round(mcmc$getLogLikelihoodTrace()[10])
  expect_equal(knownLogLikelihood, testLogLikelihood)
})
