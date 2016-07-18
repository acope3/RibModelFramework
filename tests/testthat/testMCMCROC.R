library(testthat)
library(ribModel)

context("MCMC with ROC")

# This file currently checks the logLikelihood recorded at iteration 100, between an old, hard-coded test and a current Unit Test.
# Two tests are run: one with Phi, one without Phi. The existence of the relevant input files is also checked.

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

sphi_init <- c(1,1)
numMixtures <- 2
mixDef <- "allUnique"

samples <- 10
thining <- 10
adaptiveWidth <- 10
divergence.iteration <- 0

mcmc <- initializeMCMCObject(samples = samples, thining = thining, adaptive.width = adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)

### With Phi
set.seed(446141)

genome <- initializeGenomeObject(file = fileName, expression.file = expressionFile)

geneAssignment <- sample(c(1,2), size = length(genome), replace = TRUE, prob = c(0.3, 0.7)) #c(rep(1,500), rep(2,500))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)

model <- initializeModelObject(parameter, "ROC", with.phi = TRUE) 

outFile = file.path("UnitTestingOut", "testMCMCROCLogPhi.txt")

sink(outFile)
runMCMC(mcmc, genome, model, 1, divergence.iteration)
sink()

test_that("identical MCMC-ROC input with Phi, same log likelihood", {
  knownLogLikelihood <- -825482
  testLogLikelihood <- round(mcmc$getLogLikelihoodTrace()[10])
  expect_equal(knownLogLikelihood, testLogLikelihood)
})

### Without Phi
set.seed(446141)

genome <- initializeGenomeObject(file = fileName) 

geneAssignment <- sample(c(1,2), size = length(genome), replace = TRUE, prob = c(0.3, 0.7)) #c(rep(1,500), rep(2,500))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)

model <- initializeModelObject(parameter, "ROC", with.phi = FALSE) 

outFile = file.path("UnitTestingOut", "testMCMCROCLogWithoutPhi.txt")

sink(outFile)
runMCMC(mcmc, genome, model, 1, divergence.iteration)
sink()

test_that("identical MCMC-ROC input without Phi, same log likelihood", {
  knownLogLikelihood <- -825693
  testLogLikelihood <- round(mcmc$getLogLikelihoodTrace()[10])
  expect_equal(knownLogLikelihood, testLogLikelihood)
})

