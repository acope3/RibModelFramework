library(testthat)
library(ribModel)

context("MCMC with ROC")

# TODO: This file unit checks an entire outputting format, checking not only logLikelihood but also
# other variables (good, if errors occur there) 
# in a set formatted corpus (bad, this means future edits to the printing will result in testthat errors!)
#
# If no tampering is done to the correct variables,
# a correct loglikelihood with the same seed should be equal, disregarding other variables.
#
# Thus, in the future, simply take the logLikelihood value and hard code it here, and compare via
# mcmc$getLogLikelihoodTrace(), which returns a vector. Get the average of these values
# and compare it with the hard-coded average of logLikelihoodTrace.
# This is currently not implemented due to laziness and mild helpfulness, and it is currently working.
# Once it breaks, it should be converted.

# In R, file.path is faster than paste
fileName = file.path("UnitTestingData", "testROCModelFiles", "simulatedAllUniqueR.fasta")
expressionFile = file.path("UnitTestingData", "testROCModelFiles", "simulatedAllUniqueR_phi_withPhiSet.csv")
phiFile = file.path("UnitTestingData", "testROCModelFiles", "simulatedAllUniqueR_phi.csv")

test_that("file exists: simulatedAllUniqueR.fasta", {
  expect_equal(file.exists(fileName), T)
})
test_that("file exists: simulatedAllUniqueR_phi_withPhiSet.csv", {
  expect_equal(file.exists(expressionFile), T)
})
test_that("file exists: simulatedAllUniqueR_phi.csv", {
  expect_equal(file.exists(phiFile), T)
})

set.seed(446141)

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
mcmc <- initializeMCMCObject(samples, thining, adaptive.width=adaptiveWidth, est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE) 

model <- initializeModelObject(parameter, "ROC", with.phi = with.phi) 

file1 = file.path("UnitTestingData", "testMCMC.txt")
file2 = file.path("UnitTestingData", "testMCMCNew.txt")

sink(file2)
runMCMC(mcmc, genome, model, 4, divergence.iteration)
sink()

knownLogLikelihood <- -825482
testLogLikelihood <- mcmc$getLogLikelihoodTrace()[10]
testLogLikelihood
#corpus1 <- paste0(readLines(file1), collapse=" ")
#corpus2 <- paste0(readLines(file2), collapse=" ")

#test_that("identical MCMC-ROC output same seed", {
#  expect_equal(corpus1, corpus2)
#})

test_that("identical MCMC-ROC output same log likelihood", {
  expect_equal(knownLogLikelihood, testLogLikelihood)
})
