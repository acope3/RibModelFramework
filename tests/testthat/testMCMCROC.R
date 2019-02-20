library(testthat)
library(AnaCoDa)

context("MCMC with ROC")

# This file currently checks the logPosterior recorded at iteration 100, between an old, hard-coded test and a current Unit Test.
# Two tests are run: one with Phi, one without Phi. The existence of the relevant input files is also checked.

# Possible implementation change: take the logPosterior value and hard code it here, and compare via
# mcmc$getLogPosteriorTrace(), which returns a vector. Get the average of these values
# and compare it with the hard-coded average of logPosteriorTrace.

# In R, file.path is faster than paste
fileName = file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")
expressionFile = file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR_phi_withPhiSet.csv")
selectionMainFile = file.path("UnitTestingData", "testMCMCROCFiles", "selection_mod_Ecoli_K12_MG1655_ncbi_main_liberal.csv")
selectionHtFile = file.path("UnitTestingData", "testMCMCROCFiles", "selection_mod_Ecoli_K12_MG1655_ncbi_ht_liberal.csv")
mutationMainFile = file.path("UnitTestingData", "testMCMCROCFiles", "mutation_mod_Ecoli_K12_MG1655_ncbi_main_liberal.csv")
mutationHtFile = file.path("UnitTestingData", "testMCMCROCFiles", "mutation_mod_Ecoli_K12_MG1655_ncbi_ht_liberal.csv")

# Ensure the input files exist.
test_that("file exists: simulatedAllUniqueR.fasta", {
  expect_equal(file.exists(fileName), T)
})
test_that("file exists: simulatedAllUniqueR_phi_withPhiSet.csv", {
  expect_equal(file.exists(expressionFile), T)
})
test_that("file exists: selection_mod_Ecoli_K12_MG1655_ncbi_main_liberal.csv", {
  expect_equal(file.exists(selectionMainFile), T)
})
test_that("file exists: selection_mod_Ecoli_K12_MG1655_ncbi_ht_liberal.csv", {
  expect_equal(file.exists(selectionHtFile), T)
})
test_that("file exists: mutation_mod_Ecoli_K12_MG1655_ncbi_main_liberal.csv", {
  expect_equal(file.exists(mutationMainFile), T)
})
test_that("file exists: mutation_mod_Ecoli_K12_MG1655_ncbi_ht_liberal.csv", {
  expect_equal(file.exists(mutationHtFile), T)
})

sphi_init <- c(1,1)
numMixtures <- 2
mixDef <- "allUnique"

samples <- 10
thinning <- 10
adaptiveWidth <- 10
divergence.iteration <- 0

mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, adaptive.width = adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)

### With Phi
set.seed(446141)

genome <- initializeGenomeObject(file = fileName, observed.expression.file = expressionFile, match.expression.by.id=FALSE)

geneAssignment <- sample(c(1,2), size = length(genome), replace = TRUE, prob = c(0.3, 0.7)) #c(rep(1,500), rep(2,500))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)
parameter$initSelectionCategories(c(selectionMainFile, selectionHtFile), 2)
parameter$initMutationCategories(c(mutationMainFile, mutationHtFile), 2)

model <- initializeModelObject(parameter, "ROC", with.phi = TRUE) 

outFile = file.path("UnitTestingOut", "testMCMCROCLogPhi.txt")

sink(outFile)
runMCMC(mcmc, genome, model, 1, divergence.iteration)
sink()

test_that("identical MCMC-ROC input with Phi, same log posterior", {
  knownLogPosterior <- -953149
  testLogPosterior <- round(mcmc$getLogPosteriorTrace()[10])
  expect_equal(knownLogPosterior, testLogPosterior)
})

### Without Phi
set.seed(446141)

genome <- initializeGenomeObject(file = fileName) 

geneAssignment <- sample(c(1,2), size = length(genome), replace = TRUE, prob = c(0.3, 0.7)) #c(rep(1,500), rep(2,500))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)
parameter$initSelectionCategories(c(selectionMainFile, selectiosHtFile), 2)
parameter$initMutationCategories(c(mutationMainFile, mutationHtFile), 2)

model <- initializeModelObject(parameter, "ROC", with.phi = FALSE) 

outFile = file.path("UnitTestingOut", "testMCMCROCLogWithoutPhi.txt")

sink(outFile)
runMCMC(mcmc, genome, model, 1, divergence.iteration)
sink()

test_that("identical MCMC-ROC input without Phi, same log posterior", {
  knownLogPosterior <- -879815
  testLogPosterior <- round(mcmc$getLogPosteriorTrace()[10])
  expect_equal(knownLogPosterior, testLogPosterior)
})

