
library(testthat)
library(AnaCoDa)
rm(list=ls(all.names=TRUE))
context("MCMC with ROC")

# This file currently checks the logPosterior recorded at iteration 100, between an old, hard-coded test and a current Unit Test.
# Two tests are run: one with Phi, one without Phi. The existence of the relevant input files is also checked.

# Possible implementation change: take the logPosterior value and hard code it here, and compare via
# mcmc$getLogPosteriorTrace(), which returns a vector. Get the average of these values
# and compare it with the hard-coded average of logPosteriorTrace.

# In R, file.path is faster than paste
fileName = file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")
expressionFile = file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR_phi_withPhiSet.csv")
selectionMainFile = file.path("UnitTestingData", "testMCMCROCFiles", "selection_1.csv")
selectionHtFile = file.path("UnitTestingData", "testMCMCROCFiles", "selection_2.csv")
mutationMainFile = file.path("UnitTestingData", "testMCMCROCFiles", "mutation_1.csv")
mutationHtFile = file.path("UnitTestingData", "testMCMCROCFiles", "mutation_2.csv")
mcmcSaveFile <- file.path("UnitTestingOut", "testMCMCROCobject.Rda")

# Ensure the input files exist.
test_that("file exists: simulatedAllUniqueR.fasta", {
  expect_equal(file.exists(fileName), T)
})
test_that("file exists: simulatedAllUniqueR_phi_withPhiSet.csv", {
  expect_equal(file.exists(expressionFile), T)
})
test_that("file exists: selection_1.csv", {
  expect_equal(file.exists(selectionMainFile), T)
})
test_that("file exists: selection_2.csv", {
  expect_equal(file.exists(selectionHtFile), T)
})
test_that("file exists: mutation_1.csv", {
  expect_equal(file.exists(mutationMainFile), T)
})
test_that("file exists: mutation_2.csv", {
  expect_equal(file.exists(mutationHtFile), T)
})

sphi_init <- c(1,1)
numMixtures <- 2
mixDef <- "allUnique"

samples <- 10
thinning <- 10
adaptiveWidth <- 10
divergence.iteration <- 0
knownLogPosteriorValues <- c(with.phi = -933180, without.phi = -949228)
seedValue <- 446141

## Note that length of sample object will be samples + 1
mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, adaptive.width = adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)

### With Phi
set.seed(seedValue)

genome <- initializeGenomeObject(file = fileName, observed.expression.file = expressionFile, match.expression.by.id=FALSE)

geneAssignment <- sample(c(1,2), size = length(genome), replace = TRUE, prob = c(0.3, 0.7)) #c(rep(1,500), rep(2,500))

parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)
parameter$initSelectionCategories(c(selectionMainFile, selectionHtFile), 2,F)
parameter$initMutationCategories(c(mutationMainFile, mutationHtFile), 2,F)

model <- initializeModelObject(parameter, "ROC", with.phi = TRUE) 

outFile <- file.path("UnitTestingOut", "testMCMCROCLogPhi.txt")

sink(outFile)
runMCMC(mcmc = mcmc, genome = genome, model = model, ncores = 1, divergence.iteration = divergence.iteration)
sink()


test_that("identical MCMC-ROC input with Phi, same log posterior", {
   testLogPosterior <- round(mcmc$getLogPosteriorTrace()[(samples + 1)])
   print(testLogPosterior)
   expect_equal(knownLogPosteriorValues[["with.phi"]], testLogPosterior)
})


## 2022-08-15: tests for saving and loading mcmc object added by Elizabeth Barnes and Mike Gilchrist

## Notes:
##   - When using a brand new mcmc object (not one loaded or extended), the first LogPosteriorTrace() and LogLikelihoodTrace() values are both 0.
##   - However, when an mcmc object is loaded, the first index is dropped.
## TEMPORARY RESOLUTION: adjust index value for values
## LONGTERM RESOLUTION: Alex will modify code per issue
##   -`saving and loading MCMC object results in loss of first element #388`


test_that("object can be written successfully: mcmc", {
  expect_null(writeMCMCObject(mcmc = mcmc, file = mcmcSaveFile))
})

test_that("object can be loaded successfully: mcmc", {
  expect_silent(loadMCMCObject(file = mcmcSaveFile))
})

## Loading object in test_that failes to put it in global environment
## Solution: Load explicitly here
mcmcLoaded <- loadMCMCObject(file = mcmcSaveFile)

test_that("object trace matches expected length of (samples): mcmc",{
  expect_equal(
    length(mcmcLoaded$getLogPosteriorTrace()), (samples)) ## note once bug #388 is fixed, replace samples with (samples + 1)
})

test_that("object loaded has expected log posterior", {
   testLogPosterior <- round(mcmc$getLogPosteriorTrace()[(samples)]) ## note once bug #388 is fixed, replace samples with (samples + 1)
   print(testLogPosterior)
   expect_equal(knownLogPosteriorValues[["with.phi"]], testLogPosterior)
})

### end tests by Elizabeth Barnes and Mike Gilchrist


### Without Phi
set.seed(seedValue)

genome <- initializeGenomeObject(file = fileName) 

geneAssignment <- sample(c(1,2), size = length(genome), replace = TRUE, prob = c(0.3, 0.7)) #c(rep(1,500), rep(2,500))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)
parameter$initSelectionCategories(c(selectionMainFile, selectionHtFile), 2,F)
parameter$initMutationCategories(c(mutationMainFile, mutationHtFile), 2,F)

model <- initializeModelObject(parameter, "ROC", with.phi = FALSE) 

outFile = file.path("UnitTestingOut", "testMCMCROCLogWithoutPhi.txt")

sink(outFile)
runMCMC(mcmc, genome, model, 1, divergence.iteration)
sink()
test_that("identical MCMC-ROC input without Phi, same log posterior", {
  testLogPosterior <- round(mcmc$getLogPosteriorTrace()[(samples + 1)])
  print(testLogPosterior)
  expect_equal(knownLogPosteriorValues[["without.phi"]], testLogPosterior)
})


geneAssignment <- sample(c(1,2), size = length(genome), replace = TRUE, prob = c(0.3, 0.7)) #c(rep(1,500), rep(2,500))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)
parameter$initSelectionCategories(c(selectionMainFile, selectionHtFile), 2,F)
parameter$initMutationCategories(c(mutationMainFile, mutationHtFile), 2,F)

model <- initializeModelObject(parameter, "ROC", with.phi = FALSE)

mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, adaptive.width = adaptiveWidth,
                             est.expression=FALSE, est.csp=TRUE, est.hyper=TRUE,est.mix = FALSE)


outFile = file.path("UnitTestingOut", "testMCMCROCLogWithoutPhi.txt")

sink(outFile)
runMCMC(mcmc, genome, model, 1, divergence.iteration)
sink()


numMixtures <- 1
sphi_init <- 1
geneAssignment <- rep(1,length(genome))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)
parameter$initSelectionCategories(c(selectionMainFile), 1,F)
parameter$initMutationCategories(c(mutationMainFile), 1,T)

model <- initializeModelObject(parameter, "ROC", with.phi = FALSE) 

mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, adaptive.width = adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)


outFile = file.path("UnitTestingOut", "testMCMCROCLogWithoutPhi.txt")

sink(outFile)
runMCMC(mcmc, genome, model, 1, divergence.iteration)
sink()

aa <- aminoAcids()
test_that("Making sure DeltaM does not change when fixed", {
  trace <- parameter$getTraceObject()
  for (a in aa)
  {
    if (a == "M" || a == "W" || a == "X") next
    codons <- AAToCodon(a,T)
    for (j in 1:length(codons))
    {
      dm <-  trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,codons[j],0,T)
      expect_equal(var(dm),0)
    }
  }
})

geneAssignment <- rep(1,length(genome))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)
parameter$initSelectionCategories(c(selectionMainFile), 1,T)
parameter$initMutationCategories(c(mutationMainFile), 1,F)

model <- initializeModelObject(parameter, "ROC", with.phi = FALSE) 

mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, adaptive.width = adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)


outFile = file.path("UnitTestingOut", "testMCMCROCLogWithoutPhi.txt")

sink(outFile)
runMCMC(mcmc, genome, model, 1, divergence.iteration)
sink()

## 2022-08-17: tests for saving and loading mcmc object added by Elizabeth Barnes and Mike Gilchrist

test_that("object can be written successfully: mcmc", {
  expect_null(writeMCMCObject(mcmc = mcmc, file = outFile))
})

test_that("object can be loaded successfully: mcmc", {
  expect_silent(mcmcSaved <- loadMCMCObject(file = outFile))  
})

### end tests by Elizabeth Barnes and Mike Gilchrist

aa <- aminoAcids()
test_that("Making sure DeltaEta does not change when fixed", {
  trace <- parameter$getTraceObject()
  for (a in aa)
  {
    if (a == "M" || a == "W" || a == "X") next
    codons <- AAToCodon(a,T)
    for (j in 1:length(codons))
    {
      deta <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1,codons[j],1,T)
      expect_equal(var(deta),0)
    }
  }
})


