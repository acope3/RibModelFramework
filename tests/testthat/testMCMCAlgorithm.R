library(testthat)
library(ribModel)

context("MCMCAlgorithm")

test_that("general MCMCAlgorithm functions", {
  expect_equal(testMCMCAlgorithm(), 0)
})

samples <- 1000
thining <- 1
adaptiveWidth <- 100
mcmc <- new(MCMCAlgorithm, samples, thining, adaptiveWidth, TRUE, TRUE, TRUE)

test_that("get Samples", {
  expect_equal(mcmc$getSamples(), samples)
})

test_that("get Thining", {
  expect_equal(mcmc$getThining(), thining)
})

test_that("get Adaptive Width", {
  expect_equal(mcmc$getAdaptiveWidth(), adaptiveWidth)
})

test_that("set Samples", {
  mcmc$setSamples(10)
  expect_equal(mcmc$getSamples(), 10)
})

test_that("set Thining", {
  mcmc$setThining(10)
  expect_equal(mcmc$getThining(), 10)
})

test_that("set Adaptive Width", {
  mcmc$setAdaptiveWidth(10)
  expect_equal(mcmc$getAdaptiveWidth(), 10)
})

test_that("set Log Likelihood Trace", {
  vect <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  mcmc$setLogLikelihoodTrace(vect)
  expect_equal(mcmc$getLogLikelihoodTrace(), vect)
})