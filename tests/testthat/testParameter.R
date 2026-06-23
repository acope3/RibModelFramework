library(testthat)
library(AnaCoDa)

context("Parameter")

test_that("general parameter functions", {
  expect_equal(testParameter("UnitTestingData/testMCMCROCFiles"), 0)
})

# ======================================================================
# ROC parameter prior getter/setter and dimension tests -- task #7
# ======================================================================
context("ROC parameter priors and proposal dimensions (task #7)")

fastaFile <- file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")
genome    <- initializeGenomeObject(file = fastaFile)
n_genes   <- length(genome)

parameter <- initializeParameterObject(
  genome         = genome,
  sphi           = 1,
  num.mixtures   = 1,
  gene.assignment = rep(1, n_genes)
)

# getMutationPriorMean/SD return a list-of-lists (categories x params).
# setMutationPriorMean/SD take a flat numeric vector (row-major).

# --- prior returns are lists with finite positive-SD values ---
test_that("ROC mutation prior mean is a list and all elements are finite", {
  mu_mean <- parameter$getMutationPriorMean()
  expect_true(is.list(mu_mean))
  expect_true(all(is.finite(unlist(mu_mean))))
})

test_that("ROC mutation prior SD is a list with all-positive values (default)", {
  mu_sd <- parameter$getMutationPriorStandardDeviation()
  expect_true(is.list(mu_sd))
  expect_true(all(unlist(mu_sd) > 0))
})

test_that("ROC selection prior mean is a list and all elements are finite", {
  sel_mean <- parameter$getSelectionPriorMean()
  expect_true(is.list(sel_mean))
  expect_true(all(is.finite(unlist(sel_mean))))
})

test_that("ROC selection prior SD is a list with all-positive values (default)", {
  sel_sd <- parameter$getSelectionPriorStandardDeviation()
  expect_true(is.list(sel_sd))
  expect_true(all(unlist(sel_sd) > 0))
})

# --- prior setter round-trip (flat-vector interface) ---
test_that("ROC mutation prior mean round-trips via set/get", {
  original   <- parameter$getMutationPriorMean()
  flat_orig  <- unlist(original)
  flat_new   <- flat_orig + 0.5
  parameter$setMutationPriorMean(flat_new)
  retrieved  <- unlist(parameter$getMutationPriorMean())
  expect_equal(retrieved, flat_new, tolerance = 1e-10)
  parameter$setMutationPriorMean(flat_orig)  # restore
})

test_that("ROC mutation prior SD round-trips via set/get", {
  original  <- parameter$getMutationPriorStandardDeviation()
  flat_orig <- unlist(original)
  flat_new  <- flat_orig * 2.0
  parameter$setMutationPriorStandardDeviation(flat_new)
  retrieved <- unlist(parameter$getMutationPriorStandardDeviation())
  expect_equal(retrieved, flat_new, tolerance = 1e-10)
  parameter$setMutationPriorStandardDeviation(flat_orig)  # restore
})

test_that("ROC selection prior SD round-trips via set/get", {
  original  <- parameter$getSelectionPriorStandardDeviation()
  flat_orig <- unlist(original)
  flat_new  <- flat_orig * 0.5
  parameter$setSelectionPriorStandardDeviation(flat_new)
  retrieved <- unlist(parameter$getSelectionPriorStandardDeviation())
  expect_equal(retrieved, flat_new, tolerance = 1e-10)
  parameter$setSelectionPriorStandardDeviation(flat_orig)  # restore
})

# --- current parameter dimensions ---
test_that("ROC currentMutationParameter is a list (one element per mixture)", {
  mu <- parameter$currentMutationParameter
  expect_true(is.list(mu))
  expect_equal(length(mu), 1)          # 1 mixture
  expect_true(is.numeric(mu[[1]]))
  expect_gt(length(mu[[1]]), 0)
})

test_that("ROC currentSelectionParameter is a list (one element per mixture)", {
  sel <- parameter$currentSelectionParameter
  expect_true(is.list(sel))
  expect_equal(length(sel), 1)
  expect_true(is.numeric(sel[[1]]))
  expect_gt(length(sel[[1]]), 0)
})

test_that("ROC currentMutation and currentSelection vectors have the same length", {
  mu  <- parameter$currentMutationParameter[[1]]
  sel <- parameter$currentSelectionParameter[[1]]
  expect_equal(length(mu), length(sel))
})

# --- proposed parameter structure ---
test_that("ROC proposedMutationParameter is a list with numeric elements", {
  prop <- parameter$proposedMutationParameter
  expect_true(is.list(prop))
  expect_equal(length(prop), 1)
  expect_true(is.numeric(prop[[1]]))
})

test_that("ROC proposedSelectionParameter is a list with numeric elements", {
  prop <- parameter$proposedSelectionParameter
  expect_true(is.list(prop))
  expect_equal(length(prop), 1)
  expect_true(is.numeric(prop[[1]]))
})
