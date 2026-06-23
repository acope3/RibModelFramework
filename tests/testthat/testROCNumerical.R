
library(testthat)
library(AnaCoDa)
rm(list = ls(all.names = TRUE))

# ======================================================================
# Tests for ROC CalculateProbabilitiesForCodons (tasks #5 and #9)
#
# Mirrors the structure of testFONSE.R: a hand-rolled oracle encodes
# the same math as ROCModel::calculateCodonProbabilityVector
# (src/ROCModel.cpp lines 852-896, the non-log version used by
# simulateGenome and CalculateProbabilitiesForCodons).
#
# ROC codon probability formula (non-log version):
#
#   If min(dEta_non_ref) >= 0  [unshifted branch]:
#     unnorm[i]   = exp(-dM[i] - phi * dEta[i])  for non-ref i
#     unnorm[ref] = 1
#
#   If min(dEta_non_ref) < 0  [shifted branch]:
#     shift by (dM[minIdx], dEta[minIdx]) to keep exponents bounded
#     unnorm[i]   = exp(-(dM[i]-dM[minIdx]) - phi*(dEta[i]-dEta[minIdx]))
#     unnorm[ref] = exp(dM[minIdx] + phi*dEta[minIdx])
#
#   In both cases probabilities are normalized by the sum.
#
# Note: mutation and selection passed to CalculateProbabilitiesForCodons
# are length (numCodons - 1); the reference codon (alphabetically last)
# is implicit and does NOT appear in the input vectors.
# ======================================================================

context("ROC calculateCodonProbabilityVector (tasks #5, #9)")

# ----------------------------------------------------------------------
# Hand-rolled oracle: mirrors ROCModel::calculateCodonProbabilityVector
# ----------------------------------------------------------------------
roc_prob_oracle <- function(dM, dEta, phi) {
  stopifnot(length(dM) == length(dEta), length(dM) >= 1)
  min_idx <- which.min(dEta)
  min_sel <- dEta[min_idx]
  if (min_sel < 0) {
    unnorm_nonref <- exp(-(dM - dM[min_idx]) - (dEta - min_sel) * phi)
    unnorm_ref    <- exp(dM[min_idx] + min_sel * phi)
    raw <- c(unnorm_nonref, unnorm_ref)
  } else {
    unnorm_nonref <- exp(-dM - dEta * phi)
    raw <- c(unnorm_nonref, 1.0)
  }
  raw / sum(raw)
}

# ----------------------------------------------------------------------
# Minimal ROC model setup (same data file used by testMCMCROC.R)
# ----------------------------------------------------------------------
fastaFile <- file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")
test_that("file exists: simulatedAllUniqueR.fasta", {
  expect_true(file.exists(fastaFile))
})

genome    <- initializeGenomeObject(file = fastaFile)
parameter <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1,
                                        gene.assignment = rep(1, length(genome)))
model     <- initializeModelObject(parameter, "ROC")

# Helper: call RMF (CalculateProbabilitiesForCodons) and compare against oracle
compare_roc_vs_oracle <- function(dM_full, dEta_full, phi, label, tol = 1e-12) {
  n   <- length(dM_full)
  rmf <- model$CalculateProbabilitiesForCodons(dM_full[-n], dEta_full[-n], phi)
  oracle <- roc_prob_oracle(dM_full[-n], dEta_full[-n], phi)
  test_that(paste("ROC matches oracle:", label), {
    expect_equal(length(rmf), n)
    expect_equal(sum(rmf), 1, tolerance = 1e-12)
    expect_equal(rmf, oracle, tolerance = tol)
  })
  invisible(list(rmf = rmf, oracle = oracle))
}

# ----------------------------------------------------------------------
# Task #5 -- correctness cases
# ----------------------------------------------------------------------

# Case 1: 2-codon AA, positive dEta. Unshifted branch.
compare_roc_vs_oracle(
  dM_full   = c(0.20, 0.0),
  dEta_full = c(0.01, 0.0),
  phi       = 1.0,
  label     = "2-codon AA, positive dEta, phi=1 (unshifted)"
)

# Case 2: 2-codon AA, negative dEta. Shifted branch.
compare_roc_vs_oracle(
  dM_full   = c(-0.15, 0.0),
  dEta_full = c(-0.02, 0.0),
  phi       = 1.0,
  label     = "2-codon AA, negative dEta (shifted branch)"
)

# Case 3: 6-codon AA (Leu-like), all positive dEta. Unshifted branch.
compare_roc_vs_oracle(
  dM_full   = c(0.30, -0.20, 0.50, -0.10, 0.40, 0.0),
  dEta_full = c(0.01, 0.03, 0.005, 0.02, 0.015, 0.0),
  phi       = 2.0,
  label     = "6-codon AA, positive dEta, phi=2 (unshifted)"
)

# Case 4: 6-codon AA, mixed-sign dEta. Shifted branch.
compare_roc_vs_oracle(
  dM_full   = c(0.30, -0.20, 0.50, -0.10, 0.40, 0.0),
  dEta_full = c(0.01, -0.02, 0.005, 0.02, -0.015, 0.0),
  phi       = 2.0,
  label     = "6-codon AA, mixed-sign dEta, phi=2 (shifted branch)"
)

# Case 5: dEta=0 reduces to pure mutation multinomial.
test_that("ROC dEta=0 reduces to pure mutation multinomial", {
  dM_full   <- c(0.30, -0.20, 0.50, 0.0)
  dEta_full <- rep(0.0, length(dM_full))
  n         <- length(dM_full)
  rmf       <- model$CalculateProbabilitiesForCodons(dM_full[-n], dEta_full[-n], phi = 1.0)
  raw  <- c(exp(-dM_full[-n]), 1.0)
  pure <- raw / sum(raw)
  expect_equal(rmf, pure, tolerance = 1e-12)
})

# Case 6: dM=0 reduces to pure selection multinomial.
test_that("ROC dM=0 reduces to pure selection multinomial", {
  dM_full   <- rep(0.0, 4)
  dEta_full <- c(0.01, 0.03, 0.005, 0.0)
  n         <- length(dM_full)
  phi       <- 5.0
  rmf       <- model$CalculateProbabilitiesForCodons(dM_full[-n], dEta_full[-n], phi)
  raw  <- c(exp(-dEta_full[-n] * phi), 1.0)
  pure <- raw / sum(raw)
  expect_equal(rmf, pure, tolerance = 1e-12)
})

# ----------------------------------------------------------------------
# Task #9 -- numerical stability and extreme parameter ranges
# ----------------------------------------------------------------------

# Case 7: Very small phi (phi=1e-6) with positive dEta.
# Probabilities approach pure mutation multinomial.
compare_roc_vs_oracle(
  dM_full   = c(0.20, 0.50, 0.0),
  dEta_full = c(0.01, 0.03, 0.0),
  phi       = 1e-6,
  label     = "phi=1e-6 (tiny phi, unshifted) -- stability"
)

# Case 8: Large phi (phi=1e4) with all-positive dEta. Unshifted branch.
# Reference codon (dEta=0) should dominate.
compare_roc_vs_oracle(
  dM_full   = c(0.20, 0.50, 0.0),
  dEta_full = c(0.01, 0.03, 0.0),
  phi       = 1e4,
  label     = "phi=1e4 (large phi, unshifted) -- reference codon dominates"
)

# Case 9: Large phi (phi=1e4) with negative dEta. Shifted branch.
# Most extreme regression test: pre-fix this would produce Inf/NaN via
# the exponent overflow that the shift prevents.
compare_roc_vs_oracle(
  dM_full   = c(-0.30, 0.10, 0.20, 0.0),
  dEta_full = c(-0.02, 0.01, -0.005, 0.0),
  phi       = 1e4,
  label     = "phi=1e4, mixed-sign dEta (shifted branch) -- Inf/NaN regression"
)

# Case 10: Very large positive dM -> one non-ref codon heavily disfavored by mutation.
compare_roc_vs_oracle(
  dM_full   = c(10.0, 0.0, 0.0),
  dEta_full = c(0.01, 0.005, 0.0),
  phi       = 1.0,
  label     = "large dM[1]=10 -- codon 1 probability near zero"
)
test_that("ROC large dM: disfavored codon probability < 1e-4", {
  dM   <- c(10.0, 0.0)
  dEta <- c(0.01, 0.0)
  rmf  <- model$CalculateProbabilitiesForCodons(dM[-length(dM)], dEta[-length(dEta)], phi = 1.0)
  expect_lt(rmf[1], 1e-4)
})

# Case 11: Very large phi with very small positive dEta -> probabilities valid.
test_that("ROC phi=1e6 with small positive dEta: sum=1, all finite, no NaN", {
  dM   <- c(0.10, 0.20, 0.0)
  dEta <- c(0.001, 0.002, 0.0)
  rmf  <- model$CalculateProbabilitiesForCodons(dM[-length(dM)], dEta[-length(dEta)], phi = 1e6)
  expect_equal(sum(rmf), 1, tolerance = 1e-10)
  expect_true(all(is.finite(rmf)))
  expect_true(all(rmf >= 0))
})
