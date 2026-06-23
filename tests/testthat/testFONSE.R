
library(testthat)
library(AnaCoDa)
rm(list=ls(all.names=TRUE))
context("FONSE numerical correctness")

# This file validates the FONSE codon-probability computation at the
# numerical level, independently of MCMC. For a set of hand-chosen
# (mutation, selection, phi, a1, position) tuples, it compares the output
# of FONSEModel$CalculateProbabilitiesForCodons() against a hand-rolled
# oracle that encodes the same math as src/FONSEModel.cpp lines 961-1019
# (the non-log version used by simulateGenome).
#
# FONSE codon probability (non-log version, matches the C++ simulator path):
#
#   beta(pos) = a1 + a2 * pos           (a2 is hardcoded to 4.0 in C++)
#
#   If min(dOmega_non_ref) >= 0  [unshifted branch, lines 999-1009]:
#     unnorm[i]    = exp(-dM[i] - phi * beta(pos) * dOmega[i])   for non-ref i
#     unnorm[ref]  = 1
#
#   If min(dOmega_non_ref) < 0   [shifted branch, lines 988-998]:
#     shift by (dM[minIdx], dOmega[minIdx]) to keep exponents bounded
#     unnorm[i]    = exp(-(dM[i]-dM[minIdx]) - phi*beta(pos)*(dOmega[i]-dOmega[minIdx]))
#     unnorm[ref]  = exp( dM[minIdx] + phi*beta(pos)*dOmega[minIdx])
#
# In both cases, probabilities are normalized by the sum.
#
# Note: `mutation` and `selection` passed to CalculateProbabilitiesForCodons
# are length (numCodons - 1); the reference codon (alphabetically last) is
# implicit and does NOT appear in the input vectors. Output has length
# numCodons with the reference probability at the end.
#
# Note: a1 = 4 ATP is the biological initiation cost used here (Gilchrist 2007
# and successors). Note a2 is ALSO hardcoded to 4.0 in the C++; this
# numerical coincidence is not a biological identity and is flagged as a
# separate TODO (promote a2 to free parameter).

A2_HARDCODED <- 4.0  # must match FONSEModel.cpp lines 992, 996, 1004

# ----------------------------------------------------------------------
# Hand-rolled oracle: mirror FONSEModel::calculateCodonProbabilityVector
# (non-log version) exactly.
# ----------------------------------------------------------------------
fonse_prob_oracle <- function(dM, dOmega, phi, a1, position, a2 = A2_HARDCODED) {
  stopifnot(length(dM) == length(dOmega), length(dM) >= 1)
  beta_pos <- a1 + a2 * position
  min_idx  <- which.min(dOmega)
  min_sel  <- dOmega[min_idx]
  if (min_sel < 0) {
    # shifted branch
    unnorm_nonref <- exp(-(dM - dM[min_idx]) -
                         phi * beta_pos * (dOmega - min_sel))
    unnorm_ref    <- exp(dM[min_idx] + phi * beta_pos * min_sel)
    raw <- c(unnorm_nonref, unnorm_ref)
  } else {
    # unshifted branch
    unnorm_nonref <- exp(-dM - phi * beta_pos * dOmega)
    raw <- c(unnorm_nonref, 1.0)
  }
  raw / sum(raw)
}

# ----------------------------------------------------------------------
# Minimal FONSE model setup. CalculateProbabilitiesForCodons only uses
# its arguments; the model object is just the namespace the method lives
# in, so we construct the cheapest valid model we can.
# ----------------------------------------------------------------------
fastaFile <- file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")
test_that("file exists: simulatedAllUniqueR.fasta", {
  expect_true(file.exists(fastaFile))
})

a1_literature <- 4.0  # ATP cost of translation initiation, Gilchrist 2007
genome    <- initializeGenomeObject(file = fastaFile)
parameter <- initializeParameterObject(genome = genome,
                                       sphi = 1,
                                       num.mixtures = 1,
                                       gene.assignment = rep(1, length(genome)),
                                       model = "FONSE",
                                       init.initiation.cost = a1_literature)
model     <- initializeModelObject(parameter, "FONSE")

# Helper: given full-length dM/dOmega vectors (including reference at the
# end), call RMF with the non-reference portion and return the full-length
# RMF output. Compare against the oracle.
compare_rmf_vs_oracle <- function(dM_full, dOmega_full, phi, a1, position,
                                  label, tol = 1e-12) {
  n <- length(dM_full)
  stopifnot(length(dOmega_full) == n, n >= 2)
  # RMF takes non-reference vectors (length n-1). Reference is the last.
  rmf <- model$CalculateProbabilitiesForCodons(
    dM_full[-n], dOmega_full[-n], phi, a1, position
  )
  oracle <- fonse_prob_oracle(dM_full[-n], dOmega_full[-n], phi, a1, position)
  test_that(paste("FONSE matches oracle:", label), {
    expect_equal(length(rmf), n)
    expect_equal(sum(rmf), 1, tolerance = 1e-12)
    expect_equal(rmf, oracle, tolerance = tol)
  })
  invisible(list(rmf = rmf, oracle = oracle))
}

# ----------------------------------------------------------------------
# Test cases. dOmega values are in the physically reasonable range
# for FONSE (odds ratio p/(1-p) ~ b/c, where b~0.00515/s and c~10/s,
# giving dOmega ~ 5e-4). Some tests push outside this range to cover
# edge cases (mixed-sign shifted branch, large phi).
# ----------------------------------------------------------------------

# Case 1: 2-codon AA, all-positive dOmega, small phi. Unshifted branch.
#   Biologically plausible: Phe (TTT, TTC), dOmega ~ 5e-4.
compare_rmf_vs_oracle(
  dM_full     = c(0.20, 0.0),
  dOmega_full = c(5e-4, 0.0),
  phi         = 1.0,
  a1          = a1_literature,
  position    = 100,
  label       = "2-codon AA, positive dOmega, small phi"
)

# Case 2: 2-codon AA, mixed-sign dOmega forcing the shifted branch.
#   dOmega[1] = -1e-3, dOmega[ref] = 0 (implicit). min < 0 triggers shift.
compare_rmf_vs_oracle(
  dM_full     = c(-0.10, 0.0),
  dOmega_full = c(-1e-3, 0.0),
  phi         = 1.0,
  a1          = a1_literature,
  position    = 100,
  label       = "2-codon AA, negative dOmega (shifted branch)"
)

# Case 3: 6-codon AA (Leu-like), all positive dOmega. Unshifted branch.
compare_rmf_vs_oracle(
  dM_full     = c(0.30, -0.20, 0.50, -0.10, 0.40, 0.0),
  dOmega_full = c(4e-4, 6e-4, 3e-4, 8e-4, 5e-4, 0.0),
  phi         = 1.0,
  a1          = a1_literature,
  position    = 100,
  label       = "6-codon AA, positive dOmega"
)

# Case 4: Large phi (1e4) with positive dOmega. Tests numerical stability
# of the unshifted branch at the phi range where overflow becomes
# plausible. phi * beta(pos) * dOmega = 1e4 * (4 + 4*100) * 5e-4 ~ 2020,
# i.e. exp(-2020) ~ 0 and normalization stays well-defined.
compare_rmf_vs_oracle(
  dM_full     = c(0.20, 0.50, 0.0),
  dOmega_full = c(5e-4, 3e-4, 0.0),
  phi         = 1e4,
  a1          = a1_literature,
  position    = 100,
  label       = "Large phi (1e4), positive dOmega"
)

# Case 5: Large phi with mixed-sign dOmega forcing the shifted branch.
# This is the regression test for bugfix 4af0635 (min/max inversion in
# calculateCodonProbabilityVector). Pre-fix, this case would produce
# Inf/NaN because the shift failed to fire.
compare_rmf_vs_oracle(
  dM_full     = c(-0.30, 0.10, 0.20, 0.0),
  dOmega_full = c(-8e-4, 6e-4, -3e-4, 0.0),
  phi         = 1e4,
  a1          = a1_literature,
  position    = 500,
  label       = "Large phi (1e4), mixed-sign dOmega (shifted branch, bugfix regression)"
)

# Case 6: All-zero dOmega reduces to pure mutation multinomial.
# P(i) = exp(-dM[i]) / (sum_k exp(-dM[k]) + 1)
# P(ref) = 1 / (sum_k exp(-dM[k]) + 1)
test_that("FONSE dOmega=0 reduces to pure mutation multinomial", {
  dM_full     <- c(0.30, -0.20, 0.50, 0.0)
  dOmega_full <- rep(0.0, length(dM_full))
  n <- length(dM_full)
  rmf <- model$CalculateProbabilitiesForCodons(
    dM_full[-n], dOmega_full[-n], phi = 1.0, a1_literature, position = 100
  )
  # Pure mutation reference oracle, independent of phi/beta/position.
  raw <- c(exp(-dM_full[-n]), 1.0)
  pure <- raw / sum(raw)
  expect_equal(rmf, pure, tolerance = 1e-12)
})

# ======================================================================
# FONSE MCMC Integration Tests
#
# Run a short chain and verify structural correctness (trace dimensions,
# finiteness, monotone burn-in). We do NOT hardcode a specific logPosterior
# value because FONSE's position-dependent likelihood starts very far from
# the posterior mode; 10 samples is still deep in burn-in, and run-to-run
# variance from RNG-inside-OpenMP (see testMCMCROC.R note) is large relative
# to the initial transient.
# ======================================================================
context("FONSE MCMC integration")

samples       <- 10
thinning      <- 10
adaptiveWidth <- 10

set.seed(446141)
# Reuse genome/model objects constructed above for the numerical tests.
mcmc <- initializeMCMCObject(samples = samples, thinning = thinning,
                              adaptive.width = adaptiveWidth,
                              est.expression = TRUE, est.csp = TRUE,
                              est.hyper = TRUE)

outFile <- file.path("UnitTestingOut", "testFONSEMCMCLog.txt")
sink(outFile)
runMCMC(mcmc, genome, model, 1, 0)
sink()

test_that("FONSE MCMC logPosterior trace has length samples+1", {
  expect_equal(length(mcmc$getLogPosteriorTrace()), samples + 1)
})

test_that("FONSE MCMC logPosterior is finite at all non-initial samples", {
  lp <- mcmc$getLogPosteriorTrace()
  expect_true(all(is.finite(lp[-1])))  # skip index 1 (initial 0)
})

test_that("FONSE MCMC logPosterior increases during burn-in (first < last)", {
  lp <- mcmc$getLogPosteriorTrace()
  # During burn-in, posterior increases (less negative) toward mode.
  expect_true(lp[2] < lp[samples + 1])
})

test_that("FONSE MCMC initiation cost (a1) trace has length samples+1", {
  trace     <- parameter$getTraceObject()
  initTrace <- trace$getInitiationCostTrace()
  expect_equal(length(initTrace), samples + 1)
})

test_that("FONSE MCMC initiation cost (a1) is positive and finite at last sample", {
  trace     <- parameter$getTraceObject()
  initTrace <- trace$getInitiationCostTrace()
  expect_true(all(is.finite(initTrace[-1])))
  expect_gt(initTrace[samples + 1], 0)
})

test_that("FONSE MCMC codon-specific parameter (dM, dOmega) traces have length samples+1", {
  trace    <- parameter$getTraceObject()
  muTrace  <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1, "AAA", 0, FALSE)
  selTrace <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1, "AAA", 1, FALSE)
  expect_equal(length(muTrace),  samples + 1)
  expect_equal(length(selTrace), samples + 1)
})

test_that("FONSE MCMC synthesis rate trace has length samples+1", {
  trace      <- parameter$getTraceObject()
  synthTrace <- trace$getSynthesisRateTraceForGene(1)
  expect_equal(length(synthTrace), samples + 1)
})
