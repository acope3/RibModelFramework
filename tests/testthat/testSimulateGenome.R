
library(testthat)
library(AnaCoDa)
rm(list = ls(all.names = TRUE))

# ======================================================================
# Tests for simulateGenome() -- task #6
#
# simulateGenome(genome) modifies the genome in-place: for ROC/FONSE it
# re-draws codon sequences from the current parameter distribution; for
# PA/PANSE it re-draws RFP counts. The method returns void (invisible NULL
# in R). Tests verify:
#   - The call completes without error
#   - The number of genes is unchanged
#   - Simulated sequences are valid DNA and have the same length
#   - Simulated sequences differ from the originals (probabilistic; may
#     occasionally fail for a tiny genome; seed is set for reproducibility)
# ======================================================================

context("simulateGenome -- ROC and FONSE (task #6)")

fastaFile <- file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")

# ----------------------------------------------------------------------
# Helper: given a genome, return a character vector of all gene sequences
# ----------------------------------------------------------------------
genome_seqs <- function(genome) {
  genes <- genome$getGenes(simulated = FALSE)
  vapply(genes, function(g) g$seq, character(1))
}

sim_seqs <- function(genome) {
  genes <- genome$getGenes(simulated = TRUE)
  vapply(genes, function(g) g$seq, character(1))
}

# ----------------------------------------------------------------------
# ROC simulateGenome
# ----------------------------------------------------------------------
set.seed(42)
genome_roc    <- initializeGenomeObject(file = fastaFile)
original_seqs <- genome_seqs(genome_roc)
n_genes       <- length(genome_roc)
original_lens <- nchar(original_seqs)

parameter_roc <- initializeParameterObject(
  genome         = genome_roc,
  sphi           = 1,
  num.mixtures   = 1,
  gene.assignment = rep(1, n_genes)
)
model_roc <- initializeModelObject(parameter_roc, "ROC")

test_that("ROC simulateGenome returns NULL invisibly", {
  result <- model_roc$simulateGenome(genome_roc)
  expect_null(result)
})

test_that("ROC simulateGenome preserves gene count", {
  expect_equal(length(genome_roc), n_genes)
})

test_that("ROC simulated sequences have the same length as originals", {
  sim <- sim_seqs(genome_roc)
  expect_equal(nchar(sim), original_lens)
})

test_that("ROC simulated sequences contain only valid DNA bases", {
  sim <- sim_seqs(genome_roc)
  invalid <- grep("[^ACGTacgt]", sim, value = TRUE)
  expect_equal(length(invalid), 0L)
})

test_that("ROC simulated sequences differ from originals (seed 42, 1000 genes)", {
  sim <- sim_seqs(genome_roc)
  # With 1000 genes and random parameters, at least one sequence must change.
  expect_false(all(sim == original_seqs))
})

test_that("ROC simulated sequences are multiples of 3 (valid codon frames)", {
  sim <- sim_seqs(genome_roc)
  expect_true(all(nchar(sim) %% 3 == 0))
})

# ----------------------------------------------------------------------
# FONSE simulateGenome
# ----------------------------------------------------------------------
set.seed(42)
genome_fonse    <- initializeGenomeObject(file = fastaFile)
original_seqs_f <- genome_seqs(genome_fonse)

parameter_fonse <- initializeParameterObject(
  genome          = genome_fonse,
  sphi            = 1,
  num.mixtures    = 1,
  gene.assignment = rep(1, length(genome_fonse)),
  model           = "FONSE",
  init.initiation.cost = 4
)
model_fonse <- initializeModelObject(parameter_fonse, "FONSE")

test_that("FONSE simulateGenome returns NULL invisibly", {
  result <- model_fonse$simulateGenome(genome_fonse)
  expect_null(result)
})

test_that("FONSE simulateGenome preserves gene count", {
  expect_equal(length(genome_fonse), length(genome_roc))
})

test_that("FONSE simulated sequences have the same length as originals", {
  sim <- sim_seqs(genome_fonse)
  expect_equal(nchar(sim), original_lens)
})

test_that("FONSE simulated sequences contain only valid DNA bases", {
  sim <- sim_seqs(genome_fonse)
  invalid <- grep("[^ACGTacgt]", sim, value = TRUE)
  expect_equal(length(invalid), 0L)
})

test_that("FONSE simulated sequences are multiples of 3", {
  sim <- sim_seqs(genome_fonse)
  expect_true(all(nchar(sim) %% 3 == 0))
})
