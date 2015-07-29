library(testthat)
library(ribModel)

context("Gene")

g <- new(Gene)

test_that("set ID", {
  g$id = "blah"
  expect_equal(g$id, "blah")
})

test_that("set description", {
  g$description = "blah blah"
  expect_equal(g$description, "blah blah")
})

test_that("set sequence", {
  g$seq = "ATGCTCATTCTCACTGCTGCCTCGTAG"
  expect_equal(g$seq, "ATGCTCATTCTCACTGCTGCCTCGTAG")
})

ss <- new(SequenceSummary, "ATGCTCATTCTCACTGCTGCCTCGTAG")

test_that("get Sequence Summary", {
  ss2 <- g$getSequenceSummary()
  expect_equal(ss, ss2)
}) 



test_that("get Nucleotide At", {
  expect_equal(g$getNucleotideAt(0), "A")
  expect_equal(g$getNucleotideAt(1), "T")
  expect_equal(g$getNucleotideAt(2), "G")
  expect_equal(g$getNucleotideAt(3), "C")
  expect_equal(g$getNucleotideAt(4), "T")
  expect_equal(g$getNucleotideAt(5), "C")
  expect_equal(g$getNucleotideAt(6), "A")
  expect_equal(g$getNucleotideAt(7), "T")
  expect_equal(g$getNucleotideAt(8), "T")
  expect_equal(g$getNucleotideAt(9), "C")
  expect_equal(g$getNucleotideAt(10), "T")
  expect_equal(g$getNucleotideAt(11), "C")
  expect_equal(g$getNucleotideAt(12), "A")
  expect_equal(g$getNucleotideAt(13), "C")
  expect_equal(g$getNucleotideAt(14), "T")
  expect_equal(g$getNucleotideAt(15), "G")
  expect_equal(g$getNucleotideAt(16), "C")
  expect_equal(g$getNucleotideAt(17), "T")
  expect_equal(g$getNucleotideAt(18), "G")
  expect_equal(g$getNucleotideAt(19), "C")
  expect_equal(g$getNucleotideAt(20), "C")
  expect_equal(g$getNucleotideAt(21), "T")
  expect_equal(g$getNucleotideAt(22), "C")
  expect_equal(g$getNucleotideAt(23), "G")
  expect_equal(g$getNucleotideAt(24), "T")
  expect_equal(g$getNucleotideAt(25), "A")
  expect_equal(g$getNucleotideAt(26), "G")
})

test_that("Length", {
  expect_equal(g$length(), 27)
})

test_that("Not Nucleotide", {
  
})

test_that("clear", {
  g$clear()
  expect_equal(g$id, "")
  expect_equal(g$description, "")
  expect_equal(g$seq, "")
})