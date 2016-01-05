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
  for (codon in codons()){
    expect_equal(ss$getCodonCountForCodon(codon), ss2$getCodonCountForCodon(codon))
  }
  for (codon in codons()) {
    expect_equal(ss$getRFPObservedForCodon(codon), ss2$getRFPObservedForCodon(codon))
  }
  for (aa in aminoAcids()){
    expect_equal(ss$getAACountForAA(aa), ss2$getAACountForAA(aa))
  }
  for (codon in codons()) {
    expect_equal(ss$getCodonPositionsForCodon(codon), ss2$getCodonPositionsForCodon(codon))
  }
}) 

test_that("get Observed Phi Values", {
  expect_equal(g$getObservedPhiValues(), numeric())
  g$setObservedPhiValues(c(2.34,3.234,0.123))
  expect_equal(g$getObservedPhiValues(), c(2.34,3.234,0.123))
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

test_that("Reverse Complement", {
  expect_equal(g$reverseComplement()$seq, "CTACGAGGCAGCAGTGAGAATGAGCAT")
})

test_that("AA Sequence", {
  expect_equal(g$toAASequence(), "MLILTAASX")
})

test_that("clear", {
  g$clear()
  expect_equal(g$id, "")
  expect_equal(g$description, "")
  expect_equal(g$seq, "")
})


test_that("clean Sequence", {
  g$seq <- "ATGGTAACTTAG"
  g$cleanSeq()
  expect_equal(g$seq, "ATGGTAACTTAG")
  
#  g$seq <- "ATGGTAACTNNNQQQTAG"
#  g$cleanSeq()
#  expect_equal(g$seq, "ATGGTAACTNNNTAG")
})

g$clear()

g$seq <- "ATGCTCATTCTCACTGCTGCCTCGTAG"

test_that("get AA Count", {
  expect_equal(g$getAACount("M"), 1)
  expect_equal(g$getAACount("L"), 2)
  expect_equal(g$getAACount("I"), 1)
  expect_equal(g$getAACount("T"), 1)
  expect_equal(g$getAACount("A"), 2)
  expect_equal(g$getAACount("S"), 1)
  expect_equal(g$getAACount("X"), 1)
  expect_equal(g$getAACount("G"), 0)
})

test_that("get Codon Counts", {
  expect_equal(g$getCodonCount("ATG"), 1)
  expect_equal(g$getCodonCount("CTC"), 2)
  expect_equal(g$getCodonCount("ATT"), 1)
  expect_equal(g$getCodonCount("ACT"), 1)
  expect_equal(g$getCodonCount("GCT"), 1)
  expect_equal(g$getCodonCount("GCC"), 1)
  expect_equal(g$getCodonCount("TCG"), 1)
  expect_equal(g$getCodonCount("TAG"), 1)
  expect_equal(g$getCodonCount("AAA"), 0)
})

test_that("get RFP Observed", {
  g$setRFPObserved(4, 35)
  g$setRFPObserved(16, 45)
  g$setRFPObserved(54, 2)
  g$setRFPObserved(45, 0)
  expect_equal(g$getRFPObserved("TGC"), 35)
  expect_equal(g$getRFPObserved("CAC"), 45)
  expect_equal(g$getRFPObserved("GTG"),2)
  expect_equal(g$getRFPObserved("TCC"), 0)
})

test_that("get Codon Positions", {
  expect_equal(g$getCodonPositions("ATG"), c(0))
  expect_equal(g$getCodonPositions("CTC"), c(1,3))
  expect_equal(g$getCodonPositions("ATT"), c(2))
  expect_equal(g$getCodonPositions("ACT"), c(4))
  expect_equal(g$getCodonPositions("GCT"), c(5))
  expect_equal(g$getCodonPositions("GCC"), c(6))
  expect_equal(g$getCodonPositions("TCG"), c(7))
  expect_equal(g$getCodonPositions("TAG"), c(8))
  expect_equal(g$getCodonPositions("GTG"), numeric(0))
})


