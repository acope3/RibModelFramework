library(testthat)
library(ribModel)

context("Sequence Summary")

createCodonTable(1,TRUE)

ss <- new(SequenceSummary, 
"ATGCTCATTCTCACTGCTGCCTCGTAG"
)

test_that("Get AA counts for AA", {
	expect_equal(ss$getAACountForAA("M"), 1)
	expect_equal(ss$getAACountForAA("L"), 2)
	expect_equal(ss$getAACountForAA("I"), 1)
	expect_equal(ss$getAACountForAA("T"), 1)
	expect_equal(ss$getAACountForAA("A"), 2)
	expect_equal(ss$getAACountForAA("S"), 1)
	expect_equal(ss$getAACountForAA("X"), 1)
	expect_equal(ss$getAACountForAA("G"), 0)
})

test_that("Get AA counts for AA Index", {
  expect_equal(ss$getAACountForAAIndex(10), 1)
  expect_equal(ss$getAACountForAAIndex(9), 2)
  expect_equal(ss$getAACountForAAIndex(7), 1)
  expect_equal(ss$getAACountForAAIndex(16), 1)
  expect_equal(ss$getAACountForAAIndex(0), 2)
  expect_equal(ss$getAACountForAAIndex(16), 1)
  expect_equal(ss$getAACountForAAIndex(21), 1)
  expect_equal(ss$getAACountForAAIndex(2), 0)
})

test_that("Codon Counts for Codon", {
  expect_equal(ss$getCodonCountForCodon("ATG"), 1)
  expect_equal(ss$getCodonCountForCodon("CTC"), 2)
  expect_equal(ss$getCodonCountForCodon("ATT"), 1)
  expect_equal(ss$getCodonCountForCodon("ACT"), 1)
  expect_equal(ss$getCodonCountForCodon("GCT"), 1)
  expect_equal(ss$getCodonCountForCodon("GCC"), 1)
  expect_equal(ss$getCodonCountForCodon("TCG"), 1)
  expect_equal(ss$getCodonCountForCodon("TAG"), 1)
  expect_equal(ss$getCodonCountForCodon("AAA"), 0)
})

test_that("Codon Counts for Codon Index", {
  expect_equal(ss$getCodonCountForCodonIndex(29), 1)
  expect_equal(ss$getCodonCountForCodonIndex(24), 2)
  expect_equal(ss$getCodonCountForCodonIndex(20), 1)
  expect_equal(ss$getCodonCountForCodonIndex(51), 1)
  expect_equal(ss$getCodonCountForCodonIndex(3), 1)
  expect_equal(ss$getCodonCountForCodonIndex(1), 1)
  expect_equal(ss$getCodonCountForCodonIndex(46), 1)
  expect_equal(ss$getCodonCountForCodonIndex(62), 1)
  expect_equal(ss$getCodonCountForCodonIndex(2), 0)
})

test_that("RFP Observed for codon", {
  ss$setRFPObserved(4, 35)
  ss$setRFPObserved(16, 45)
  ss$setRFPObserved(54, 2)
  ss$setRFPObserved(45, 0)
  expect_equal(ss$getRFPObservedForCodon("TGC"), 35)
  expect_equal(ss$getRFPObservedForCodon("CAC"), 45)
  expect_equal(ss$getRFPObservedForCodon("GTG"),2)
  expect_equal(ss$getRFPObservedForCodon("TCC"), 0)
})

test_that("RFP Observed test for codon Index", {
  ss$setRFPObserved(0, 45)
  ss$setRFPObserved(1, 52)
  ss$setRFPObserved(2, 63)
  ss$setRFPObserved(60, 23)
  expect_equal(ss$getRFPObservedForCodonIndex(0), 45)
  expect_equal(ss$getRFPObservedForCodonIndex(1), 52)
  expect_equal(ss$getRFPObservedForCodonIndex(2), 63)
  expect_equal(ss$getRFPObservedForCodonIndex(60), 23)
})

test_that("Codon Positions by Codon", {
  expect_equal(ss$getCodonPositionsForCodon("ATG"), c(0))
  expect_equal(ss$getCodonPositionsForCodon("CTC"), c(1,3))
  expect_equal(ss$getCodonPositionsForCodon("ATT"), c(2))
  expect_equal(ss$getCodonPositionsForCodon("ACT"), c(4))
  expect_equal(ss$getCodonPositionsForCodon("GCT"), c(5))
  expect_equal(ss$getCodonPositionsForCodon("GCC"), c(6))
  expect_equal(ss$getCodonPositionsForCodon("TCG"), c(7))
  expect_equal(ss$getCodonPositionsForCodon("TAG"), c(8))
  expect_equal(ss$getCodonPositionsForCodon("GTG"), numeric(0))
})

test_that("Codon Positions", {
  expect_equal(ss$getCodonPositionsForCodonIndex(29), c(0))
  expect_equal(ss$getCodonPositionsForCodonIndex(24), c(1,3))
  expect_equal(ss$getCodonPositionsForCodonIndex(20), c(2))
  expect_equal(ss$getCodonPositionsForCodonIndex(51), c(4))
  expect_equal(ss$getCodonPositionsForCodonIndex(3), c(5))
  expect_equal(ss$getCodonPositionsForCodonIndex(1), c(6))
  expect_equal(ss$getCodonPositionsForCodonIndex(46), c(7))
  expect_equal(ss$getCodonPositionsForCodonIndex(62), c(8))
  expect_equal(ss$getCodonPositionsForCodonIndex(54), numeric(0))
})

test_that("Clear", {
  ss$clear()
  
  # i-1 is because of the difference in indexing from R to C++
  
  for (i in 1:64) {
    expect_equal(ss$getCodonCountForCodonIndex(i-1), 0)
    expect_equal(ss$getRFPObservedForCodonIndex(i-1), 0)
  }
  CT <- getInstance()
  map <- CT$getAAMap()
  size <- length(map)
  for (i in 1:size) {
    expect_equal(ss$getAACountForAAIndex(i-1), 0)
  }
})

ss$processSequence("ATGCTCATTCTCACTGCTGCCTCGTAG")

test_that("Process Sequence", {
  
  expect_equal(ss$getAACountForAA("I"), 1)
  expect_equal(ss$getAACountForAA("T"), 1)
  expect_equal(ss$getAACountForAAIndex(16), 1)
  expect_equal(ss$getAACountForAAIndex(21), 1)
  expect_equal(ss$getCodonCountForCodon("ATT"), 1)
  expect_equal(ss$getCodonCountForCodon("ACT"), 1)
  expect_equal(ss$getCodonCountForCodon("GCT"), 1)
  expect_equal(ss$getCodonCountForCodon("GCC"), 1)
  expect_equal(ss$getCodonPositionsForCodonIndex(24), c(1,3))
  expect_equal(ss$getCodonPositionsForCodonIndex(20), c(2))
})


test_that("Complement Nucleotide", {
  expect_equal(complimentNucleotide("A"), "T")
  expect_equal(complimentNucleotide("T"), "A")
  expect_equal(complimentNucleotide("C"), "G")
  expect_equal(complimentNucleotide("G"), "C")
  expect_equal(complimentNucleotide("Q"), "C")
})

