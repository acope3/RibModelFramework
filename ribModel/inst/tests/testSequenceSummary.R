library(testthat)
library(ribModel)

context("Sequence Summary")

ss <- new(SequenceSummary, 
"ATGCTCATTCTCACTGCTGCCTCGTAG"
)

test_that("AA counts for AA", {
	expect_equal(ss$getAACountForAA("M"), 1)
	expect_equal(ss$getAACountForAA("L"), 2)
	expect_equal(ss$getAACountForAA("I"), 1)
	expect_equal(ss$getAACountForAA("T"), 1)
	expect_equal(ss$getAACountForAA("A"), 2)
	expect_equal(ss$getAACountForAA("S"), 1)
	expect_equal(ss$getAACountForAA("X"), 1)
	expect_equal(ss$getAACountForAA("G"), 0)
})

test_that("AA counts for AA Index", {
  expect_equal(ss$getAACountForAAIndex(10), 1)
  expect_equal(ss$getAACountForAAIndex(9), 2)
  expect_equal(ss$getAACountForAAIndex(7), 1)
  expect_equal(ss$getAACountForAAIndex(16), 1)
  expect_equal(ss$getAACountForAAIndex(0), 2)
  expect_equal(ss$getAACountForAAIndex(15), 1)
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

test_that("RFP Observed test", {
  ss$SetRFPObserved(0, 45)
  ss$SetRFPObserved(1, 52)
  ss$setRFPObserved(2, 63)
  ss$setRFPObserved(60, 23)
  expect_equal(ss$getRFPObserved(0), 45)
  expect_equal(ss$getRFPObserved(1), 52)
  expect_equal(ss$getRFPObserved(2), 63)
  expect_equal(ss$getRFPObserved(60), 23)
})

test_that("Codon Positions", {
  # TODO: These return vectors. Fix it.
  expect_equal(ss$getCodonPositions(29), c(0))
  expect_equal(ss$getCodonPositions(24), c(1,3))
  expect_equal(ss$getCodonPositions(20), c(2))
  expect_equal(ss$getCodonPositions(51), c(4))
  expect_equal(ss$getCodonPositions(3), c(5))
  expect_equal(ss$getCodonPositions(1), c(6))
  expect_equal(ss$getCodonPositions(46), c(7))
  expect_equal(ss$getCodonPositions(62), c(8))
  expect_equal(ss$getCodonPositions(54), c())
})

test_that("AA to AA index", {
  expect_equal(ss$AAToAAIndex("A"), 0)
  expect_equal(ss$AAToAAIndex("C"), 1)
  expect_equal(ss$AAToAAIndex("D"), 2)
  expect_equal(ss$AAToAAIndex("E"), 3)
  expect_equal(ss$AAToAAIndex("F"), 4)
  expect_equal(ss$AAToAAIndex("G"), 5)
  expect_equal(ss$AAToAAIndex("H"), 6)
  expect_equal(ss$AAToAAIndex("I"), 7)
  expect_equal(ss$AAToAAIndex("K"), 8)
  expect_equal(ss$AAToAAIndex("L"), 9)
  expect_equal(ss$AAToAAIndex("M"), 10)
  expect_equal(ss$AAToAAIndex("N"), 11)
  expect_equal(ss$AAToAAIndex("P"), 12)
  expect_equal(ss$AAToAAIndex("Q"), 13)
  expect_equal(ss$AAToAAIndex("R"), 14)
  expect_equal(ss$AAToAAIndex("S"), 15)
  expect_equal(ss$AAToAAIndex("T"), 16)
  expect_equal(ss$AAToAAIndex("V"), 17)
  expect_equal(ss$AAToAAIndex("W"), 18)
  expect_equal(ss$AAToAAIndex("Y"), 19)
  expect_equal(ss$AAToAAIndex("Z"), 20)
  expect_equal(ss$AAToAAIndex("X"), 21)
})

test_that("Clear", {
  ss$clear()
  for (i in 0:64) {
    expect_equal(ss$getNumCodons(i), 0)
    expect_equal(ss$getRFPObserved(i), 0)
  }
  
  for (i in 0:22) {
    expect_equal(ss$getAAcountforAAIndex(i), 0)
  }
})

test_that("Process Sequence", {
  ss$processSequence("ATGCTCATTCTCACTGCTGCCTCGTAG")
  expect_equal(ss$getAACountForAA("I"), 1)
  expect_equal(ss$getAACountForAA("T"), 1)
  expect_equal(ss$getAACountForAAIndex(15), 1)
  expect_equal(ss$getAACountForAAIndex(21), 1)
  expect_equal(ss$getCodonCountForCodon("ATT"), 1)
  expect_equal(ss$getCodonCountForCodon("ACT"), 1)
  expect_equal(ss$getCodonCountForCodon("GCT"), 1)
  expect_equal(ss$getCodonCountForCodon("GCC"), 1)
  expect_equal(ss$getCodonPositions(24), c(1,3))
  expect_equal(ss$getCodonPositions(20), c(2))
})


