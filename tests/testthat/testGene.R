library(testthat)
library(ribModel)

context("Gene")

g <- new(Gene)

test_that("set & get ID", {
  g$id = "blah"
  expect_equal(g$id, "blah")
})

test_that("set & get description", {
  g$description = "blah blah"
  expect_equal(g$description, "blah blah")
})

test_that("set & get sequence", {
  g$seq = "ATGCTCATTCTCACTGCTGCCTCGTAG"
  expect_equal(g$seq, "ATGCTCATTCTCACTGCTGCCTCGTAG")
})


test_that("Length", {
  expect_equal(g$length(), 27)
})

g <- new(Gene, "ATGCTCATTCTCACTGCTGCCTCGTAG", "2", "New Test Gene")
#TODO: problem. This used to say g$seq <- "kjdklsjfkdj" (string). It
#does set the string correctly, but the values in SS are NOT CLEARED.
test_that("get AA Count", {
  expect_equal(g$getAACount("M"), 1)
  expect_equal(g$getAACount("L"), 2)
  expect_equal(g$getAACount("I"), 1)
  expect_equal(g$getAACount("T"), 1)
  expect_equal(g$getAACount("A"), 2)
  expect_equal(g$getAACount("S"), 1)
  expect_equal(g$getAACount("X"), 1)
  expect_equal(g$getAACount("G"), 0)
  
  #Checking invalid cases
  expect_equal(g$getAACount("g"), 0)
  expect_equal(g$getAACount("AA"), 0)
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
  
  #Checking invalid cases
  expect_equal(g$getCodonCount("atg"), 0)
  expect_equal(g$getCodonCount("ATGG"), 0)
})

test_that("get RFP Observed", {
  expect_equal(g$getRFPObserved("TGC"), 0)
  expect_equal(g$getRFPObserved("CAC"), 0)
  expect_equal(g$getRFPObserved("GTG"),0)
  expect_equal(g$getRFPObserved("TCC"), 0)
  
  #Checking invalid cases
  expect_equal(g$getRFPObserved("atg"), 0)
  expect_equal(g$getRFPObserved("ATGG"), 0)
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
  
  #Checking invalid cases
  expect_equal(g$getCodonPositions("atg"), numeric(0))
  expect_equal(g$getCodonPositions("ATGG"), numeric(0))
})


#TODO NOTE:
#See if there is a way to expect error messages instead of retrun values, or 
#if you can check for both.