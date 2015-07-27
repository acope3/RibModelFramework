library(testthat)
library(ribModel)

context("Sequence Summary")

ss <- new(SequenceSummary, 
"ATGCTCATTCTCACTGCTGCCTCGTAG"
)

test_that("AA counts for AA", {
	expect_equal(ss$getAAcount("M"), 1)
	expect_equal(ss$getAAcount("L"), 2)
	expect_equal(ss$getAAcount("I"), 1)
	expect_equal(ss$getAAcount("T"), 1)
	expect_equal(ss$getAAcount("A"), 2)
	expect_equal(ss$getAAcount("S"), 1)
	expect_equal(ss$getAAcount("X"), 1)
	expect_equal(ss$getAAcount("G"), 0)
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