library(testthat)
library(ribModel)

context("Codon Table")

CT <- new(CodonTable, 1, FALSE)

test_that("Get TableId",{
  for (i in 1:6)
  {
    CT <- new(CodonTable, i, FALSE)
    expect_equal(CT$getTableId(), i)
  }
  for (i in 9:14)
  {
    CT <- new(CodonTable, i, FALSE)
    expect_equal(CT$getTableId(), i)
  }
  for (i in 21:25)
  {
    CT <- new(CodonTable, i, FALSE)
    expect_equal(CT$getTableId(), i)
  }
  for (i in 7:8)
  {
    CT <- new(CodonTable, i, FALSE)
    expect_equal(CT$getTableId(), 1)
  }
  CT <- new(CodonTable, 15, FALSE)
  expect_equal(CT$getTableId(), 1)
  CT <- new(CodonTable, 16, FALSE)
  expect_equal(CT$getTableId(), 16)
  for (i in 17:20)
  {
    CT <- new(CodonTable, i, FALSE)
    expect_equal(CT$getTableId(), 1)
  }
})

CT <- new(CodonTable, 1, FALSE)
CT$setupCodonTable()
test_that("Get SplitAA",{
  expect_equal(CT$getSplitAA(), FALSE)
  otherCT <- new(CodonTable, 1, TRUE)
  expect_equal(otherCT$getSplitAA(), TRUE)
})

test_that("Default Constructor",{
  defaultCT <- new(CodonTable)
  expect_equal(defaultCT$getTableId(), 1)
  expect_equal(defaultCT$getSplitAA(), TRUE)
})

test_that("Get Codon Array",{
  codonArray <- getCodonArray()
  myArray <- c("GCA", "GCC", "GCG", "GCT", "TGC", "TGT", "GAC", "GAT", "GAA", "GAG",
               "TTC", "TTT", "GGA", "GGC", "GGG", "GGT", "CAC", "CAT", "ATA", "ATC",
               "ATT", "AAA", "AAG", "CTA", "CTC", "CTG", "CTT", "TTA", "TTG", "ATG",
               "AAC", "AAT", "CCA", "CCC", "CCG", "CCT", "CAA", "CAG", "AGA", "AGG",
               "CGA", "CGC", "CGG", "CGT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC",
               "ACG", "ACT", "GTA", "GTC", "GTG", "GTT", "TGG", "TAC", "TAT", "AGC",
               "AGT", "TAA", "TAG", "TGA")
  expect_equal(codonArray, myArray)
})


test_that("Index to Codon",{
  codonArray <- getCodonArray()
  for (index in 1:64) {
    expect_equal(CT$indexToCodon(index, FALSE), codonArray[index])
  }
})