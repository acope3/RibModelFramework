library(testthat)
library(ribModel)

context("Codon Table")

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


test_that("Index to Codon",{
  codonArray <- getCodonArray()
  for (index in 1:64) {
    expect_equal(CT$indexToCodon(index, FALSE), codonArray[index])
  }
})


test_that("Get Ser2",{
  expect_equal(getSer2(), "Z")
})


test_that("Get Ser1",{
  expect_equal(getSer1(), "J")
})


test_that("Get Thr4_1",{
  expect_equal(getThr4_1(), "T")
})


test_that("Get Thr4_2",{
  expect_equal(getThr4_2(), "B")
})


test_that("Get Leu1",{
  expect_equal(getLeu1(), "U")
})


test_that("Get Amino Acid Array",{
  AAList <- getAminoAcidArray()
  trueList <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "U", "M", "N", "P", "Q",
                "R", "J", "Z", "S", "T", "T", "B", "V", "W", "Y", "X")
  expect_equal(AAList, trueList)
})


test_that("Get Amino Acid Array Without Split",{
  AAList <- getAminoAcidArrayWithoutSplit()
  trueList <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q",
                "R", "S", "T", "V", "W", "Y", "X")
  expect_equal(AAList, trueList)
})


test_that("Get Num Codons Per AA For Table",{
  numCodonsList <- getNumCodonsPerAAForTable()
  trueList <- list(c(4,2,2,2,2,4,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,3))
  trueList[[2]] <- c(4,2,2,2,2,4,2,2,2,6,0,2,2,4,2,4,0,2,4,4,0,0,4,2,2,4)
  trueList[[3]] <- c(4,2,2,2,2,4,2,2,2,2,0,2,2,4,2,6,0,2,4,0,4,4,4,2,2,2)
  trueList[[4]] <- c(4,2,2,2,2,4,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,2,2,2)
  trueList[[5]] <- c(4,2,2,2,2,4,2,2,2,6,0,2,2,4,2,4,0,4,4,4,0,0,4,2,2,2)
  trueList[[6]] <- c(4,2,2,2,2,4,2,3,2,6,0,1,2,4,4,6,0,2,4,4,0,0,4,1,2,1)
  trueList[[7]] <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  trueList[[8]] <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  trueList[[9]] <- c(4,2,2,2,2,4,2,3,1,6,0,1,3,4,2,4,0,4,4,4,0,0,4,2,2,2)
  trueList[[10]] <- c(4,3,2,2,2,4,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,2)
  trueList[[11]] <- c(4,2,2,2,2,4,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,3)
  trueList[[12]] <- c(4,2,2,2,2,4,2,3,2,5,0,1,2,4,2,6,1,2,4,4,0,0,4,1,2,3)
  trueList[[13]] <- c(4,2,2,2,2,6,2,2,2,6,0,2,2,4,2,4,0,2,4,4,0,0,4,2,2,2)
  trueList[[14]] <- c(4,2,2,2,2,4,2,3,1,6,0,1,3,4,2,4,0,4,4,4,0,0,4,2,3,1)
  trueList[[15]] <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  trueList[[16]] <- c(4,2,2,2,2,4,2,3,2,6,1,1,2,4,2,6,0,2,4,4,0,0,4,1,2,2)
  trueList[[17]] <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  trueList[[18]] <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  trueList[[19]] <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  trueList[[20]] <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  trueList[[21]] <- c(4,2,2,2,2,4,2,2,1,6,0,2,3,4,2,4,0,4,4,4,0,0,4,2,2,2)
  trueList[[22]] <- c(4,2,2,2,2,4,3,2,2,6,1,1,2,4,2,6,0,2,3,4,0,0,4,1,2,3)
  trueList[[23]] <- c(4,2,2,2,2,4,2,3,2,5,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,4)
  trueList[[24]] <- c(4,2,2,2,2,4,2,3,3,6,0,1,2,4,2,4,0,3,4,4,0,0,4,2,2,2)
  trueList[[25]] <- c(4,2,2,2,2,5,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,2)
  expect_equal(numCodonsList[[3]], trueList[[3]])
})


test_that("Get Codon Table Definition",{
  definitionList <- getCodonTableDefinition()
  trueList <- c("1. The Standard Code","2. The Vertebrate Mitochondrial Code",
                "3. The Yeast Mitochondrial Code",
                "4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code",
                "5. The Invertebrate Mitochondrial Code",
                "6. The Ciliate, Dasycladacean and Hexamita Nuclear Code",
                "7. Invalid Codon Table","8. Invalid Codon Table",
                "9. The Echinoderm and Flatworm Mitochondrial Code",
                "10. The Euplotid Nuclear Code",
                "11. The Bacterial, Archaeal and Plant Plastid Code",
                "12. The Alternative Yeast Nuclear Code",
                "13. The Ascidian Mitochondrial Code",
                "14. The Alternative Flatworm Mitochondrial Code",
                "15. Invalid Codon Table",
                "16. Chlorophycean Mitochondrial Code",
                "17. Invalid Codon Table","18. Invalid Codon Table",
                "19. Invalid Codon Table","20. Invalid Codon Table",
                "21. Trematode Mitochondrial Code", "22. Scenedesmus obliquus Mitochondrial Code",
                "23. Thraustochytrium Mitochondrial Code", "24. Pterobranchia Mitochondrial Code",
                "25. Candidate Division SR1 and Gracilibacteria Code")
  expect_equal(definitionList, trueList)
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
