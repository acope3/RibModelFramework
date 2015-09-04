library(testthat)
library(ribModel)

context("Codon Table")

CT <- new(CodonTable, 1, FALSE)

test_that("Get TableId",{
  expect_equal(CT$getTableId(), 1)
  badCT <- new(CodonTable, 7, FALSE)
  expect_equal(badCT$getTableId(), 1)
  #ToDo: check all values
})

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

CT$setupCodonTable()
codonIndexList <- list(c(1,2,3,4))
codonIndexList[[2]] <- c(5,6)
codonIndexList[[3]] <- c(7,8)
codonIndexList[[4]] <- c(9,10)
codonIndexList[[5]] <- c(11,12)
codonIndexList[[6]] <- c(13,14,15,16)
codonIndexList[[7]] <- c(17,18)
codonIndexList[[8]] <- c(19,20,21)
codonIndexList[[9]] <- c(22,23)
codonIndexList[[10]] <- c(24,25,26,27,28,29)
codonIndexList[[11]] <- c(30)
codonIndexList[[12]] <- c(31,32)
codonIndexList[[13]] <- c(33,34,35,36)
codonIndexList[[14]] <- c(37,38)
codonIndexList[[15]] <- c(39,40,41,42,43,44)
codonIndexList[[16]] <- c(45,46,47,48,60,61)
codonIndexList[[17]] <- c(49,50,51,52)
codonIndexList[[18]] <- c(53,54,55,56)
codonIndexList[[19]] <- c(57)
codonIndexList[[20]] <- c(58,59)
codonIndexListWithoutRef <- list(c(1,2,3))
codonIndexListWithoutRef[[2]] <- c(5)
codonIndexListWithoutRef[[3]] <- c(7)
codonIndexListWithoutRef[[4]] <- c(9)
codonIndexListWithoutRef[[5]] <- c(11)
codonIndexListWithoutRef[[6]] <- c(13,14,15)
codonIndexListWithoutRef[[7]] <- c(17)
codonIndexListWithoutRef[[8]] <- c(19,20)
codonIndexListWithoutRef[[9]] <- c(22)
codonIndexListWithoutRef[[10]] <- c(24,25,26,27,28)
codonIndexListWithoutRef[[11]] <- numeric(0)
codonIndexListWithoutRef[[12]] <- c(31)
codonIndexListWithoutRef[[13]] <- c(33,34,35)
codonIndexListWithoutRef[[14]] <- c(37)
codonIndexListWithoutRef[[15]] <- c(39,40,41,42,43)
codonIndexListWithoutRef[[16]] <- c(45,46,47,48,60)
codonIndexListWithoutRef[[17]] <- c(49,50,51)
codonIndexListWithoutRef[[18]] <- c(53,54,55)
codonIndexListWithoutRef[[19]] <- numeric(0)
codonIndexListWithoutRef[[20]] <- c(58)


test_that("Get Codon Index Listing For Table 1 (no split)",{
  listing <- CT$getCodonIndexListing()
  expect_equal(listing, codonIndexList)
})



test_that("Get Codon Index Listing Without Reference For Table 1 (no split)", {
  listing <- CT$getCodonIndexListingWithoutReference()
  expect_equal(listing, codonIndexListWithoutRef)
})

test_that("Get AA Listing For Table 1 (no split)",{
  listing <- CT$getAAListing()
  trueList <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R",
                "S", "T", "V", "W", "Y")
  expect_equal(listing, trueList)
})

test_that("Get For Param Vector Listing for Table 1 (no split)",{
  listing <- CT$getForParamVectorListing()
  trueList <- c("GCA", "GCC", "GCG", "TGC", "GAC", "GAA", "TTC", "GGA", "GGC", "GGG", "CAC",
                "ATA", "ATC", "AAA", "CTA", "CTC", "CTG", "CTT", "TTA", "AAC", "CCA", "CCC",
                "CCG", "CAA", "AGA", "AGG", "CGA", "CGC", "CGG", "TCA", "TCC", "TCG", "TCT",
                "AGC", "ACA", "ACC", "ACG", "GTA", "GTC", "GTG", "TAC")
  expect_equal(listing, trueList)
})

test_that("Get Codon To AA Map for Table 1 (no split)",{
  map <- CT$getCodonToAAMap()
  trueMap <- c(AAA = "K", AAC = "N", AAG = "K", AAT = "N", ACA = "T", ACC = "T", ACG = "T",
               ACT = "T", AGA = "R", AGC = "S", AGG = "R", AGT = "S", ATA = "I", ATC = "I",
               ATG = "M", ATT = "I", CAA = "Q", CAC = "H", CAG = "Q", CAT = "H", CCA = "P",
               CCC = "P", CCG = "P", CCT = "P", CGA = "R", CGC = "R", CGG = "R", CGT = "R", 
               CTA = "L", CTC = "L", CTG = "L", CTT = "L", GAA = "E", GAC = "D", GAG = "E", 
               GAT = "D", GCA = "A", GCC = "A", GCG = "A", GCT = "A", GGA = "G", GGC = "G", 
               GGG = "G", GGT = "G", GTA = "V", GTC = "V", GTG = "V", GTT = "V", TAC = "Y", 
               TAT = "Y", TCA = "S", TCC = "S", TCG = "S", TCT = "S", TGC = "C", TGG = "W", 
               TGT = "C", TTA = "L", TTC = "F", TTG = "L", TTT = "F")
  expect_equal(map, trueMap)
})

test_that("Get AA Map for Table 1 (no split)",{
  map <- CT$getAAMap()
  trueMap <- c(A = 1, C = 2, D = 3, E = 4, F = 5, G = 6, H = 7, I = 8, K = 9, L = 10, M = 11,
               N = 12, P = 13, Q = 14, R = 15, S = 16, T = 17, V = 18, W = 19, Y = 20)
  expect_equal(map, trueMap)
})

test_that("Get AA To Num Codons Map for Table 1 (no split)",{
  map <- CT$getAAToNumCodonsMap()
  trueMap <- c(A = 4, C = 2, D = 2, E = 2, F = 2, G = 4, H = 2, I = 3, K = 2, L = 6, M = 1,
               N = 2, P = 4, Q = 2, R = 6, S = 6, T = 4, V = 4, W = 1, Y = 2)
  expect_equal(map,trueMap)
})

test_that("Get Num Codons For AA for Table 1 (no split)",{
  expect_equal(CT$getNumCodonsForAA("A", FALSE), 4)
  expect_equal(CT$getNumCodonsForAA("C", FALSE), 2)
  expect_equal(CT$getNumCodonsForAA("D", FALSE), 2)
  expect_equal(CT$getNumCodonsForAA("E", FALSE), 2)
  expect_equal(CT$getNumCodonsForAA("F", FALSE), 2)
  expect_equal(CT$getNumCodonsForAA("G", FALSE), 4)
  expect_equal(CT$getNumCodonsForAA("H", FALSE), 2)
  expect_equal(CT$getNumCodonsForAA("I", FALSE), 3)
  expect_equal(CT$getNumCodonsForAA("K", FALSE), 2)
  expect_equal(CT$getNumCodonsForAA("L", FALSE), 6)
  expect_equal(CT$getNumCodonsForAA("M", FALSE), 1)
  expect_equal(CT$getNumCodonsForAA("N", FALSE), 2)
  expect_equal(CT$getNumCodonsForAA("P", FALSE), 4)
  expect_equal(CT$getNumCodonsForAA("Q", FALSE), 2)
  expect_equal(CT$getNumCodonsForAA("R", FALSE), 6)
  expect_equal(CT$getNumCodonsForAA("S", FALSE), 6)
  expect_equal(CT$getNumCodonsForAA("T", FALSE), 4)
  expect_equal(CT$getNumCodonsForAA("V", FALSE), 4)
  expect_equal(CT$getNumCodonsForAA("W", FALSE), 1)
  expect_equal(CT$getNumCodonsForAA("Y", FALSE), 2)
  
  expect_equal(CT$getNumCodonsForAA("A", TRUE), 3)
  expect_equal(CT$getNumCodonsForAA("C", TRUE), 1)
  expect_equal(CT$getNumCodonsForAA("D", TRUE), 1)
  expect_equal(CT$getNumCodonsForAA("E", TRUE), 1)
  expect_equal(CT$getNumCodonsForAA("F", TRUE), 1)
  expect_equal(CT$getNumCodonsForAA("G", TRUE), 3)
  expect_equal(CT$getNumCodonsForAA("H", TRUE), 1)
  expect_equal(CT$getNumCodonsForAA("I", TRUE), 2)
  expect_equal(CT$getNumCodonsForAA("K", TRUE), 1)
  expect_equal(CT$getNumCodonsForAA("L", TRUE), 5)
  expect_equal(CT$getNumCodonsForAA("M", TRUE), 0)
  expect_equal(CT$getNumCodonsForAA("N", TRUE), 1)
  expect_equal(CT$getNumCodonsForAA("P", TRUE), 3)
  expect_equal(CT$getNumCodonsForAA("Q", TRUE), 1)
  expect_equal(CT$getNumCodonsForAA("R", TRUE), 5)
  expect_equal(CT$getNumCodonsForAA("S", TRUE), 5)
  expect_equal(CT$getNumCodonsForAA("T", TRUE), 3)
  expect_equal(CT$getNumCodonsForAA("V", TRUE), 3)
  expect_equal(CT$getNumCodonsForAA("W", TRUE), 0)
  expect_equal(CT$getNumCodonsForAA("Y", TRUE), 1)
})

#TODO: index wrapper needed
test_that("Get Num Codons For AA Index",{
  expect_equal(CT$getNumCodonsForAAIndex(1, FALSE), 4)
  expect_equal(CT$getNumCodonsForAAIndex(2, FALSE), 2)
  expect_equal(CT$getNumCodonsForAAIndex(3, FALSE), 2)
  expect_equal(CT$getNumCodonsForAAIndex(4, FALSE), 2)
  expect_equal(CT$getNumCodonsForAAIndex(5, FALSE), 2)
  expect_equal(CT$getNumCodonsForAAIndex(6, FALSE), 4)
  expect_equal(CT$getNumCodonsForAAIndex(7, FALSE), 2)
  expect_equal(CT$getNumCodonsForAAIndex(8, FALSE), 3)
  expect_equal(CT$getNumCodonsForAAIndex(9, FALSE), 2)
  expect_equal(CT$getNumCodonsForAAIndex(10, FALSE), 6)
  expect_equal(CT$getNumCodonsForAAIndex(11, FALSE), 1)
  expect_equal(CT$getNumCodonsForAAIndex(12, FALSE), 2)
  expect_equal(CT$getNumCodonsForAAIndex(13, FALSE), 4)
  expect_equal(CT$getNumCodonsForAAIndex(14, FALSE), 2)
  expect_equal(CT$getNumCodonsForAAIndex(15, FALSE), 6)
  expect_equal(CT$getNumCodonsForAAIndex(16, FALSE), 6)
  expect_equal(CT$getNumCodonsForAAIndex(17, FALSE), 4)
  expect_equal(CT$getNumCodonsForAAIndex(18, FALSE), 4)
  expect_equal(CT$getNumCodonsForAAIndex(19, FALSE), 1)
  expect_equal(CT$getNumCodonsForAAIndex(20, FALSE), 2)
  
  expect_equal(CT$getNumCodonsForAAIndex(1, TRUE), 3)
  expect_equal(CT$getNumCodonsForAAIndex(2, TRUE), 1)
  expect_equal(CT$getNumCodonsForAAIndex(3, TRUE), 1)
  expect_equal(CT$getNumCodonsForAAIndex(4, TRUE), 1)
  expect_equal(CT$getNumCodonsForAAIndex(5, TRUE), 1)
  expect_equal(CT$getNumCodonsForAAIndex(6, TRUE), 3)
  expect_equal(CT$getNumCodonsForAAIndex(7, TRUE), 1)
  expect_equal(CT$getNumCodonsForAAIndex(8, TRUE), 2)
  expect_equal(CT$getNumCodonsForAAIndex(9, TRUE), 1)
  expect_equal(CT$getNumCodonsForAAIndex(10, TRUE), 5)
  expect_equal(CT$getNumCodonsForAAIndex(11, TRUE), 0)
  expect_equal(CT$getNumCodonsForAAIndex(12, TRUE), 1)
  expect_equal(CT$getNumCodonsForAAIndex(13, TRUE), 3)
  expect_equal(CT$getNumCodonsForAAIndex(14, TRUE), 1)
  expect_equal(CT$getNumCodonsForAAIndex(15, TRUE), 5)
  expect_equal(CT$getNumCodonsForAAIndex(16, TRUE), 5)
  expect_equal(CT$getNumCodonsForAAIndex(17, TRUE), 3)
  expect_equal(CT$getNumCodonsForAAIndex(18, TRUE), 3)
  expect_equal(CT$getNumCodonsForAAIndex(19, TRUE), 0)
  expect_equal(CT$getNumCodonsForAAIndex(20, TRUE), 1)
})

test_that("Get For Param Vector Codon",{
  index <- 1
  list <- CT$getForParamVectorListing()
  while (index < length(list)) {
    expect_equal(CT$getForParamVectorCodon(index), list[index])
    index <- index + 1
  }
})

test_that("AA To AA Index",{
  AAList <- CT$getAAListing()
  index <- 1
  for (aa in AAList) {
    expect_equal(CT$AAToAAIndex(aa), index)
    index <- index + 1
  }
  expect_equal(CT$AAToAAIndex("a"), 1) #prove that lower case AA works
  expect_equal(CT$AAToAAIndex("bob"), 0) #prove that invalid cases (nonsense) return 0
  expect_equal(CT$AAToAAIndex("J"), 0) #Valid for some other codon tables, but not 1 
})

test_that("AA Index To Codon Range",{
  for (i in 1:20){
    expect_equal(CT$AAIndexToCodonRange(i, FALSE), codonIndexList[[i]])
  }
  for (i in 1:20){
    expect_equal(CT$AAIndexToCodonRange(i,TRUE), codonIndexListWithoutRef[[i]]) #problem in C++ function 
  }
})