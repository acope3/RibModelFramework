library(testthat)
library(ribModel)

context("Codon Table")
CT <- new(CodonTable, 1, FALSE)
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
    expect_equal(CT$AAIndexToCodonRange(i,TRUE), codonIndexListWithoutRef[[i]])
  }
  
  expect_equal(CT$AAIndexToCodonRange(21, FALSE), numeric(0)) #Test with bad index
  expect_equal(CT$AAIndexToCodonRange(0, FALSE), numeric(0)) #Test with bad index
  expect_equal(CT$AAIndexToCodonRange(-1, FALSE), numeric(0)) #Test with bad index
})

test_that("Index to Codon",{
  expect_equal(CT$indexToCodon(1, FALSE), "GCA")
})

test_that("AA To Codon Range",{
  aaListing <- CT$getAAListing()
  i <- 1
  for (aa in aaListing){
    expect_equal(CT$AAToCodonRange(aa, FALSE), codonIndexList[[i]])
    i <- i + 1
  }
  i <- 1
  for (aa in aaListing){
    expect_equal(CT$AAToCodonRange(aa,TRUE), codonIndexListWithoutRef[[i]])
    i <- i + 1
  }
  
  expect_equal(CT$AAToCodonRange("a", FALSE), codonIndexList[[1]]) #Lower case AA test
  expect_equal(CT$AAToCodonRange("BOB", FALSE), numeric(0)) #Test with bad AA given
})

test_that("AA To Codon",{
  expect_equal(CT$AAToCodon("A", FALSE), c("GCA", "GCC", "GCG", "GCT"))
  expect_equal(CT$AAToCodon("C", FALSE), c("TGC", "TGT"))
  expect_equal(CT$AAToCodon("D", FALSE), c("GAC", "GAT"))
  expect_equal(CT$AAToCodon("E", FALSE), c("GAA", "GAG"))
  expect_equal(CT$AAToCodon("F", FALSE), c("TTC", "TTT"))
  expect_equal(CT$AAToCodon("G", FALSE), c("GGA", "GGC", "GGG", "GGT"))
  expect_equal(CT$AAToCodon("H", FALSE), c("CAC", "CAT"))
  expect_equal(CT$AAToCodon("I", FALSE), c("ATA", "ATC", "ATT"))
  expect_equal(CT$AAToCodon("K", FALSE), c("AAA", "AAG"))
  expect_equal(CT$AAToCodon("L", FALSE), c("CTA", "CTC", "CTG", "CTT", "TTA", "TTG"))
  expect_equal(CT$AAToCodon("M", FALSE), c("ATG"))
  expect_equal(CT$AAToCodon("N", FALSE), c("AAC", "AAT"))
  expect_equal(CT$AAToCodon("P", FALSE), c("CCA", "CCC", "CCG", "CCT"))
  expect_equal(CT$AAToCodon("Q", FALSE), c("CAA", "CAG"))
  expect_equal(CT$AAToCodon("R", FALSE), c("AGA", "AGG", "CGA", "CGC", "CGG", "CGT"))
  expect_equal(CT$AAToCodon("S", FALSE), c("TCA", "TCC", "TCG", "TCT", "AGC", "AGT"))
  expect_equal(CT$AAToCodon("T", FALSE), c("ACA", "ACC", "ACG", "ACT"))
  expect_equal(CT$AAToCodon("V", FALSE), c("GTA", "GTC", "GTG", "GTT"))
  expect_equal(CT$AAToCodon("W", FALSE), c("TGG"))
  expect_equal(CT$AAToCodon("Y", FALSE), c("TAC", "TAT"))
  
  
  expect_equal(CT$AAToCodon("A", TRUE), c("GCA", "GCC", "GCG"))
  expect_equal(CT$AAToCodon("C", TRUE), c("TGC"))
  expect_equal(CT$AAToCodon("D", TRUE), c("GAC"))
  expect_equal(CT$AAToCodon("E", TRUE), c("GAA"))
  expect_equal(CT$AAToCodon("F", TRUE), c("TTC"))
  expect_equal(CT$AAToCodon("G", TRUE), c("GGA", "GGC", "GGG"))
  expect_equal(CT$AAToCodon("H", TRUE), c("CAC"))
  expect_equal(CT$AAToCodon("I", TRUE), c("ATA", "ATC"))
  expect_equal(CT$AAToCodon("K", TRUE), c("AAA"))
  expect_equal(CT$AAToCodon("L", TRUE), c("CTA", "CTC", "CTG", "CTT", "TTA"))
  expect_equal(CT$AAToCodon("M", TRUE), character(0))
  expect_equal(CT$AAToCodon("N", TRUE), c("AAC"))
  expect_equal(CT$AAToCodon("P", TRUE), c("CCA", "CCC", "CCG"))
  expect_equal(CT$AAToCodon("Q", TRUE), c("CAA"))
  expect_equal(CT$AAToCodon("R", TRUE), c("AGA", "AGG", "CGA", "CGC", "CGG"))
  expect_equal(CT$AAToCodon("S", TRUE), c("TCA", "TCC", "TCG", "TCT", "AGC"))
  expect_equal(CT$AAToCodon("T", TRUE), c("ACA", "ACC", "ACG"))
  expect_equal(CT$AAToCodon("V", TRUE), c("GTA", "GTC", "GTG"))
  expect_equal(CT$AAToCodon("W", TRUE), character(0))
  expect_equal(CT$AAToCodon("Y", TRUE), c("TAC"))
  
  expect_equal(CT$AAToCodon("a", TRUE), c("GCA", "GCC", "GCG")) #test that lower case works
  expect_equal(CT$AAToCodon("z", TRUE), character(0)) #check invalid
  
})