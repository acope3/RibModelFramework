library(testthat)
library(ribModel)

context("Genome")

genome <- new(Genome)

test_that("check Index", {
    expect_true(genome$checkIndex(5,1,10))
    expect_true(genome$checkIndex(5,5,10))
    expect_true(genome$checkIndex(10,5,10))
    expect_false(genome$checkIndex(2,5,10))
    expect_false(genome$checkIndex(12,5,10))
})

test_that("add Gene And Get Genes", {

    # “Actual” genes
    gene <- new(Gene, "ATGCTCATTCTCACTGCTGCCTCGTAG", "001", "just a test")
    genome$addGene(gene, FALSE)
    gene2 <- new(Gene, "ATGCTCATTTAG", "002", "just a second test")
    genome$addGene(gene2, FALSE)
    expect_equal(genome$getGenes(FALSE), c(gene, gene2))

    # “Simulated” genes
    gene <- new(Gene, "ATGCTCATTCTCACTGCTGCCTCGTAG", "001", "just a test")
    genome$addGene(gene, TRUE)
    gene2 <- new(Gene, "ATGCTCATTTAG", "002", "just a second test")
    genome$addGene(gene2, TRUE)
    expect_equal(genome$getGenes(TRUE), c(gene, gene2))

})

test_that("clear", {
   genome$clear()
   expect_equal(genome$getGenes(FALSE), list())
   expect_equal(genome$getGenes(TRUE), list())
   expect_equal(genome$getNumGenesWithPhi(), 0)

})

test_that("get Genome Size", {
    expect_equal(genome$getGenomeSize(), 0)
    # “Actual” genes
    gene <- new(Gene, "ATGCTCATTCTCACTGCTGCCTCGTAG", "001", "just a test")
    genome$addGene(gene, FALSE)
    gene2 <- new(Gene, "ATGCTCATTTAG", "002", "just a second test")
    genome$addGene(gene2, FALSE)
    expect_equal(genome$getGenomeSize(), 2)
})

test_that("get Codon Counts Per Gene", {
    expect_equal(genome$getCodonCountsPerGene("ATG"), c(1,1))
    expect_equal(genome$getCodonCountsPerGene("CTC"), c(2,1))
    expect_equal(genome$getCodonCountsPerGene("ATT"), c(1,1))
    expect_equal(genome$getCodonCountsPerGene("ACT"), c(1,0))
    expect_equal(genome$getCodonCountsPerGene("GCT"), c(1,0))
    expect_equal(genome$getCodonCountsPerGene("GCC"), c(1,0))
    expect_equal(genome$getCodonCountsPerGene("TCG"), c(1,0))
    expect_equal(genome$getCodonCountsPerGene("TAG"), c(1,1))
    expect_equal(genome$getCodonCountsPerGene("CAT"), c(0,0))
})

test_that("read Fasta", {
    genome$clear()
    genome$readFasta("../../data/test.fasta", FALSE)
    compareGenome <- new(Genome)
    gene1 <- new(Gene, "ATGGCCACTATTGGGTCTTAG", "TEST001", "TEST001 Test Gene")
    gene2 <- new(Gene, "ATGACCGTAATTTTTTACTAG", "TEST002", "TEST002 Test Gene")
    gene3 <- new(Gene, "ATGGTCTACTTTCTGACATAG", "TEST003", "TEST003 Test Gene")
    compareGenome$addGene(gene1, FALSE)
    compareGenome$addGene(gene2, FALSE)
    compareGenome$addGene(gene3, FALSE)
    expect_equal(genome$getGenes(FALSE), compareGenome$getGenes(FALSE))
})

test_that("write Fasta", {
    genome$writeFasta("../../data/testWrite.fasta", FALSE)
    fromFile <- new(Genome)
    fromFile$readFasta("../../data/testWrite.fasta", FALSE)
    expect_equal(genome$getGenes(FALSE), fromFile$getGenes(FALSE))
})

test_that("read RFP File", {
    genome$clear()
    genome$readRFPFile("../../data/IKJ.csv")
    compareGenome <- new(Genome)
    gene1 <- new(Gene, "TAGATGGCCGCGGCGGGC", "YAL001C", "No description for RFP Model")
    gene1$setRFPObserved(1, 5)
    gene1$setRFPObserved(2, 4)
    gene1$setRFPObserved(13, 2)
    compss1 <- gene1$getSequenceSummary()
    gene2 <- new(Gene, "TTTTTTTAGCTTCTTATG", "YAL002W", "No description for RFP Model")
    gene2$setRFPObserved(11, 4)
    gene2$setRFPObserved(62, 3)
    gene2$setRFPObserved(26, 1)
    compss2 <- gene2$getSequenceSummary()
    gene3 <- new(Gene, "TAGATGATGATGATGATG", "YAL003W", "No description for RFP Model")
    gene3$setRFPObserved(29, 13)
    compss3 <- gene3$getSequenceSummary()
    compareGenome$addGene(gene1, FALSE)
    compareGenome$addGene(gene2, FALSE)
    compareGenome$addGene(gene3, FALSE)
    expect_equal(genome$getGenes(FALSE), compareGenome$getGenes(FALSE))
    expect_equal(genome$getGenomeSize(), 3)
    g1 <- genome$getGeneByIndex(1, FALSE)
    g2 <- genome$getGeneByIndex(2, FALSE)
    g3 <- genome$getGeneByIndex(3, FALSE)
    ss1 <- g1$getSequenceSummary()
    ss2 <- g2$getSequenceSummary()
    ss3 <- g3$getSequenceSummary()
    codonList <- codons()
    for (codon in codonList){
      expect_equal(compss1$getRFPObservedForCodon(codon),ss1$getRFPObservedForCodon(codon))
    }
    for (codon in codonList){
      expect_equal(compss2$getRFPObservedForCodon(codon),ss2$getRFPObservedForCodon(codon))
    }
    for (codon in codonList){
      expect_equal(compss3$getRFPObservedForCodon(codon),ss3$getRFPObservedForCodon(codon))
    }
})

test_that("write RFP File", {
  genome$writeRFPFile("../../data/testWriteRFP.csv", FALSE)
  fromFile <- new(Genome)
  fromFile$readRFPFile("../../data/testWriteRFP.csv")
  expect_equal(genome$getGenomeSize(), fromFile$getGenomeSize())
  g1 <- genome$getGeneByIndex(1, FALSE)
  g2 <- genome$getGeneByIndex(2, FALSE)
  g3 <- genome$getGeneByIndex(3, FALSE)
  ss1 <- g1$getSequenceSummary()
  ss2 <- g2$getSequenceSummary()
  ss3 <- g3$getSequenceSummary()
  ffg1 <- fromFile$getGeneByIndex(1, FALSE)
  ffg2 <- fromFile$getGeneByIndex(2, FALSE)
  ffg3 <- fromFile$getGeneByIndex(3, FALSE)
  compss1 <- ffg1$getSequenceSummary()
  compss2 <- ffg2$getSequenceSummary()
  compss3 <- ffg3$getSequenceSummary()
  codonList <- codons()
  for (codon in codonList){
    expect_equal(compss1$getRFPObservedForCodon(codon),ss1$getRFPObservedForCodon(codon))
    expect_equal(compss1$getCodonCountForCodon(codon), ss1$getCodonCountForCodon(codon))
  }
  for (codon in codonList){
    expect_equal(compss2$getRFPObservedForCodon(codon),ss2$getRFPObservedForCodon(codon))
    expect_equal(compss2$getCodonCountForCodon(codon), ss2$getCodonCountForCodon(codon))
  }
  for (codon in codonList){
    expect_equal(compss3$getRFPObservedForCodon(codon),ss3$getRFPObservedForCodon(codon))
    expect_equal(compss3$getCodonCountForCodon(codon), ss3$getCodonCountForCodon(codon))
  }
})

test_that("read Observed Phi Values", {
  genome$clear()
  genome$readFasta("../../data/test.fasta", FALSE)
  
  #Test for Valid values by Id
  genome$readObservedPhiValues("../../data/testReadObservedPhiValuesValid.csv", TRUE)
  gene1 <- genome$getGeneById("TEST001", FALSE)
  gene2 <- genome$getGeneById("TEST002", FALSE)
  gene3 <- genome$getGeneById("TEST003", FALSE)
  expect_equal(gene1$getObservedPhiValues(), c(0.4834983,1.4839493))
  expect_equal(gene2$getObservedPhiValues(), c(0.5388484,0.2222321))
  expect_equal(gene3$getObservedPhiValues(), c(0.4328382,2.3838239))
  
  #Test for valid values by index
  genome$clear()
  genome$readFasta("../../data/test.fasta", FALSE)
  genome$readObservedPhiValues("../../data/testReadObservedPhiValuesValid.csv", FALSE)
  gene1 <- genome$getGeneByIndex(1, FALSE)
  gene2 <- genome$getGeneByIndex(2, FALSE)
  gene3 <- genome$getGeneByIndex(3, FALSE)
  expect_equal(gene1$getObservedPhiValues(), c(0.4834983,1.4839493))
  expect_equal(gene2$getObservedPhiValues(), c(0.5388484,0.2222321))
  expect_equal(gene3$getObservedPhiValues(), c(0.4328382,2.3838239))
  
  
  #Test for when phi Values are missing
  genome$clear()
  genome$readFasta("../../data/test.fasta", FALSE)
  genome$readObservedPhiValues("../../data/testReadObservedPhiValuesMissing.csv", TRUE)
  gene1 <- genome$getGeneById("TEST001", FALSE)
  gene2 <- genome$getGeneById("TEST002", FALSE)
  gene3 <- genome$getGeneById("TEST003", FALSE)
  expect_equal(gene1$getObservedPhiValues(), numeric(0))
  expect_equal(gene2$getObservedPhiValues(), numeric(0))
  expect_equal(gene3$getObservedPhiValues(), numeric(0))
  
  
  genome$clear()
  genome$readFasta("../../data/test.fasta", FALSE)
  genome$readObservedPhiValues("../../data/testReadObservedPhiValuesMissing.csv", FALSE)
  gene1 <- genome$getGeneByIndex(1, FALSE)
  gene2 <- genome$getGeneByIndex(2, FALSE)
  gene3 <- genome$getGeneByIndex(3, FALSE)
  expect_equal(gene1$getObservedPhiValues(), numeric(0))
  expect_equal(gene2$getObservedPhiValues(), numeric(0))
  expect_equal(gene3$getObservedPhiValues(), numeric(0))
  
  
  #Test for valid phi files with "filler" values
  genome$clear()
  genome$readFasta("../../data/test.fasta", FALSE)
  genome$readObservedPhiValues("../../data/testReadObservedPhiValuesFilled.csv", TRUE)
  gene1 <- genome$getGeneById("TEST001", FALSE)
  gene2 <- genome$getGeneById("TEST002", FALSE)
  gene3 <- genome$getGeneById("TEST003", FALSE)
  expect_equal(gene1$getObservedPhiValues(), c(0.4834983,-1))
  expect_equal(gene2$getObservedPhiValues(), c(0.5388484,-1))
  expect_equal(gene3$getObservedPhiValues(), c(-1.323,2.3838239))
  
  
  genome$clear()
  genome$readFasta("../../data/test.fasta", FALSE)
  genome$readObservedPhiValues("../../data/testReadObservedPhiValuesFilled.csv", FALSE)
  gene1 <- genome$getGeneByIndex(1, FALSE)
  gene2 <- genome$getGeneByIndex(2, FALSE)
  gene3 <- genome$getGeneByIndex(3, FALSE)
  expect_equal(gene1$getObservedPhiValues(), c(0.4834983,-1))
  expect_equal(gene2$getObservedPhiValues(), c(0.5388484,-1))
  expect_equal(gene3$getObservedPhiValues(), c(-1.323,2.3838239))
  
  
  #Testing missing genes
  genome$clear()
  genome$readFasta("../../data/test.fasta", FALSE)
  genome$readObservedPhiValues("../../data/TestReadObservedPhiValuesGeneMissing.csv", TRUE)
  gene1 <- genome$getGeneByIndex(1, FALSE)
  gene2 <- genome$getGeneByIndex(2, FALSE)
  gene3 <- genome$getGeneByIndex(3, FALSE)
  expect_equal(gene1$getObservedPhiValues(), c(0.324394,4.3849298))
  expect_equal(gene2$getObservedPhiValues(), c(-1,-1))
  expect_equal(gene3$getObservedPhiValues(), c(3.348394,1.3943493))
  
  #Testing missing genes by index
  genome$clear()
  genome$readFasta("../../data/test.fasta", FALSE)
  genome$readObservedPhiValues("../../data/TestReadObservedPhiValuesMissingGeneByIndex.csv", FALSE)
  gene1 <- genome$getGeneByIndex(1, FALSE)
  gene2 <- genome$getGeneByIndex(2, FALSE)
  gene3 <- genome$getGeneByIndex(3, FALSE)
  expect_equal(gene1$getObservedPhiValues(), c(0.324394,4.3849298))
  expect_equal(gene2$getObservedPhiValues(), c(3.348394,1.3943493))
  expect_equal(gene3$getObservedPhiValues(), c(-1,-1))
  
})

test_that("get Num Genes With Phi", {
  expect_equal(genome$getNumGenesWithPhi(), 2)
})

genome$clear()
gene <- new(Gene, "ATGCTCATTCTCACTGCTGCCTCGTAG", "001", "just a test")
gene2 <- new(Gene, "ATGCTCATTTAG", "002", "just a second test")
genome$addGene(gene, FALSE)
genome$addGene(gene2, FALSE)
genome$addGene(gene, TRUE)
genome$addGene(gene2, TRUE)

test_that("get Gene By Index", {
    expect_equal(genome$getGeneByIndex(1, FALSE), gene)
    expect_equal(genome$getGeneByIndex(2, FALSE), gene2)
    expect_equal(genome$getGeneByIndex(1, TRUE), gene)
    expect_equal(genome$getGeneByIndex(2, TRUE), gene2)
})

test_that("get Gene By Id", {
    expect_equal(genome$getGeneById("001", FALSE), gene)
    expect_equal(genome$getGeneById("002", FALSE), gene2)
    expect_equal(genome$getGeneById("001", TRUE), gene)
    expect_equal(genome$getGeneById("002", TRUE), gene2)
})

test_that("get Genome For Gene Indicies", {
    tmp <- new(Genome)
    blank <- new(Genome)
    secondHalf <- new(Genome)
    secondHalf$addGene(gene2, FALSE)
    tmp <- genome$getGenomeForGeneIndicies(c(1,2), FALSE)
    expect_equal(tmp$getGenomeSize(), genome$getGenomeSize())
    expect_equal(tmp$getGeneByIndex(1, FALSE), genome$getGeneByIndex(1, FALSE))
    expect_equal(tmp$getGeneByIndex(2, FALSE), genome$getGeneByIndex(2, FALSE))
    tmp <- genome$getGenomeForGeneIndicies(c(3), FALSE)
    expect_equal(tmp$getGenomeSize(), blank$getGenomeSize())
    tmp <- genome$getGenomeForGeneIndicies(c(2), FALSE)
    expect_equal(tmp$getGenomeSize(), secondHalf$getGenomeSize())
    expect_equal(tmp$getGeneByIndex(1, FALSE), secondHalf$getGeneByIndex(1, FALSE))
    
    #simulated
    secondHalf$clear()
    secondHalf$addGene(gene2, FALSE)
    tmp <- genome$getGenomeForGeneIndicies(c(1,2), TRUE)
    expect_equal(tmp$getGenomeSize(), genome$getGenomeSize())
    expect_equal(tmp$getGeneByIndex(1, FALSE), genome$getGeneByIndex(1, FALSE))
    expect_equal(tmp$getGeneByIndex(2, FALSE), genome$getGeneByIndex(2, FALSE))
    tmp <- genome$getGenomeForGeneIndicies(c(3), TRUE)
    expect_equal(tmp$getGenomeSize(), blank$getGenomeSize())
    tmp <- genome$getGenomeForGeneIndicies(c(2), TRUE)
    expect_equal(tmp$getGenomeSize(), secondHalf$getGenomeSize())
    expect_equal(tmp$getGeneByIndex(1, FALSE), secondHalf$getGeneByIndex(1, FALSE))
})
