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
    expect_equal(genome, fromFile)
})

test_that("read RFP File", {
    genome$clear()
    genome$readRFPFile("../../data/testRFPFile.csv")
    compareGenome <- new(Genome)
    gene1 <- new(Gene, "TAGATGGCCGCGGCGGGC", "YAL001C", "YAL001C No description for RFP Model")
    gene2 <- new(Gene, "TTTTTTTAGCTTCTTATG", "YAL002W", "YAL002W No description for RFP Model")
    gene3 <- new(Gene, "TAGATGATGATGATGATG", "YAL003W", "YAL003W No description for RFP Model")
    compareGenome$addGene(gene1, FALSE)
    compareGenome$addGene(gene2, FALSE)
    compareGenome$addGene(gene3, FALSE)
    expect_equal(genome$getGenes(FALSE), compareGenome$getGenes(FALSE))
})

genome$clear()
gene <- new(Gene, "ATGCTCATTCTCACTGCTGCCTCGTAG", "001", "just a test")
gene2 <- new(Gene, "ATGCTCATTTAG", "002", "just a second test")
genome$addGene(gene, FALSE)
genome$addGene(gene2, FALSE)


test_that("get Gene By Index", {
    expect_equal(genome$getGeneByIndex(1), gene)
    expect_equal(genome$getGeneByIndex(2), gene2)
})

test_that("get Gene By Id", {
    expect_equal(genome$getGeneById("001"), gene)
    expect_equal(genome$getGeneById("002"), gene2)
})

test_that("get Genome For Gene Indicies", {
    blank <- new(Genome)
    secondHalf <- new(Genome)
    secondHalf$addGene(gene2, FALSE)
    expect_equal(genome$getGenomeForGeneIndicies(c(1,2)), genome)
    expect_equal(genome$getGenomeForGeneIndicies(c(3)), blank)
    expect_equal(genome$getGenomeForGeneIndicies(c(2)), secondHalf)
})
