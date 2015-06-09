library(ribModel)
genome <- new(Genome)
genome$readFasta("../ribModel/data/Skluyveri_ChrA_andCleft.fasta", F)


sphi_init <- 2;
numMixtures <- 2;
mixDef <- "allUnique";
geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457))
parameter <- new(ROCParameter, genome$getGenomeSize(), sphi_init, numMixtures, geneAssignment, T, mixDef)
parameter$initializeExpressionByGenome(genome, sphi_init)
files <- c("../ribModel/data/Skluyveri_CSP_ChrA.csv", "../ribModel/data/Skluyveri_CSP_ChrCleft.csv")
parameter$initMutationSelectionCategories(files, 2, "Mutation")
parameter$initMutationSelectionCategories(files, 2, "Selection")

#check to see if the initCovarianceMatrix is working for ROCParameter
A = matrix(c(2, 4, 3, 1), nrow = 2, ncol = 2, byrow = TRUE)
 parameter$initCovarianceMatrix( A, "A" )
cm <- parameter$getCovarianceMatrixForAA("a")
cm$printCovarianceMatrix()
A
