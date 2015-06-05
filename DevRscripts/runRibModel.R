library(ribModel)
genome <- new(Genome)
genome$readFasta("../ribModel/data/Skluyveri_ChrA_ChrB_andCleft.fasta", F)


sphi_init <- 2;
numMixtures <- 2;
mixDef <- "allUnique";
geneAssignment <- c(rep(0,448), rep(0,513), rep(1,457))
parameter <- new(ROCParameter, genome$getGenomeSize(), sphi_init, numMixtures, geneAssignment, T, mixDef)
parameter$initializeExpressionByGenome(genome, sphi_init)
files <- c("../ribModel/data/Skluyveri_CSP_ChrA.csv", "../ribModel/data/Skluyveri_CSP_ChrCleft.csv")
parameter$initMutationSelectionCategories(files, 2, 0)
parameter$initMutationSelectionCategories(files, 2, 1)



samples <- 1000
thining <- 10
adaptiveWidth <- 10
useSamples <- 50

model <- new(ROCModel)
mcmc <- new(MCMCAlgorithm, samples, thining, adaptiveWidth, T, T, T)
mcmc$run(genome, model, parameter);

plot(mcmc)

plot(parameter, what = "MixtureProbability")
plot(parameter, what = "SPhi")
plot(parameter, what = "ExpectedPhi")
plot(parameter, what = "Expression", geneIndex = 905)

plot(parameter$getExpressionTraceByCategoryForGene( 2, 905 ), type="l")


plot(parameter, what = "Mutation", category = 1)
plot(parameter, what = "Selection", category = 1)

plot(model, genome, parameter, samples = 100, category = 1)


