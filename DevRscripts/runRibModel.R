library(ribModel)
genome <- new(Genome)
genome$readFasta("../ribModel/data/Skluyveri_ChrA_ChrB_andCleft.fasta", F)

sphi_init <- 2;
numMixtures <- 2;
mixDef <- "allUnique";
geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457))

parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)

samples <- 5000
thining <- 10
adaptiveWidth <- 10
useSamples <- 10

model <- new(ROCModel)
mcmc <- new(MCMCAlgorithm, samples, thining, adaptiveWidth, TRUE, TRUE, TRUE)
mcmc$run(genome, model, parameter);

plot(mcmc)

plot(parameter, what = "MixtureProbability")
plot(parameter, what = "SPhi")
plot(parameter, what = "ExpectedPhi")
plot(parameter, what = "Expression", geneIndex = 905)
plot(parameter, what = "Mutation", mixture = 2)
plot(parameter, what = "Selection", mixture = 2)

plot(model, genome, parameter, samples = samples*0.1, mixture = 1, main = "S. kluyveri Chr (A,B,Cleft) Codon Usage Plot")

