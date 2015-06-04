library(ribModel)
genome <- new(Genome)
genome$readFasta("../ribModel/data/Skluyveri.fasta")


sphi_init <- 2;
numMixtures <- 2;
mixDef <- "allUnique";
geneAssignment <- c(rep(0,448), rep(1,457))
parameter <- new(ROCParameter, 905, sphi_init, numMixtures, geneAssignment, T, mixDef)
parameter$initializeExpressionByGenome(genome, sphi_init)
# I do not initialize CSP right now


samples <- 100
thining <- 10
adaptiveWidth <- 100
useSamples <- 50

model <- new(ROCModel)
mcmc <- new(MCMCAlgorithm, samples, thining, adaptiveWidth, T, T, T)
mcmc$run(genome, model, parameter);

plot(NULL, NULL, xlim = c(0,samples), ylim=c(0,1), xlab = "samples", ylab="Mixture Probability")
lines(parameter$getCategoryProbabilitiesTrace(0), col="black")
lines(parameter$getCategoryProbabilitiesTrace(1), col="red")

