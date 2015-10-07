## run variables
with.phi <- FALSE
genome <- "../data/realGenomes/Skluyveri.fasta"
#genome <- "../data/twoMixtures/simulatedAllUniqueR.fasta"
#genome <- "../data/genome_2000.fasta"
#genome <- "../data/twoMixtures/simulatedKluyveri_full.fasta"
expression.file <- "../data/simulatedAllUniqueR_phi.csv"
observed.phi.file <- "../data/realGenomes/Skluyveri_phi.csv"

observed.mutation <- "c(\"../data/realGenomes/Skluyveri_mutation_ChrA.csv\", \"../data/realGenomes/Skluyveri_mutation_ChrCleft.csv\")"
observed.selection <- "c(\"../data/realGenomes/Skluyveri_selection_ChrA.csv\", \"../data/realGenomes/Skluyveri_selection_ChrCleft.csv\")"

sphi.init <- "c(1,1)"
numMixtures <- 2
mixdef <- "allUnique"
#geneAssignment <- "c(rep(1,4860), rep(2,457))" # simulated Skluyveri
geneAssignment <- "c(rep(1,961), rep(2,457), rep(1, 3903))" # S.kluyveri full genome
#geneAssignment <- "c(rep(1,448), rep(1,513), rep(2,457))"
#geneAssignment <- "c(rep(1,500), rep(2,500))"
#geneAssignment <- "rep(1,2000)"
samples <- 100
thining <- 10
adaptiveWidth <- 10
divergence.iteration <- 0
model <- "ROC"

## queue variables
run.name <- "_kluyveri_full"
name <- "clandere_kluyveri_"
n.cores <- 4
queue = "long*"
number.of.runs <- 1

for (index in 1:number.of.runs) {
    cat("rm(list=ls()) \n",
        "library(ribModel) \n",
        "with.phi <- ", with.phi, "\n\n",
        "if (with.phi) { \n",
        "genome <- initializeGenomeObject(file =\"", genome, "\", expression.file = \"", expression.file, "\") \n",
        "} else { \n",
        "genome <- initializeGenomeObject(file = \"", genome, "\") \n",
        "} \n\n",
        "sphi_init <- ", sphi.init , "\n",
        "numMixtures <- ", numMixtures, "\n",
        "mixDef <- \"", mixdef, "\"\n",
        "geneAssignment <- ", geneAssignment, "\n", 
        "parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef) \n\n",
        "samples <- ", samples, "\n", 
        "thining <- ", thining, "\n",
        "adaptiveWidth <- ", adaptiveWidth, "\n",
        "divergence.iteration <- ", divergence.iteration, "\n",
        "mcmc <- initializeMCMCObject(samples, thining, adaptive.width=adaptiveWidth, est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE) \n\n", 
        "setRestartSettings(mcmc, \"", paste(index, run.name, ".rst", sep=""), "\", adaptiveWidth*200, TRUE) \n\n",
        "model <- initializeModelObject(parameter, \"", model, "\", with.phi = with.phi) \n\n",
        "start <- Sys.time() \n",
        "system.time(runMCMC(mcmc, genome, model, ", n.cores, ", divergence.iteration)) \n",
        "end <- Sys.time() \n",
        "end - start \n\n\n",
        "pdf(\"", paste(index, run.name, "_global.pdf", sep=""), "\", width = 11, height = 12) \n",
        "plot(mcmc) \n",
        "loglik.trace <- mcmc$getLogLikelihoodTrace() \n",
        "acf(loglik.trace) \n",
        "convergence.test(mcmc, n.samples = 500, plot=T) \n\n",
        "trace <- parameter$getTraceObject() \n",
        "plot(trace, what = \"MixtureProbability\") \n",
        "plot(trace, what = \"Sphi\") \n",
        "plot(trace, what = \"Mphi\") \n",
        "if (with.phi) { \n",
          "plot(trace, what = \"Aphi\") \n",
          "plot(trace, what = \"Sepsilon\") \n",
        "} \n",
        "plot(trace, what = \"ExpectedPhi\") \n\n",
        "mixtureAssignment <- unlist(lapply(1:genome$getGenomeSize(),  function(geneIndex){parameter$getEstimatedMixtureAssignmentForGene(samples, geneIndex)})) \n",
        "expressionValues <- unlist(lapply(1:genome$getGenomeSize(), function(geneIndex){ \n",
        "  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixtureAssignment[geneIndex]) \n",
        "  parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples, geneIndex, expressionCategory) \n",
        "})) \n",
        "expressionValues <- log10(expressionValues) \n",
        "obs.phi <- log10(read.table(\"", observed.phi.file, "\", sep=\",\", header=T)[, 2]) \n",
        "plot(NULL, NULL, xlim=range(obs.phi) + c(-0.1, 0.1), ylim=range(expressionValues, na.rm = T) + c(-0.1, 0.1), \n",
             "main = \"Synthesis Rate\", xlab = \"True values\", ylab = \"Estimated values\") \n",
        "for(k in 1:numMixtures){ \n",
        "  upper.panel.plot(obs.phi[mixtureAssignment == k], expressionValues[mixtureAssignment == k], col=ribModel:::.mixtureColors[k]) \n",
        "} \n",
        "legend(\"topleft\", legend = paste(\"Mixture Element\", 1:numMixtures), \n",
        "       col = ribModel:::.mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = \"n\") \n\n",
        "plot(parameter, what = \"Mutation\") \n",
        "plot(parameter, what = \"Selection\") \n",
        "dev.off() \n\n",
        "observed.mutation <- ", observed.mutation, "\n",
        "observed.selection <- ", observed.selection, "\n",        
        "for(k in 1:numMixtures){ \n",
        "   pdf(paste(\"", paste(index, run.name, sep=""), "_mixture_\", k, \".pdf\", sep=\"\"), width = 11, height = 12) \n",
        "   mixture <- k \n", 
        "   plot(trace, what = \"Mutation\", mixture = mixture) \n",
        "   plot(trace, what = \"Selection\", mixture = mixture) \n",
        "   plot(model, genome, parameter, samples = samples*0.1, mixture = mixture, main = \"Codon Usage Plot\") \n",
        "   names.aa <- aminoAcids() \n",
        "   selection.ci <- c() \n",
        "   mutation.ci <- c() \n",
        "   selection <- c() \n",
        "   mutation <- c() \n",
        "   codon.storage <- c() \n",
        "   csp.m <- read.table(observed.mutation[k], sep=\",\", header=T) \n",
        "   csp.e <- read.table(observed.selection[k], sep=\",\", header=T) \n",
        "   csp <- rbind(csp.m,csp.e) \n",
        "   idx.eta <- 41:80 \n",
        "   idx.mu <- 1:40 \n",
        "   for(aa in names.aa) \n",
        "   { \n",
        "     if(aa == \"M\" || aa == \"W\" || aa == \"X\") next \n",
        "     codons <- AAToCodon(aa, T) \n",
        "     codon.storage <- c(codon.storage,codons) \n",
        "     for(i in 1:length(codons)) \n",
        "     { \n",
        "       selection <- c(selection, parameter$getSelectionPosteriorMeanForCodon(mixture, samples*0.1, codons[i])) \n",
        "       selection.ci <- c(selection.ci, parameter$getSelectionVarianceForCodon(mixture, samples*0.1, codons[i], TRUE)) \n",
        "       mutation <- c(mutation, parameter$getMutationPosteriorMeanForCodon(mixture, samples*0.1, codons[i])) \n",
        "       mutation.ci <- c(mutation.ci, parameter$getMutationVarianceForCodon(mixture, samples*0.1, codons[i], TRUE)) \n",
        "     } \n",
        "   } \n",
        "   plot(NULL, NULL, xlim=range(csp[idx.mu, 3], na.rm = T), ylim=range(mutation), \n",
        "        main = \"Mutation\", xlab = \"True values\", ylab = \"Estimated values\") \n",
        "   upper.panel.plot(x = csp[idx.mu, 3], y = mutation, sd.y = 1.96*sqrt(mutation.ci)) \n",
        "   plot(NULL, NULL, xlim=range(csp[idx.eta, 3], na.rm = T), ylim=range(selection), \n",
        "        main = \"Selection\", xlab = \"True values\", ylab = \"Estimated values\") \n",
        "   upper.panel.plot(x = csp[idx.eta, 3], y = selection, sd.y = 1.96*sqrt(selection.ci)) \n",
        "   dev.off() \n",
        "} \n",
        file=paste(index, run.name, ".R", sep=""), sep="")
    
    qsubCommands <- paste("#$ -N ", name, index, "\n", sep="")
    qsubCommands <- paste(qsubCommands, "#$ -pe threads ", n.cores, "\n", sep="")
    qsubCommands <- paste(qsubCommands, "#$ -q ", queue, "\n", sep="")
    qsubCommands <- paste(qsubCommands, "#$ -cwd \n", sep="")
    qsubCommands <- paste(qsubCommands, "#$ -M clandere@vols.utk.edu \n", sep="")
    qsubCommands <- paste(qsubCommands, "module load R/3.1.0 \n", sep="")
    qsubCommands <- paste(qsubCommands, "R CMD BATCH ", index, run.name, ".R", sep="")
    
    qsubFile <- paste(index, run.name, ".sge", sep="")
    cat(qsubCommands, file=qsubFile, append=FALSE)
    #system(paste("chmod u+x ", qsubFile, sep=""))
    #system(paste("/opt/grid/bin/lx26-amd64/qsub ", qsubFile, sep=""))
}

