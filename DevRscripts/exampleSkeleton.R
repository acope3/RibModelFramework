## run variables
with.phi <- FALSE
genome <- "../ribModel/data/realGenomes/Skluyveri.fasta"
#genome <- "../ribModel/data/twoMixtures/simulatedAllUniqueR.fasta"
#genome <- "../ribModel/data/genome_2000.fasta"
#genome <- "../ribModel/data/twoMixtures/simulatedKluyveri_full.fasta"
expression.file <- "../ribModel/data/simulatedAllUniqueR_phi.csv"
sphi.init <- c(1,1)
numMixtures <- 2
mixdef <- "allUnique"
#geneAssignment <- c(rep(1,4860), rep(2,457)) # simulated Skluyveri
geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457), rep(1, 3903)) # S.kluyveri full genome
#geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457))
#geneAssignment <- rep(1,328) # human brain
#geneAssignment <- c(rep(1,500), rep(2,500))
#geneAssignment <- rep(1,2000)
samples <- 1000
thining <- 10
adaptiveWidth <- 10
divergence.iteration <- 0
model <- "ROC"

## queue variables
run.name <- "_kluyveri_full"
name <- "clandere_kluyveri_"
n.cores <- 8
queue = "long*"
for (index in ) {
    cat("rm(list=ls()) \n",
        "library(ribModel) \n",
        "with.phi <- ", with.phi, "\n",
        "if (with.phi) { \n",
        "genome <- initializeGenomeObject(file = ", genome, ", expression.file = "../ribModel/data/simulatedAllUniqueR_phi.csv") \n",
        "} else { \n",
        "genome <- initializeGenomeObject(file = ", genome, ") \n",
        "} \n",
        "sphi_init <- ", sphi.init , "\n",
        "numMixtures <- ", numMixtures, "\n",
        "mixDef <- ", mixdef, "\n",
        "geneAssignment <- " geneAssignment, "\n", 
        "parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)"
        "samples <- ", samples, "\n", 
        "thining <- ", thining, "\n",
        "adaptiveWidth <- ", adaptiveWidth, "\n",
        "divergence.iteration <- ", divergence.iteration, "\n",
        "mcmc <- initializeMCMCObject(samples, thining, adaptive.width=adaptiveWidth, est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE) \n", 
        "model <- initializeModelObject(parameter, ", model, ", with.phi = with.phi) \n",
        "setRestartSettings(mcmc, ", paste(index, run.name, ".rst", sep=""), ", adaptiveWidth*200, TRUE), \n",
        "start <- Sys.time() \n",
        "system.time(runMCMC(mcmc, genome, model, ", n.cores,", divergence.iteration)) \n",
        "end <- Sys.time() \n",
        "end - start \n",
        "pdf("paste(index, run.name, "_global.pdf", sep="")", width = 11, height = 12) \n",
        "plot(mcmc) \n",
        "loglik.trace <- mcmc$getLogLikelihoodTrace() \n",
        "acf(loglik.trace) \n",
        "convergence.test(mcmc, n.samples = 500, plot=T) \n",
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
    system(paste("chmod u+x ", pbsFile, sep=""))
    system(paste("/opt/grid/bin/lx26-amd64/qsub ", qsubFile, sep=""))
}

