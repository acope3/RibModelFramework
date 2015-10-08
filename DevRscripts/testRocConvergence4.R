rm(list=ls())
library(ribModel)

with.phi <- FALSE

if (with.phi) {
  genome <- initializeGenomeObject(file = "../ribModel/data/simulatedAllUniqueR.fasta", expression.file = "../ribModel/data/simulatedAllUniqueR_phi.csv")
} else {
  #genome <- initializeGenomeObject(file = "../ribModel/data/simulatedAllUniqueR.fasta")
  genome <- initializeGenomeObject(file = "../ribModel/data/simulatedOneMix.fasta")
}

#initialize parameter object
sphi_init <- 1
numMixtures <- 1
mixDef <- "allUnique"
#geneAssignment <- c(rep(1,500), rep(2,500))
geneAssignment <- rep(1,1000)
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)

# initialize MCMC object
samples <- 10000
thining <- 20
adaptiveWidth <- 10
divergence.iteration <- 50
mcmc <- initializeMCMCObject(samples, thining, adaptive.width=adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)

# get model object
model <- initializeModelObject(parameter, "ROC", with.phi = with.phi)

#run mcmc on genome with parameter using model
start <- Sys.time()
system.time(
    runMCMC(mcmc, genome, model, 4, divergence.iteration)
)
end <- Sys.time()
end - start

#collect traces
loglik.trace <- mcmc$getLogLikelihoodTrace()
trace <- parameter$getTraceObject()
expected.phi.trace <- trace$getExpectedPhiTrace()
sphi.trace <- trace$getSphiTrace()

initialExpressionValues <- unlist(lapply(1:genome$getGenomeSize(), function(geneIndex){
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(geneAssignment[geneIndex])
  trace$getSynthesisRateTraceForGene(geneIndex)[1]
}))

estsimatedExpressionValues <- unlist(lapply(1:genome$getGenomeSize(), function(geneIndex){
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(geneAssignment[geneIndex])
  parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples, geneIndex, expressionCategory)
}))

colour <- c("black", "chartreuse4", "red", "blue", "blue4", "blueviolet", "darkgoldenrod2", "darkgreen", "darkorchid1", "deeppink2",
            "khaki4", "midnightblue", "lightsteelblue", "ivory4", "gray47", "orangered3", "slateblue4", "yellow3", "tomato3", "turquoise3",
            "plum", "orangered", "red4", "navy", "gold", "darkred", "darkmagenta", "burlywood4", "mediumseagreen", "cornflowerblue",
            "cyan2", "darkcyan", "darkolivegreen3", "darkturquoise", "hotpink", "lightpink4", "mediumaquamarine", "springgreen2", "chartreuse", "azure4")


run.name <- "single_mixture_convergence_test_4"

pdf(paste(run.name, ".pdf", sep=""))
plot(mcmc)
acf(loglik.trace)
plot(trace, what = "MixtureProbability")
plot(trace, what = "Sphi")

plot(trace, what = "ExpectedPhi")
mixture <- 1
plot(trace, what = "Mutation", mixture = mixture)
plot(trace, what = "Selection", mixture = mixture)


mixtureAssignment <- unlist(lapply(1:genome$getGenomeSize(),  function(geneIndex){parameter$getEstimatedMixtureAssignmentForGene(samples, geneIndex)}))
expressionValues <- unlist(lapply(1:genome$getGenomeSize(), function(geneIndex){
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixtureAssignment[geneIndex])
  parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples, geneIndex, expressionCategory)
}))
expressionValues <- log10(expressionValues)
#obs.phi <- log10(read.table("../ribModel/data/simulatedAllUniqueR_phi.csv", sep=",", header=T)[, 2])
obs.phi <- log10(read.table("../ribModel/data/simulatedOneMix_phi.csv", sep=",", header=T)[, 2])
plot(NULL, NULL, xlim=range(obs.phi) + c(-0.1, 0.1), ylim=range(expressionValues, na.rm = T) + c(-0.1, 0.1), 
     main = "Synthesis Rate", xlab = "true values", ylab = "estimated values")
upper.panel.plot(obs.phi[mixtureAssignment == 1], expressionValues[mixtureAssignment == 1], col="black")
#upper.panel.plot(obs.phi[mixtureAssignment == 2], expressionValues[mixtureAssignment == 2], col="red")
legend("topleft", legend = paste("Mixture Element", 1:numMixtures), 
       col = ribModel:::.mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")



#true.sel <- read.csv(file="../ribModel/data/simulated_selection0.csv", header=T)
true.sel <- read.csv(file="../ribModel/data/simulatedOneMix_selection.csv", header=T)
names.aa <- aminoAcids()
start.value <- c()
end.value <- c()
post.estm <- c()
post.var <- c()
true.value <- c()
for(aa in names.aa)
{
  if(aa == "M" || aa == "W" || aa == "X") next
  codons <- AAToCodon(aa, T)
  for(i in 1:length(codons))
  {
      cur.trace <- trace$getSelectionParameterTraceByMixtureElementForCodon(1, codons[i])
      start.value <- c(start.value, cur.trace[1])
      end.value <- c(end.value, cur.trace[length(cur.trace)])
      true.value <- c(true.value, true.sel[as.character(true.sel[, 2]) == codons[i], 3])
      post.estm <- c(post.estm, parameter$getSelectionPosteriorMeanForCodon(1, length(cur.trace) * 0.2, codons[i]) )
      post.var <- c(post.var, 1.96*sqrt(parameter$getSelectionVarianceForCodon(1, length(cur.trace) * 0.2, codons[i], TRUE)) )
  }
}
up <- post.estm + post.var
low <- post.estm - post.var
limits <- range(c(start.value, post.estm, true.value, low, up))
plot(NULL, NULL, xlim=limits, ylim=limits, ylab=expression(Delta~eta["estm."]), xlab=expression(Delta~eta["true"]), axes = FALSE, main = "Translational Inefficiency")
upper.panel.plot(true.value, post.estm)
box()
abline(0,1)
arrows(x0=true.value, y0=start.value, x1=true.value, y1=post.estm, col=colour)
points(true.value, start.value, col=colour, pch=20)
points(true.value, post.estm, col=colour, pch=1)
epsilon <- 0.02
segments(true.value, low, true.value, up, col=colour)
segments(true.value-epsilon, up ,true.value+epsilon, up, col=colour)
segments(true.value-epsilon, low ,true.value+epsilon, low, col=colour)
axis(1)
axis(2)
legend("topleft", legend = c("I.C.", "Posterior"), col=c("black", "black"), pch=c(20, 1), bty = "n")
start.sel.values <- start.value
cat("Percent quantiles overlap with true value: ", sum(true.value < up & true.value > low) / 40, "\n")
cat("Geweke Score: ", convergence.test(mcmc, n.samples = 500, plot=F)$z, ", using 500 samples\n")
cat("Geweke Score: ", convergence.test(mcmc, n.samples = 5000, plot=F)$z, ", using 5,000 samples\n")
cat("Geweke Score: ", convergence.test(mcmc, n.samples = 50000, plot=F)$z, ", using 50,000 samples\n")

#true.sel <- read.csv(file="../ribModel/data/simulated_mutation0.csv", header=T)
true.sel <- read.csv(file="../ribModel/data/simulatedOneMix_mutation.csv", header=T)
names.aa <- aminoAcids()
start.value <- c()
end.value <- c()
post.estm <- c()
post.var <- c()
true.value <- c()
for(aa in names.aa)
{
  if(aa == "M" || aa == "W" || aa == "X") next
  codons <- AAToCodon(aa, T)
  for(i in 1:length(codons))
  {
    cur.trace <- trace$getMutationParameterTraceByMixtureElementForCodon(1, codons[i])
    start.value <- c(start.value, cur.trace[1])
    end.value <- c(end.value, cur.trace[length(cur.trace)])
    true.value <- c(true.value, true.sel[as.character(true.sel[, 2]) == codons[i], 3])
    post.estm <- c(post.estm, parameter$getMutationPosteriorMeanForCodon(1, length(cur.trace) * 0.2, codons[i]) )
    post.var <- c(post.var, 1.96*sqrt(parameter$getMutationVarianceForCodon(1, length(cur.trace) * 0.2, codons[i], TRUE)) )
  }
}
up <- post.estm + post.var
low <- post.estm - post.var
limits <- range(c(start.value, post.estm, true.value, low, up))
plot(NULL, NULL, xlim=limits, ylim=limits, ylab=expression(Delta~"M"["estm."]), xlab=expression(Delta~"M"["true"]), axes = FALSE, main = "Mutation Bias")
upper.panel.plot(true.value, post.estm)
box()
abline(0,1)
arrows(x0=true.value, y0=start.value, x1=true.value, y1=post.estm, col=colour)
points(true.value, start.value, col=colour, pch=20)
points(true.value, post.estm, col=colour, pch=1)
epsilon <- 0.02
segments(true.value, low, true.value, up, col=colour)
segments(true.value-epsilon, up ,true.value+epsilon, up, col=colour)
segments(true.value-epsilon, low ,true.value+epsilon, low, col=colour)
axis(1)
axis(2)
legend("topleft", legend = c("I.C.", "Posterior"), col=c("black", "black"), pch=c(20, 1), bty = "n")
start.mut.values <- start.value
cat("Percent quantiles overlap with true value: ", sum(true.value < up & true.value > low) / 40, "\n")
cat("Geweke Score: ", convergence.test(mcmc, n.samples = 500, plot=F)$z, ", using 500 samples\n")
cat("Geweke Score: ", convergence.test(mcmc, n.samples = 5000, plot=F)$z, ", using 5,000 samples\n")
cat("Geweke Score: ", convergence.test(mcmc, n.samples = 50000, plot=F)$z, ", using 50,000 samples\n")

dev.off()


save(list=c("start.mut.values", "start.sel.values", "initialExpressionValues", "estsimatedExpressionValues", 
            "loglik.trace", "expected.phi.trace", "sphi.trace"), file = paste(run.name, ".rda", sep=""))

