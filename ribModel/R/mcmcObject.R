
initializeMCMCObject <- function(samples, thining=1, adaptive.width=100, est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
{
  mcmc <- new(MCMCAlgorithm, samples, thining, adaptive.width, est.expression, est.csp, est.hyper)
  return(mcmc)
}

runMCMC <- function(mcmc, genome, model)
{
  UseMethod("runMCMC", mcmc)
}
runMCMC.Rcpp_MCMCAlgorithm <- function(mcmc, genome, model)
{
  mcmc$run(genome, model)
}
setRestartSettings <- function(mcmc, filename, samples, write.multiple=TRUE)
{
  UseMethod("setRestartSettings", mcmc)
}
setRestartSettings.Rcpp_MCMCAlgorithm <- function(mcmc, filename, samples, write.multiple=TRUE)
{
  mcmc$setRestartFileSettings(filename, samples, write.multiple)
}
