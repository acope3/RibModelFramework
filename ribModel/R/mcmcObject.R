
initializeMCMCObject <- function(samples, thining=1, adaptive.width=100, est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
{
  mcmc <- new(MCMCAlgorithm, samples, thining, adaptive.width, est.expression, est.csp, est.hyper)
  return(mcmc)
}

runMCMC <- function(mcmc, genome, model, ncores = 1)
{
  UseMethod("runMCMC", mcmc)
}
runMCMC.Rcpp_MCMCAlgorithm <- function(mcmc, genome, model, ncores = 1)
{
  mcmc$run(genome, model, ncores)
}
setRestartSettings <- function(mcmc, filename, samples, write.multiple=TRUE)
{
  UseMethod("setRestartSettings", mcmc)
}
setRestartSettings.Rcpp_MCMCAlgorithm <- function(mcmc, filename, samples, write.multiple=TRUE)
{
  mcmc$setRestartFileSettings(filename, samples, write.multiple)
}
