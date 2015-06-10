
initializeMCMCObject <- function(samples, thining=1, adaptive.width=100, est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
{
  mcmc <- new(MCMCAlgorithm, samples, thining, adaptive.width, est.expression, est.csp, est.hyper)
  return(mcmc)
}

runMCMC <- function(mcmc, genome, model, parameter)
{
  UseMethod("runMCMC", mcmc)
}
runMCMC.Rcpp_MCMCAlgorithm <- function(mcmc, genome, model, parameter)
{
  mcmc$run(genome, model, parameter)
}