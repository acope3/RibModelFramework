
initializeMCMCObject <- function(samples, thining=1, adaptive.width=100, est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
{
  mcmc <- new(MCMCAlgorithm, samples, thining, adaptive.width, est.expression, est.csp, est.hyper)
  return(mcmc)
}

runMCMC <- function(mcmc, genome, model, ncores = 1, divergence.iteration = 0)
{
  UseMethod("runMCMC", mcmc)
}
runMCMC.Rcpp_MCMCAlgorithm <- function(mcmc, genome, model, ncores = 1, divergence.iteration = 0)
{
  mcmc$run(genome, model, ncores, divergence.iteration)
}
setRestartSettings <- function(mcmc, filename, samples, write.multiple=TRUE)
{
  UseMethod("setRestartSettings", mcmc)
}
setRestartSettings.Rcpp_MCMCAlgorithm <- function(mcmc, filename, samples, write.multiple=TRUE)
{
  mcmc$setRestartFileSettings(filename, samples, write.multiple)
}

convergence.test <- function(mcmc, n.samples = 10, frac1 = 0.1, frac2 = 0.5, plot = FALSE)
{
  UseMethod("convergence.test", mcmc)
}
convergence.test.Rcpp_MCMCAlgorithm <- function(mcmc, n.samples = 10, frac1 = 0.1, frac2 = 0.5, plot = FALSE)
{
  # TODO: extend to work with multiple chains once we have that capability.
  
  loglik.trace <- mcmc$getLogLikelihoodTrace()
  trace.length <- length(loglik.trace)
  start <- max(1, trace.length - n.samples)
  
  # the start and end parameter do NOT work, using subsetting to achieve goal
  mcmcobj <- coda::mcmc(data=loglik.trace[start:trace.length])
  diag <- coda::geweke.diag(mcmcobj, frac1=frac1, frac2=frac2)
  if(plot){ 
    coda::geweke.plot(mcmcobj, frac1=frac1, frac2=frac2)
  }else{
    return(diag)
  }
}



