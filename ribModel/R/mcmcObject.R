
initializeMCMCObject <- function(samples, thining=1, adaptive.width=100, est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE, est.mix=TRUE)
{
  mcmc <- new(MCMCAlgorithm, samples, thining, adaptive.width, est.expression, est.csp, est.hyper)
  mcmc$setEstimateMixtureAssignment(est.mix)
  return(mcmc)
}

runMCMC <- function(mcmc, genome, model, ncores = 1, divergence.iteration = 0)
{
  UseMethod("runMCMC", mcmc)
}

# NOT EXPOSED
runMCMC.Rcpp_MCMCAlgorithm <- function(mcmc, genome, model, ncores = 1, divergence.iteration = 0)
{
  mcmc$run(genome, model, ncores, divergence.iteration)
}
setRestartSettings <- function(mcmc, filename, samples, write.multiple=TRUE)
{
  UseMethod("setRestartSettings", mcmc)
}

# NOT EXPOSED
setRestartSettings.Rcpp_MCMCAlgorithm <- function(mcmc, filename, samples, write.multiple=TRUE)
{
  mcmc$setRestartFileSettings(filename, samples, write.multiple)
}

convergence.test <- function(object, nsamples = 10, frac1 = 0.1, frac2 = 0.5, thin = 1, plot = FALSE, ...)
{
  UseMethod("convergence.test", object)
}

# NOT EXPOSED
convergence.test.Rcpp_MCMCAlgorithm <- function(object, nsamples = 10, frac1 = 0.1, frac2 = 0.5, thin = 1, plot = FALSE, ...)
{
  # TODO: extend to work with multiple chains once we have that capability.
  
  loglik.trace <- object$getLogLikelihoodTrace()
  trace.length <- length(loglik.trace)
  start <- max(1, trace.length - nsamples)
  
  # the start and end parameter do NOT work, using subsetting to achieve goal
  mcmcobj <- coda::mcmc(data=loglik.trace[start:trace.length], thin = thin)
  diag <- coda::geweke.diag(mcmcobj, frac1=frac1, frac2=frac2)
  if(plot){ 
    coda::geweke.plot(mcmcobj, frac1=frac1, frac2=frac2, ...)
  }else{
    return(diag)
  }
}


writeMCMCObject <- function(mcmc, file)
{
  loglikeTrace <- mcmc$getLogLikelihoodTrace()
  samples <- mcmc$getSamples()
  thining <- mcmc$getThining()
  adaptiveWidth <- mcmc$getAdaptiveWidth()
  save(list = c("loglikeTrace", "samples", "thining", "adaptiveWidth"), file=file)
}

loadMCMCObject <- function(file)
{
  mcmc <- new(MCMCAlgorithm)
  tempEnv <- new.env();
  load(file = file, envir = tempEnv)
  mcmc$setSamples(tempEnv$samples)
  mcmc$setThining(tempEnv$thining)
  mcmc$setAdaptiveWidth(tempEnv$adaptiveWidth)
  mcmc$setLogLikelihoodTrace(tempEnv$loglikeTrace)
  return(mcmc)
}


