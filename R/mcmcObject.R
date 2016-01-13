#' Initialize MCMC 
#' 
#' @param samples Number of samples that should be produced when running the 
#' MCMC algorithm. No default value.
#' 
#' @param thining Number that the samples should be thinned by. If set to 
#' 1, every step would be saved as a sample. Default value is 1.
#' 
#' @param adaptive.width Number that determines how often the acceptance/rejection
#' window should be altered. Default value is 100.
#' 
#' @param est.expression Boolean that tells whether or not synthesis rate values
#' should be estimated in the MCMC algorithm run. Default value is TRUE.
#' 
#' @param est.csp Boolean that tells whether or not codon specific values
#' should be estimated in the MCMC algorithm run. Default value is TRUE.
#' 
#' @param est.hyper Boolean that tells whether or not hyper parameters
#' should be estimated in the MCMC algorithm run. Default value is TRUE.
#' 
#' @param est.mix Boolean that tells whether or not the genes' mixture element
#' should be estimated in the MCMC algorithm run. Default value is TRUE.
#' 
#' @return mcmc Returns an intialized MCMC object. 
#' 
#' @description \code{initializeMCMCObject} initializes a MCMC object to 
#' perform a model fitting for a parameter and model object.
#' 
#' @details \code{initializeMCMCObject} sets up the MCMC object 
#' (monte carlo markov chain) and returns the object so a model fitting can be done.
#' It is important to note that est.expression and est.hyper will affect one another
#' negatively if their values differ.
#' 
initializeMCMCObject <- function(samples, thining=1, adaptive.width=100, 
                                 est.expression=TRUE, est.csp=TRUE, 
                                 est.hyper=TRUE, est.mix=TRUE){
  
  #TODO: error check given values.
  mcmc <- new(MCMCAlgorithm, samples, thining, adaptive.width, est.expression, 
              est.csp, est.hyper)
  mcmc$setEstimateMixtureAssignment(est.mix)
  return(mcmc)
}


#' Run MCMC 
#' 
#' @param mcmc MCMC object that will run the model fitting algorithm.
#' 
#' @param genome Genome that the model fitting will run on. Should be 
#' the same genome associated with the parameter and model objects.
#' 
#' @param model Model to run the fitting on. Should be associated with
#' the given genome.
#' 
#' @param ncores Number of cores to perform the model fitting with. Default
#' value is 1.
#' 
#' @param divergence.iteration Number of steps that the initial conditions
#' can diverge from the original conditions given. Default value is 0.
#' 
#' @return This function has no return value.
#' 
#' @description \code{runMCMC} will run a monte carlo markov chain algorithm
#' for the given mcmc, genome, and model objects to perform a model fitting.
#' 
#' @details \code{runMCMC} will run for the number of samples times the number
#' thining given when the mcmc object is initialized. Updates are provided every 100
#' steps, and the state of the chain is saved every thining steps.
#' 
runMCMC <- function(mcmc, genome, model, ncores = 1, divergence.iteration = 0){
  
  #TODO: error check values
  UseMethod("runMCMC", mcmc)
}


#Called from "runMCMC."
runMCMC.Rcpp_MCMCAlgorithm <- function(mcmc, genome, model, ncores = 1, 
                                       divergence.iteration = 0){
  mcmc$run(genome, model, ncores, divergence.iteration)
}







#' Set Restart Settings 
#' 
#' @param mcmc MCMC object that will run the model fitting algorithm.
#' 
#' @param filename Filename for the restart files to be written.
#' 
#' @param samples Number of samples that should occur before a file is written.
#' 
#' @param write.multiple Boolean that determines if multiple restart files
#' are written. Default value is TRUE.
#' 
#' @return This function has no return value.
#' 
#' @description \code{setRestartSettings} sets the needed information (what the file 
#' is called, how often the file should be written) to write
#' information to restart the MCMC algorithm from a given point.
#' 
#' @details \code{setRestartSettings} writes a restart file every set amount of samples
#' that occur. Also, if write.multiple is true, instead of overwriting the previous restart
#' file, the sample number is prepended onto the file name and multiple rerstart files
#' are generated for a run.
#' 
setRestartSettings <- function(mcmc, filename, samples, write.multiple=TRUE){
  UseMethod("setRestartSettings", mcmc)
}


setRestartSettings.Rcpp_MCMCAlgorithm <- function(mcmc, filename, samples, 
                                                  write.multiple=TRUE){
  mcmc$setRestartFileSettings(filename, samples, write.multiple)
}
#TODO: Why is this seperated into 2 functions?




#' Convergence Test
#' 
#' @param object
#' 
#' @param nsamples
#' 
#' @param frac1
#' 
#' @param frac2
#' 
#' @param thin
#' 
#' @param plot
#' 
#' @param ...
#' 
#' @return 
#' 
#' @description \code{convergence.test} 
#' 
#' @details \code{convergence.test}
#' 
convergence.test <- function(object, nsamples = 10, frac1 = 0.1, frac2 = 0.5, 
                    thin = 1, plot = FALSE, ...){
  UseMethod("convergence.test", object)
}


convergence.test.Rcpp_MCMCAlgorithm <- function(object, nsamples = 10, frac1 = 0.1, 
                                       frac2 = 0.5, thin = 1, plot = FALSE, ...){
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


#' Write MCMC Object 
#' 
#' @param mcmc MCMC object that has run the model fitting algorithm.
#' 
#' @param file A filename that where the data will be stored.
#' The file should end with the extension "Rdat".
#' 
#' @return This function has no return value.
#' 
#' @description \code{writeMCMCObject} stores the MCMC information from the 
#' model fitting run in a file.
#' 
#' @details None.
#' 
writeMCMCObject <- function(mcmc, file){
  loglikeTrace <- mcmc$getLogLikelihoodTrace()
  samples <- mcmc$getSamples()
  thining <- mcmc$getThining()
  adaptiveWidth <- mcmc$getAdaptiveWidth()
  save(list = c("loglikeTrace", "samples", "thining", "adaptiveWidth"), file=file)
}


#' Load MCMC Object
#' 
#' @param file A filename that where the data will be stored.
#' 
#' @return This function has no return value.
#' 
#' @description \code{loadMCMCObject} creates a new MCMC object and fills it with
#' the information in the file given.
#' 
#' @details This MCMC object is not intended to be used to do another model fitting, only
#' to graph the stored results.
#' 
loadMCMCObject <- function(files){
  mcmc <- new(MCMCAlgorithm)
  samples <- 0
  loglikeTrace <- numeric(0)
  for (i in 1:length(files)){
    tempEnv <- new.env();
    load(file = files[i], envir = tempEnv)
    samples <- samples + tempEnv$samples
    max <- tempEnv$samples + 1
    curLoglikelihoodTrace <- tempEnv$loglikeTrace
    loglikeTrace <- c(loglikeTrace, curLoglikelihoodTrace[2:max])
   }
    mcmc$setSamples(samples)
    mcmc$setThining(tempEnv$thining) #not needed?
    mcmc$setAdaptiveWidth(tempEnv$adaptiveWidth) #not needed?
    mcmc$setLogLikelihoodTrace(loglikeTrace)

  return(mcmc)
}


