#' Plot MCMC algorithm
#' 
#' @param x An Rcpp_MCMC object initialized with \code{initializeMCMCObject}.
#' 
#' @param zoom.window A vector describing the start and end of the zoom window.
#' 
#' @param ... Arguments to be passed to methods, such as graphical parameters.
#' 
#' @return This function has no return value.
#' 
#' @description This function will plot the logLikelihood trace, and if the Hmisc package is installed, it will 
#'  plot a subplot of the logLikelihood trace with the first few samples removed.
plot.Rcpp_MCMCAlgorithm <- function(x, what=c("LogPosterior","LogLikelihood"),zoom.window = NULL, ...)
{
  if(what[1] == "LogPosterior")
  {
    trace <- x$getLogPosteriorTrace()
    ylab= "log(Posterior Probability)"
  }
  else if(what[1] == "LogLikelihood")
  {
    trace <- x$getLogLikelihoodTrace()
    ylab= "log(Likelihood Probability)"
  }
  trace <- trace[-1]
  
  trace.length <- length(trace)
  
  zoomStart <- round(0.9*trace.length)
  zoomEnd <- trace.length
  logL <- mean(trace[zoomStart:trace.length])
  #TODO change main title
  plot(trace, type="l", main=paste0("logPP:", logL), xlab="Sample", ylab=ylab)
  grid (NULL,NULL, lty = 6, col = "cornsilk2")
  trace[trace == -Inf] <- NA
  
  # TODO (Cedric): get rid of that line once problem with first element beeing 0 is solved
  trace <- trace[-1]
  
  if(!(is.null(zoom.window))) {
    zoomStart <- zoom.window[1]
    zoomEnd <- zoom.window[2]
  }
  else{
    warning("No window was given, zooming in at last 10% of trace")
  }
  
  Hmisc::subplot(
    plot(zoomStart:zoomEnd, trace[zoomStart:zoomEnd], type="l", xlab=NA, ylab=NA, las=2, cex.axis=0.55), 
    0.8*(round(0.9*trace.length)), (min(trace, na.rm = T)+max(trace, na.rm = T))/2, size=c(3,2))
}

#' Plots ACF for MCMC traces
#' @param mcmc object of class MCMC
#' @param type "LogPosterior" or "LogLikelihood", defaults to "LogPosterior"
#' @param samples number of samples to calculate ACF for,
#' Ex) if samples == 300 and hae 1000 samples total, ACF will be calculated from samples 700 to 1000
#' @param lag.max Maximum amount of lag to calculate ACF

mcmcACF <- function(mcmc,type="LogPosterior",samples=500,lag.max = 40)
{
  if(type == "Posterior")
  {
    trace <- mcmc$getLogPosteriorTrace()
  }else{
    trace <- mcmc$getLogLikelihoodTrace()
  }
  trace <- trace[length(trace)-samples:length(trace)]
  trace.acf <- acf(x = trace,lag.max = lag.max,plot = FALSE)
  header <- paste(type,"Trace Autocorrelation",sep=" ")
  plot(x = trace.acf,xlab = "Lag time",ylab = "Autocorrelation",main = header)
}
