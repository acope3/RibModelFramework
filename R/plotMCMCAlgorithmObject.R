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
plot.Rcpp_MCMCAlgorithm <- function(x, zoom.window = NULL, ...)
{
  loglik.trace <- x$getLogLikelihoodTrace()
  
  loglik.trace <- loglik.trace[-1]
  
  trace.length <- length(loglik.trace)
  
  zoomStart <- round(0.9*trace.length)
  zoomEnd <- trace.length
  logL <- mean(loglik.trace[zoomStart:trace.length])
  #TODO change main title
  plot(loglik.trace, type="l", main=paste0("logPP:", logL), xlab="Sample", ylab="log(Posterior Probability)")
  grid (NULL,NULL, lty = 6, col = "cornsilk2")
  loglik.trace[loglik.trace == -Inf] <- NA
  
  # TODO (Cedric): get rid of that line once problem with first element beeing 0 is solved
  loglik.trace <- loglik.trace[-1]
  
  if(!(is.null(zoom.window))) {
    zoomStart <- zoom.window[1]
    zoomEnd <- zoom.window[2]
  }
  else{
    print("No window was given, zooming in at last 10% of trace")
  }
  
  Hmisc::subplot(
    plot(zoomStart:zoomEnd, loglik.trace[zoomStart:zoomEnd], type="l", xlab=NA, ylab=NA, las=2, cex.axis=0.55), 
    0.8*(round(0.9*trace.length)), (min(loglik.trace, na.rm = T)+max(loglik.trace, na.rm = T))/2, size=c(3,2))
}

