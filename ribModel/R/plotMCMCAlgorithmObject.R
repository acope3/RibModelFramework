

plot.Rcpp_MCMCAlgorithm <- function(x, ...)
{
  loglik.trace <- x$getLogLikelihoodTrace()
  
  loglik.trace <- loglik.trace[-1]
  
  trace.length <- length(loglik.trace)
  zoomStart <- round(0.9*trace.length)
  logL <- mean(loglik.trace[zoomStart:trace.length])
  #TODO change main title
  plot(loglik.trace, type="l", main=paste("logPP:", logL), xlab="Sample", ylab="log(Posterior Probability)")
  grid (NULL,NULL, lty = 6, col = "cornsilk2")
  loglik.trace[loglik.trace == -Inf] <- NA
  
  # TODO (Cedric): get rid of that line once problem with first element beeing 0 is solved
  loglik.trace <- loglik.trace[-1]
  
  
  Hmisc::subplot(
    plot(zoomStart:trace.length, loglik.trace[zoomStart:trace.length], type="l", xlab=NA, ylab=NA, las=2, cex.axis=0.55), 
    0.8*zoomStart, (min(loglik.trace, na.rm = T)+max(loglik.trace, na.rm = T))/2, size=c(3,2))
}

