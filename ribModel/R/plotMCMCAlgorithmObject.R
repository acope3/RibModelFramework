

plot.Rcpp_MCMCAlgorithm <- function(mcmc, ...)
{
  loglik.trace <- mcmc$getLogLikelihoodTrace()
  trace.length <- length(loglik.trace)
  zoomStart <- round(0.9*trace.length)
  logL <- mean(loglik.trace[zoomStart:trace.length])
  #TODO change main title
  plot(loglik.trace, type="l", main=paste("logL:", logL), xlab="Sample", ylab="log(Likelihood)")
  grid (NULL,NULL, lty = 6, col = "cornsilk2")
  loglik.trace[loglik.trace == -Inf] <- NA
  Hmisc::subplot(
    plot(zoomStart:trace.length, loglik.trace[zoomStart:trace.length], type="l", xlab=NA, ylab=NA, las=2, cex.axis=0.55), 
    0.8*zoomStart, (min(loglik.trace, na.rm = T)+max(loglik.trace, na.rm = T))/2, size=c(3,2))
}

