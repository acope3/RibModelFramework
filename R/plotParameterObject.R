#' Plot ROC Parameter 
#' 
#' @param x ROC parameter object.
#' 
#' @param what Which aspect of the ROC parameter to plot. Default value is
#' "Mutation".
#' 
#' @param samples Number of samples to plot using the posterior mean. Default
#' value is 100.
#' 
#' @param ... Arguments to be passed to methods, such as graphical parameters.
#' 
#' @return This function has no return value.
#' 
#' @description \code{plot} graphs the mutation or selection parameter for a ROC
#' parameter object for each mixture element. 
#' 
#' @details Graphs are based off the last # samples for the posterior mean.
#' 
#' 
plot.Rcpp_ROCParameter <- function(x, what = "Mutation", samples = 100, ...)
{
  plotParameterObject(x, what = what, samples = 100, ...)
}
plot.Rcpp_FONSEParameter <- function(x, what = "Mutation", samples = 100, ...)
{
  plotParameterObject(x, what = what, samples = 100, ...)
}

plotParameterObject <- function(x, what = "Mutation", samples = 100, ...){
  
  numMixtures <- x$numMixtures
  csp.params <- data.frame(matrix(0, ncol=numMixtures, nrow = 40))
  
  names.aa <- aminoAcids()
  paramType <- ifelse(what == "Mutation", 0, 1)
  cat("ParamType: ", paramType, "\n")
  for(mixture in 1:numMixtures){
    param.storage <- vector("numeric", 0)
    param.name.storage <- vector("numeric", 0)
    # get codon specific parameter
    for(aa in names.aa){
      if(aa == "M" || aa == "W" || aa == "X") next
      codons <- AAToCodon(aa, T)
      for(i in 1:length(codons)){
        param.storage <- c(param.storage, x$getCodonSpecificPosteriorMean(mixture, samples, codons[i], paramType))
      }
    }
    csp.params[, mixture] <- param.storage
  }
  #rownames(csp.params) <- param.name.storage
  colnames(csp.params) <- paste("Mixture\nElement", 1:numMixtures)
  pairs(csp.params, upper.panel = upper.panel.plot, lower.panel=NULL)
}



#TODO: should RFP's ploting be here as well?

upper.panel.plot <- function(x, y, sd.x=NULL, sd.y=NULL, ...){
  abline(0, 1, col = "blue", lty = 2)
  points(x, y, ...)
  if(!is.null(sd.y)){
    y.up <- sd.y[,2]
    y.low <- sd.y[,1]
    epsilon <- range(x, na.rm = T) * 0.1
    segments(x, y.low, x, y.up, ...)
  }
  if(!is.null(sd.x)){
    x.up <- sd.x[,2]
    x.low <- sd.x[,1]
    epsilon <- range(y, na.rm = T) * 0.1
    segments(x.low, y, x.up, y, ...)
  }  
  
  lm.line <- lm(y~x, na.action = "na.exclude")
  abline(lm.line, col="blue", lwd = 2)
  
  R2 <- summary(lm.line)$r.squared
  
  b <- lm.line$coef[2]
  rho <- ifelse(b > 0, sqrt(R2), -sqrt(R2)) #make sure rho has correct sign
  
  xlim <- range(x, na.rm = T)
  ylim <- range(y, na.rm = T)
  
  width <- xlim[2] - xlim[1]
  height <- ylim[2] - ylim[1]
  
  std.error <- summary(lm.line)$coefficients[4]
  slope <- round(summary(lm.line)$coefficients[2], 3)
  intercept <- round(summary(lm.line)$coefficients[1], 3)
  t <- (slope - 1)/std.error
  
  if(t > qt(1-0.05/2, lm.line$df.residual - 1)){
    eq <- paste("y = ", sprintf("%.3f", intercept), " + ", sprintf("%.3f", slope), "x *", sep = "")
    text(xlim[1] + width * 0.1, ylim[2] - height * 0.2, eq)
  }else{
    eq <- paste("y = ", sprintf("%.3f", intercept), " + ", sprintf("%.3f", slope), "x", sep = "")
    text(xlim[1] + width * 0.1, ylim[2] - height * 0.2, eq)
  } 
  text(xlim[2] - width * 0.04, ylim[1] + height * 0.05,
       parse(text = paste("rho == ", sprintf("%.4f", rho), sep = "")),
       pos = 2, cex = 1.0, font = 2)
}


lower.panel.plot <- function(x, y, ...)
{
  
}


confidenceInterval.plot <- function(x, y, sd.x=NULL, sd.y=NULL, ...){
  points(x, y, ...)
  if(!is.null(sd.y)){
    y.up <- sd.y[,2]
    y.low <- sd.y[,1]
    epsilon <- range(x, na.rm = T) * 0.1
    segments(x, y.low, x, y.up, ...)
  }
  if(!is.null(sd.x)){
    x.up <- sd.x[,2]
    x.low <- sd.x[,1]
    epsilon <- range(y, na.rm = T) * 0.1
    segments(x.low, y, x.up, y, ...)
  }  
  
  lm.line <- lm(y~x, na.action = "na.exclude")
  
  
  b <- lm.line$coef[2]
  
  xlim <- range(x, na.rm = T)
  ylim <- range(y, na.rm = T)
  
  width <- xlim[2] - xlim[1]
  height <- ylim[2] - ylim[1]
  
  std.error <- summary(lm.line)$coefficients[4]
  slope <- round(summary(lm.line)$coefficients[2], 3)
  intercept <- round(summary(lm.line)$coefficients[1], 3)
  t <- (slope - 1)/std.error
  

}