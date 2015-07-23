plot.Rcpp_ROCParameter <- function(parameter, what = "Mutation", sample = 100, ...)
{
  
  numMixtures <- parameter$numMixtures
  csp.params <- data.frame(matrix(0, ncol=numMixtures, nrow = 40))
  
  names.aa <- aminoAcids()
  for(mixture in 1:numMixtures)
  {
    
    param.storage <- vector("numeric", 0)
    param.name.storage <- vector("numeric", 0)
    # get codon specific parameter
    for(aa in names.aa)
    {
      if(aa == "M" || aa == "W" || aa == "X") next
      codons <- AAToCodon(aa, T)
      for(i in 1:length(codons))
      {
        if(what == "Mutation")
        {
          param.storage <- c(param.storage, parameter$getMutationPosteriorMeanForCodon(mixture, samples, codons[i]))
        }else{
          param.storage <- c(param.storage, parameter$getSelectionPosteriorMeanForCodon(mixture, samples, codons[i]))
        }
        param.name.storage <- c(param.name.storage, paste(aa, codons[i], sep="."))
      }
    }
    csp.params[, mixture] <- param.storage
  }
  rownames(csp.params) <- param.name.storage
  colnames(csp.params) <- paste("Mixture\nElement", 1:numMixtures)
  pairs(csp.params, upper.panel = upper.panel.plot, lower.panel=NULL, main = ...)
}


upper.panel.plot <- function(x, y, ...)
{
  abline(0, 1, col = "blue", lty = 2)
  points(x, y, ...)
  
  lm.line <- lm(y~x)
  abline(lm.line, col="blue", lwd = 2)
  
  R2 <- summary(lm.line)$r.squared
  
  b <- lm.line$coef[2]
  rho <- ifelse(b > 0, sqrt(R2), -sqrt(R2)) #make sure rho has correct sign
  
  xlim <- range(x)
  ylim <- range(y)
  
  width <- xlim[2] - xlim[1]
  height <- ylim[2] - ylim[1]
  
  std.error <- summary(lm.line)$coefficients[4]
  slope <- round(summary(lm.line)$coefficients[2], 3)
  intercept <- round(summary(lm.line)$coefficients[1], 3)
  t <- (slope - 1)/std.error
  
  if(t > qt(1-0.05/2, lm.line$df.residual - 1))
  {
    eq <- paste("y = ", sprintf("%.3f", intercept), " + ", sprintf("%.3f", slope), "x *", sep = "")
    text(xlim[1] + width * 0.1, ylim[2] - height * 0.2, eq)
  }else{
    eq <- paste("y = ", sprintf("%.3f", intercept), " + ", sprintf("%.3f", slope), sep = "")
    text(xlim[1] + width * 0.1, ylim[2] - height * 0.2, eq)
  } 
  text(xlim[2] - width * 0.04, ylim[1] + height * 0.05,
       parse(text = paste("rho == ", sprintf("%.4f", rho), sep = "")),
       pos = 2, cex = 1.0, font = 2)
}
lower.panel.plot <- function(x, y, ...)
{
  
}