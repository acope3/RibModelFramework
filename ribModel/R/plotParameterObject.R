# Plot functions for parameter object
# The generic plot function expects the parameter object
# and a string to the the function what has to be ploted.
# additional arguments are geneIndex, and aa to index function like
# getExpressionTraceForGene or 

plot.Rcpp_ROCParameter <- function(parameter, what=c("codonSpecificParameter", "MixtureProbability" ,"SPhi", "ExpectedPhi", "Expression"), 
                                   geneIndex=1, ...)
{
  if(what[1] == "codonSpecificParameter")
  {
    plotCodonSpecificParameters(parameter)
  }
  if(what[1] == "MixtureProbability")
  {
    plotMixtureProbability(parameter)
  }
  if(what[1] == "SPhi")
  {
    plotSPhiTrace(parameter)
  }
  if(what[1] == "ExpectedPhi")
  {
    plotExpectedPhiTrace(parameter)
  }
  if(what[1] == "Expression")
  {
    plotExpressionTrace(parameter, geneIndex)
  }
}

plotCodonSpecificParameters <- function(parameter, main="Codon specific parameter traces")
{
  opar <- par(no.readonly = T) 
  
  ### Trace plot.
  nf <- layout(matrix(c(rep(1, 4), 2:21), nrow = 6, ncol = 4, byrow = TRUE),
               rep(1, 4), c(2, 8, 8, 8, 8, 8), respect = FALSE)  
  #   nf <- layout(matrix(c(rep(1, 5), 2:21), nrow = 5, ncol = 5, byrow = TRUE),
  #                rep(1, 5), c(2, 8, 8, 8, 8), respect = FALSE)
  
  ### Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, main)
  text(0.5, 0.4, date(), cex = 0.6)
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  
  for(i.aa in names.aa){
    id.tmp <- grepl(paste(i.aa, i.aa, sep="."), names.b, fixed=T) & id.plot
    trace <- lapply(1:length(bMat), function(i){ bMat[[i]][id.tmp] })
    trace <- do.call("rbind", trace)
    if(length(trace) == 0) next
    ylim <- range(trace, na.rm=T)
    
    main.aa <- oneLetterAAtoThreeLetterAA(i.aa)
    plot(NULL, NULL, xlim = xlim, ylim = ylim,
         xlab = "Samples", ylab = ylab, main = main.aa)
    plot.order <- order(apply(trace, 2, sd), decreasing = TRUE)
    for(i.codon in plot.order){
      lines(x = x, y = trace[, i.codon], col = .CF.PT$color[i.codon])
    } 
  }
  
  
  par(opar)
}

plotExpressionTrace <- function(parameter, geneIndex)
{
  plot(parameter$getExpressionTraceForGene(geneIndex), type= "l", xlab = "Sample", ylab = expression(phi))
}

plotExpectedPhiTrace <- function(parameter)
{
  plot(parameter$getExpectedPhiTrace( ), type="l", xlab = "Sample", ylab = expression("E["~phi~"]"), 
       main = expression("Trace of the Expected value of "~phi))
  abline(h=1, col="red", lwd=1.5, lty=2)
}
plotSPhiTrace <- function(parameter)
{
  opar <- par(no.readonly = T) 
  par(oma=c(1,1,2,1), mgp=c(2,1,0), mar = c(3,4,2,1), mfrow=c(2, 1))
  sphi <- parameter$getSPhiTrace()
  mphi <- -(sphi * sphi) / 2;
  plot(mphi, type="l", xlab = "Sample", ylab = expression("m"[phi]))
  plot(sphi, type="l", xlab = "Sample", ylab = expression("s"[phi]))
  par(opar)
}
plotMixtureProbability <- function(parameter)
{
  samples <- length(parameter$getCategoryProbabilitiesTrace(0))
  numMixtures <- parameter$numMixtures
  
  plot(NULL, NULL, xlim = c(0, samples), ylim=c(0, 1), xlab = "Samples", ylab="Mixture Probability")
  for(i in 1:numMixtures)
  {
    lines(parameter$getCategoryProbabilitiesTrace(i-1), col= .mixtureColors[i])    
  }
  legend("topleft", legend = paste("Mixture Element", 1:numMixtures), 
         col = .mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")
}