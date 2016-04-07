# Plot functions for trace object
# The generic plot function expects the trace object
# and a string to the the function what has to be ploted.
# additional arguments are geneIndex, and category to index function like
# getExpressionTraceForGene or plotCodonSpecificParameters

#' Plot Trace Object
#' @param x An Rcpp trace object initialized with \code{initializeTraceObject}.
#' @param what A string containing one of the following to graph: \code{Mutation, Selection, MixtureProbability, Sphi,
#'  Mphi, Aphi, Spesilon, ExpectedPhi, Expression}.
#' @param geneIndex When plotting expression, the index of the gene to be plotted.
#' @param mixture The mixture for which to plot values.
#' 
#' @return This function has no return value.
#' 
#' @description Plots different traces, specified with the \code{what} parameter.
plot.Rcpp_Trace <- function(x, what=c("Mutation", "Selection", "MixtureProbability" ,"Sphi", "Mphi", "Aphi", "Sepsilon", "ExpectedPhi", "Expression"), 
                                   geneIndex=1, mixture = 1, ...)
{
  if(what[1] == "Mutation")
  {
    plotCodonSpecificParameters(x, mixture, "mutation", main="Mutation Parameter Traces")
  }
  if(what[1] == "Selection")
  {
    plotCodonSpecificParameters(x, mixture, "selection", main="Selection Parameter Traces")
  }
  if(what[1] == "Alpha")
  {
    plotCodonSpecificParameters(x, mixture, "alpha", main="Alpha Parameter Traces", ROC=FALSE)
  }
  if(what[1] == "LambdaPrime")
  {
    plotCodonSpecificParameters(x, mixture, "lambdaPrime", main="LambdaPrime Parameter Traces", ROC=FALSE)
  }  
  if(what[1] == "MixtureProbability")
  {
    plotMixtureProbability(x)
  }
  if(what[1] == "Sphi")
  {
    plotHyperParameterTrace(x, what = what[1])
  }
  if(what[1] == "Mphi") 
  {
    plotHyperParameterTrace(x, what = what[1])
  }
  if(what[1] == "Aphi")
  {
    plotHyperParameterTrace(x, what = what[1])
  }
  if(what[1] == "Sepsilon") 
  {
    plotHyperParameterTrace(x, what = what[1])
  }
  if(what[1] == "ExpectedPhi")
  {
    plotExpectedPhiTrace(x)
  }
  if(what[1] == "Expression")
  {
    plotExpressionTrace(x, geneIndex)
  }
}

# NOT EXPOSED
plotCodonSpecificParameters <- function(trace, mixture, type="mutation", main="Mutation Parameter Traces", ROC=TRUE)
{
  opar <- par(no.readonly = T) 
  
  ### Trace plot.
  if (ROC)
  {
    nf <- layout(matrix(c(rep(1, 4), 2:21), nrow = 6, ncol = 4, byrow = TRUE),
               rep(1, 4), c(2, 8, 8, 8, 8, 8), respect = FALSE)  
  }else
  {    
    nf <- layout(matrix(c(rep(1, 4), 2:25), nrow = 7, ncol = 4, byrow = TRUE),
                    rep(1, 4), c(2, 8, 8, 8, 8, 8, 8), respect = FALSE) 
  }
  ### Plot title.
  if (ROC){
  par(mar = c(0, 0, 0, 0))
  }else{
    par(mar = c(1,1,1,1))
  }
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, main)
  text(0.5, 0.4, date(), cex = 0.6)
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  
  # TODO change to groupList -> checks for ROC like model is not necessary!
  names.aa <- aminoAcids()
  with.ref.codon <- ifelse(ROC, TRUE, FALSE)
  
  for(aa in names.aa)
  { 
    codons <- AAToCodon(aa, with.ref.codon)
    if(length(codons) == 0) next
    if (!ROC){
      if(aa == "X") next
    }
    cur.trace <- vector("list", length(codons))
    paramType <- 0
    if(type == "mutation"){
      ylab <- expression(Delta~"M")
      paramType <- 0
    }else if (type == "selection"){
      ylab <- expression(Delta~eta)
      paramType <- 1
    }else if (type == "alpha"){
      ylab <- expression(alpha)
      paramType <- 0
    }else if (type == "lambdaPrime"){
      ylab <- expression(lambda~"'")
      paramType <- 1
    }else{
      stop("parameter type not recognized")
    }
    for(i in 1:length(codons))
    {
      cur.trace[[i]] <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], paramType, with.ref.codon)
    }      

    cur.trace <- do.call("cbind", cur.trace)
    if(length(cur.trace) == 0) next
    x <- 1:dim(cur.trace)[1]
    xlim <- range(x)
    ylim <- range(cur.trace, na.rm=T)
    
    main.aa <- aa #TODO map to three leter code
    plot(NULL, NULL, xlim = xlim, ylim = ylim,
         xlab = "Samples", ylab = ylab, main = main.aa)
    plot.order <- order(apply(cur.trace, 2, sd), decreasing = TRUE)
    for(i.codon in plot.order){
      lines(x = x, y = cur.trace[, i.codon], col = .codonColors[[codons[i.codon]]])
    }
    colors <- unlist(.codonColors[codons])
    legend("topleft", legend = codons, col = colors, 
           lty = rep(1, length(codons)), bty = "n", cex = 0.75)
  }
  par(opar)
} 

# NOT EXPOSED
plotExpressionTrace <- function(trace, geneIndex)
{
  plot(log10(trace$getSynthesisRateTraceForGene(geneIndex)), type= "l", xlab = "Sample", ylab = expression("log"[10]~"("~phi~")"))
}

# NOT EXPOSED
plotExpectedPhiTrace <- function(trace)
{
  par(mar=c(5,5,4,2))
  plot(trace$getExpectedSynthesisRateTrace()[-1], type="l", xlab = "Sample", ylab = expression(bar(phi)), 
       main = expression("Trace of the Expected value of "~phi))
  abline(h=1, col="red", lwd=1.5, lty=2)
}

# NOT EXPOSED
plotHyperParameterTrace <- function(trace, what = c("Sphi", "Mphi", "Aphi", "Sepsilon"))
{
#  opar <- par(no.readonly = T) 
#  par(oma=c(1,1,2,1), mgp=c(2,1,0), mar = c(3,4,2,1), mfrow=c(2, 1))
  if (what[1] == "Sphi")
  {
    sphi <- trace$getStdDevSynthesisRateTraces();
    numMixtures <- length(sphi)
    sphi <- do.call("cbind", sphi)
    ylimit <- range(sphi) + c(-0.1, 0.1)
    xlimit <- c(1, nrow(sphi))
    plot(NULL, NULL, type="l", xlab = "Sample", ylab = expression("s"[phi]), xlim  = xlimit, ylim = ylimit)
    for(i in 1:ncol(sphi))
    {
      lines(sphi[-1,i], col = .mixtureColors[i])
    }
    legend("topleft", legend = paste("Mixture Element", 1:numMixtures), 
           col = .mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")
  }
  if (what[1] == "Mphi")
  {
    sphi <- trace$getStdDevSynthesisRateTraces();
    numMixtures <- length(sphi)
    sphi <- do.call("cbind", sphi)
    mphi <- -(sphi * sphi) / 2;
    ylimit <- range(mphi) + c(-0.1, 0.1)
    xlimit <- c(1, nrow(mphi))
    plot(NULL, NULL, type="l", xlab = "Sample", ylab = expression("m"[phi]), xlim  = xlimit, ylim = ylimit)
    for(i in 1:ncol(mphi))
    {
      lines(mphi[-1,i], col= .mixtureColors[i])
    }
    legend("topleft", legend = paste("Mixture Element", 1:numMixtures), 
           col = .mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")    

  }
  if (what[1] == "Aphi") 
  {
    aphi <- trace$getSynthesisOffsetTrace();
    aphi <- do.call("cbind", aphi)
    ylimit <- range(aphi) + c(-0.1, 0.1)
    xlimit <- c(1, nrow(aphi))
    plot(NULL, NULL, type="l", xlab = "Sample", ylab = expression("A"[phi]), xlim  = xlimit, ylim = ylimit)
    for(i in 1:ncol(aphi))
    {
      lines(aphi[-1,i], col = .mixtureColors[i])
    }
    legend("topleft", legend = paste("Observed Data", 1:numMixtures), 
           col = .mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")        
  }
  if (what[1] == "Sepsilon")
  {
    sepsilon <- trace$getObservedSynthesisNoiseTrace();
    sepsilon <- do.call("cbind", sepsilon)
    ylimit <- range(sepsilon) + c(-0.1, 0.1)
    xlimit <- c(1, nrow(sepsilon))
    plot(NULL, NULL, type="l", xlab = "Sample", ylab = expression("s"[epsilon]), xlim  = xlimit, ylim = ylimit)
    for(i in 1:ncol(sepsilon))
    {
      lines(sepsilon[-1,i], col = .mixtureColors[i])
    }
    legend("topleft", legend = paste("Observed Data", 1:numMixtures), 
           col = .mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")  
  }
  #par(opar)
}

# NOT EXPOSED
plotMixtureProbability <- function(trace)
{
  samples <- length(trace$getMixtureProbabilitiesTraceForMixture(1))
  numMixtures <- trace$getNumberOfMixtures()
  
  plot(NULL, NULL, xlim = c(0, samples), ylim=c(0, 1), xlab = "Samples", ylab="Mixture Probability")
  for(i in 1:numMixtures)
  {
    lines(trace$getMixtureProbabilitiesTraceForMixture(i)[-1], col = .mixtureColors[i])    
  }
  legend("topleft", legend = paste("Mixture Element", 1:numMixtures), 
         col = .mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")
}