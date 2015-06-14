# Plot functions for trace object
# The generic plot function expects the trace object
# and a string to the the function what has to be ploted.
# additional arguments are geneIndex, and category to index function like
# getExpressionTraceForGene or plotCodonSpecificParameters

plot.Rcpp_ROCTrace <- function(trace, what=c("Mutation", "Selection", "MixtureProbability" ,"SPhi", "ExpectedPhi", "Expression"), 
                                   geneIndex=1, mixture = 1, ...)
{
  if(what[1] == "Mutation")
  {
    plotCodonSpecificParameters(trace, mixture, "mutation", main="Mutation Parameter Traces")
  }
  if(what[1] == "Selection")
  {
    plotCodonSpecificParameters(trace, mixture, "selection", main="Selection Parameter Traces")
  }  
  if(what[1] == "MixtureProbability")
  {
    plotMixtureProbability(trace)
  }
  if(what[1] == "SPhi")
  {
    plotSPhiTrace(trace)
  }
  if(what[1] == "ExpectedPhi")
  {
    plotExpectedPhiTrace(trace)
  }
  if(what[1] == "Expression")
  {
    plotExpressionTrace(trace, geneIndex)
  }
}

plotCodonSpecificParameters <- function(trace, mixture, type="mutation", main="Mutation Parameter Traces")
{
  opar <- par(no.readonly = T) 
  
  ### Trace plot.
  nf <- layout(matrix(c(rep(1, 4), 2:21), nrow = 6, ncol = 4, byrow = TRUE),
               rep(1, 4), c(2, 8, 8, 8, 8, 8), respect = FALSE)  
  
  ### Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, main)
  text(0.5, 0.4, date(), cex = 0.6)
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  
  names.aa <- aminoAcids()
  for(aa in names.aa)
  {  
    codons <- AAToCodon(aa, T)
    if(length(codons) == 0) next
    cur.trace <- vector("list", length(codons))
    if(type == "mutation")
    {
      ylab <- expression(Delta~"M")
      for(i in 1:length(codons))
      {
        cur.trace[[i]] <- trace$getMutationParameterTraceByMixtureElementForCodon(mixture, codons[i])
      }
    }else{
      ylab <- expression(Delta~eta)
      for(i in 1:length(codons))
      {
        cur.trace[[i]] <- trace$getSelectionParameterTraceByMixtureElementForCodon(mixture, codons[i])
      }
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
      lines(x = x, y = cur.trace[, i.codon], col = ribModel:::.codonColors[[codons[i.codon]]])
    }
    colors <- unlist(ribModel:::.codonColors[codons])
    legend("topleft", legend = codons, col = colors, 
           lty = rep(1, length(codons)), bty = "n", cex = 0.75)
  }
  par(opar)
}

plotExpressionTrace <- function(trace, geneIndex)
{
  plot(log10(trace$getExpressionTraceForGene(geneIndex)), type= "l", xlab = "Sample", ylab = expression("log"[10]~"("~phi~")"))
}

plotExpectedPhiTrace <- function(trace)
{
  plot(trace$getExpectedPhiTrace( ), type="l", xlab = "Sample", ylab = expression("E["~phi~"]"), 
       main = expression("Trace of the Expected value of "~phi))
  abline(h=1, col="red", lwd=1.5, lty=2)
}
plotSPhiTrace <- function(trace)
{
  opar <- par(no.readonly = T) 
  par(oma=c(1,1,2,1), mgp=c(2,1,0), mar = c(3,4,2,1), mfrow=c(2, 1))
  sphi <- trace$getSPhiTrace()
  mphi <- -(sphi * sphi) / 2;
  plot(mphi, type="l", xlab = "Sample", ylab = expression("m"[phi]))
  plot(sphi, type="l", xlab = "Sample", ylab = expression("s"[phi]))
  par(opar)
}
plotMixtureProbability <- function(trace)
{
  samples <- length(trace$getMixtureProbabilitiesTraceForMixture(1))
  numMixtures <- trace$getNumberOfMixtures()
  
  plot(NULL, NULL, xlim = c(0, samples), ylim=c(0, 1), xlab = "Samples", ylab="Mixture Probability")
  for(i in 1:numMixtures)
  {
    lines(trace$getMixtureProbabilitiesTraceForMixture(i), col= ribModel:::.mixtureColors[i])    
  }
  legend("topleft", legend = paste("Mixture Element", 1:numMixtures), 
         col = ribModel:::.mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")
}