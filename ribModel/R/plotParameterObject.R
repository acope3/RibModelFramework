# Plot functions for parameter object
# The generic plot function expects the parameter object
# and a string to the the function what has to be ploted.
# additional arguments are geneIndex, and category to index function like
# getExpressionTraceForGene or plotCodonSpecificParameters

plot.Rcpp_ROCParameter <- function(parameter, what=c("Mutation", "Selection", "MixtureProbability" ,"SPhi", "ExpectedPhi", "Expression"), 
                                   geneIndex=1, category = 1, ...)
{
  if(what[1] == "Mutation")
  {
    plotCodonSpecificParameters(parameter, category, "mutation", main="Mutation Parameter Traces")
  }
  if(what[1] == "Selection")
  {
    plotCodonSpecificParameters(parameter, category, "selection", main="Selection Parameter Traces")
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

plotCodonSpecificParameters <- function(parameter, category, type="mutation", main="Mutation Parameter Traces")
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
    trace <- vector("list", length(codons))
    if(type == "mutation")
    {
      ylab <- expression(Delta~"M")
      for(i in 1:length(codons))
      {
        trace[[i]] <- parameter$getMutationParameterTraceByCategoryForCodon(category, codons[i])
      }
    }else{
      ylab <- expression(Delta~eta)
      for(i in 1:length(codons))
      {
        trace[[i]] <- parameter$getSelectionParameterTraceByCategoryForCodon(category, codons[i])
      }
    }
    trace <- do.call("cbind", trace)
    if(length(trace) == 0) next
    x <- 1:dim(trace)[1]
    xlim <- range(x)
    ylim <- range(trace, na.rm=T)
    
    main.aa <- aa #TODO map to three leter code
    plot(NULL, NULL, xlim = xlim, ylim = ylim,
         xlab = "Samples", ylab = ylab, main = main.aa)
    plot.order <- order(apply(trace, 2, sd), decreasing = TRUE)
    for(i.codon in plot.order){
      lines(x = x, y = trace[, i.codon], col = ribModel:::.codonColors[[codons[i.codon]]])
    }
    colors <- unlist(ribModel:::.codonColors[codons])
    legend("topleft", legend = codons, col = colors, 
           lty = rep(1, length(codons)), bty = "n", cex = 0.75)
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
  samples <- length(parameter$getCategoryProbabilitiesTraceForCategory(1))
  numMixtures <- parameter$numMixtures
  
  plot(NULL, NULL, xlim = c(0, samples), ylim=c(0, 1), xlab = "Samples", ylab="Mixture Probability")
  for(i in 1:numMixtures)
  {
    lines(parameter$getCategoryProbabilitiesTraceForCategory(i), col= ribModel:::.mixtureColors[i])    
  }
  legend("topleft", legend = paste("Mixture Element", 1:numMixtures), 
         col = .mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")
}