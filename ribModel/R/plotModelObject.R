    


plot.Rcpp_ROCModel <- function(model, genome, parameter, samples = 100, category = 1, ...)
{
  input_list <- as.list(list(...))
  
  if("main" %in% names(input_list)){
    main <- input_list$main
    input_list$main <- NULL
  }else{
    main <- ""
  }
  
  mat <- matrix(c(rep(1, 4), 2:21, rep(22, 4)),
                nrow = 7, ncol = 4, byrow = TRUE)
  mat <- cbind(rep(23, 7), mat, rep(24, 7))
  nf <- layout(mat, c(3, rep(8, 4), 2), c(3, 8, 8, 8, 8, 8, 3), respect = FALSE)
  ### Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, main)
  text(0.5, 0.4, date(), cex = 0.6)
  
  num.genes <- genome$getGenomeSize()
  
  expressionCategory <- round(unlist(lapply(1:num.genes,  function(geneIndex){parameter$getMixtureAssignmentPosteriorMean(samples, geneIndex-1)})))
  genes.in.category <- which(expressionCategory == category)
  # TODO expressionCategory is actually mixtureElement, mapping to expressionCategory is missing!
  num.genes <- length(genes.in.category)
  expressionValues <- vector("numeric", num.genes)
  i <- 1
  for(geneIndex in genes.in.category)
  {
    expressionValues[i] <- parameter$getExpressionPosteriorMean(samples, geneIndex-1, expressionCategory[geneIndex])
    i <- i + 1
  }
  expressionValues <- log10(expressionValues)
  #TODO wait for funciton from gabe to subset genome
  
  names.aa <- aminoAcids()
  for(aa in names.aa)
  {
    if(aa == "M" || aa == "W" || aa == "X") next
    plotSinglePanel(parameter, model, genome, expressionValues, samples, category, aa)
    box()
    main.aa <- aa #TODO map to three letter code
    text(0, 1, main.aa, cex = 1.5)
    if(aa %in% c("A", "F", "K", "Q", "V")){
      axis(2, las=1)
    }
    if(aa %in% c("T", "V", "Y", "Z")){
      axis(1)
    }
    if(aa %in% c("A", "C", "D", "E")){
      axis(3)
    }
    if(aa %in% c("E", "I", "P", "T")){
      axis(4, las=1)
    }
    axis(1, tck = 0.02, labels = FALSE)
    axis(2, tck = 0.02, labels = FALSE)
    axis(3, tck = 0.02, labels = FALSE)
    axis(4, tck = 0.02, labels = FALSE)    
  }
  
  ## adding a histogram of phi values to plot
  hist.values <- hist(expressionValues, plot=FALSE, nclass=30)
  plot(hist.values, axes = FALSE, main="", xlab = "", ylab = "")
  axis(1)
  axis(4, las=1)
  
  ### Plot xlab.
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.2, expression("log"[10]~"(Protein Synthesis Rate"~phi~")"))  
  #text(0.5, 0.5, "Production Rate (log10)")
  
  ### Plot ylab.
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.5, "Propotion", srt = 90)
}

plotSinglePanel <- function(parameter, model, genome, expressionValues, samples, category, aa)
{
  codons <- AAToCodon(aa, T)
  
  selection <- vector("numeric", length(codons))
  mutation <- vector("numeric", length(codons))
  for(i in 1:length(codons))
  {
    selection[i] <- parameter$getSelectionPosteriorMeanForCodon(category, samples, codons[i])
    mutation[i] <- parameter$getMutationPosteriorMeanForCodon(category, samples, codons[i])
  }
  
  rank.index <- order(expressionValues)
  codonProbability <- lapply(rank.index,  
                             function(geneIndex){
                               model$CalculateProbabilitiesForCodons(mutation, selection, expressionValues[geneIndex])
                               
                             })
  codons <- AAToCodon(aa, F)
  codonCounts <- vector("list", length(codons))
  for(i in 1:length(codons))
  {
    codonCounts[[i]] <- genome$getCodonCountsPerGene(codons[i])
  }
  codonCounts <- do.call("cbind", codonCounts)
  codonCounts <- codonCounts / rowSums(codonCounts)
  quantiles <- quantile(expressionValues, probs = seq(0.05, 0.95, 0.05), na.rm = T)

  
  
  plot(NULL, NULL, xlim=range(expressionValues, na.rm = T), ylim=c(-0.05,1.05), 
       xlab = "", ylab="", axes = FALSE)
  for(i in 1:length(quantiles))
  {
    if(i == 1){
      tmp.id <- expressionValues < quantiles[i]
    }else if(i == length(quantiles)){
      tmp.id <- expressionValues > quantiles[i]
    }else{
      tmp.id <- expressionValues > quantiles[i] & expressionValues < quantiles[i + 1]
    }
    
    means <- colMeans(codonCounts[tmp.id,], na.rm = T)
    sd <- apply(codonCounts[tmp.id,], 2, sd)
    
    for(k in 1:length(codons))
    {
      points(median(expressionValues[tmp.id]), means[k], col=ribModel:::.codonColors[k], pch=19, cex = 0.5)
      lines(rep(median(expressionValues[tmp.id]),2), c(means[k]-sd[k], means[k]+sd[k]), col=ribModel:::.codonColors[k], lwd=0.8)
    }
  }
  
  codonProbability <- do.call("rbind", codonProbability)
  for(i in 1:length(codons))
  {
    lines(expressionValues[rank.index], codonProbability[, i], col=ribModel:::.codonColors[i])
  }
  legend("topleft", legend = codons, col=ribModel:::.codonColors[1:length(codons)], bty = "n", lty=1, cex=0.75)  
}
