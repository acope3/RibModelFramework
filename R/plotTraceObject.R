# Plot functions for trace object
# The generic plot function expects the trace object
# and a string to the the function what has to be ploted.
# additional arguments are geneIndex, and category to index function like
# getExpressionTraceForGene or plotCodonSpecificParameters

#' Plot Trace Object
#' @param x An Rcpp trace object initialized with \code{initializeTraceObject}.
#' @param what A string containing one of the following to graph: \code{Mutation, Selection, Alpha, LambdaPrime, MeanWaitingTime, VarWaitingTime
#' MixtureProbability, Sphi, Mphi, Aphi, Spesilon, ExpectedPhi, Expression}.
#' @param geneIndex When plotting expression, the index of the gene to be plotted.
#' @param mixture The mixture for which to plot values.
#' @param log.10.scale A logical value determining if figures should be plotted on the log.10.scale (default=F). Should not be applied to mutation and selection parameters estimated by ROC/FONSE.
#' @param aa.names A vector of single letter amino acid names used to set order of plotting
#'
#' @param ... Optional, additional arguments.
#' For this function, may be a logical value determining if the trace is ROC-based or not.
#'
#' @return This function has no return value.
#' 
#' @description Plots different traces, specified with the \code{what} parameter.
#'
plot.Rcpp_Trace <- function(x, what=c("Mutation", "Selection", "MixtureProbability" ,"Sphi", "Mphi", "Aphi", "Sepsilon", "ExpectedPhi", "Expression","NSEProb","NSERate","InitiationCost","PartitionFunction"),
                            geneIndex=1,
                            mixture = 1,
                            log.10.scale=F,
                            aa.names = aminoAcids(),
                            legacy.layout = FALSE,
                            ...
                            )
{
  if(what[1] == "Mutation")
  {
    plotCodonSpecificParameters(x, mixture, "Mutation", main="Mutation Parameter Traces", aa.names = aa.names, legacy.layout = legacy.layout)
  }
  if(what[1] == "Selection")
  {
    plotCodonSpecificParameters(x, mixture, "Selection", main="Selection Parameter Traces", aa.names = aa.names, legacy.layout = legacy.layout)
  }
  if(what[1] == "Alpha")
  {
    plotCodonSpecificParameters(x, mixture, "Alpha", main="Alpha Parameter Traces", ROC.or.FONSE=FALSE, log.10.scale=log.10.scale, aa.names = aa.names, legacy.layout = legacy.layout)
  }
  if(what[1] == "Lambda")
  {
    plotCodonSpecificParameters(x, mixture, "Lambda", main="Lambda Parameter Traces", ROC.or.FONSE=FALSE, log.10.scale=log.10.scale, aa.names = aa.names, legacy.layout = legacy.layout)
  }

  if(what[1] == "MeanWaitingTime")
  {
    plotCodonSpecificParameters(x, mixture, "MeanWaitingTime", main="Mean Waiting Time Parameter Traces", ROC.or.FONSE=FALSE, log.10.scale=log.10.scale, aa.names = aa.names, legacy.layout = legacy.layout)
  }
  if(what[1] == "VarWaitingTime")
  {
    plotCodonSpecificParameters(x, mixture, "VarWaitingTime", main="Variance Waiting Time Parameter Traces", ROC.or.FONSE=FALSE, aa.names = aa.names, legacy.layout = legacy.layout)
  }
  if(what[1] == "NSEProb")
  {
    plotCodonSpecificParameters(x, mixture, "NSEProb", main="Nonsense Error Probability Parameter Traces", ROC.or.FONSE=FALSE, log.10.scale=log.10.scale, aa.names = aa.names, legacy.layout = legacy.layout)
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
  if(what[1] == "InitiationCost")
  {
    plotFONSEHyperParameterTrace(x,what=what[1])
  }
  if(what[1] == "PartitionFunction")
  {
    plotPANSEHyperParameterTrace(x,what=what[1])
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
  if(what[1] == "AcceptanceRatio")
  {
    plotAcceptanceRatios(x)
  }
  if(what[1] == "NSERate")
  {
    plotCodonSpecificParameters(x, mixture, "NSERate", main="NSERate", ROC.or.FONSE=FALSE, log.10.scale=log.10.scale, aa.names = aa.names, legacy.layout = legacy.layout)
  }
}

.buildTraceLayout <- function(n.slots = 20L) {
    n.rows.data <- ceiling(n.slots / 4L)
    x.lbl <- 2L * n.slots + 2L
    y.lbl <- 2L * n.slots + 3L

    data.mat <- matrix(0L, nrow = n.rows.data, ncol = 8L)
    pnum <- 2L
    for (r in seq_len(n.rows.data)) {
        for (cp in seq_len(4L)) {
            j <- (cp - 1L) * 2L + 1L
            if (pnum <= 2L * n.slots) {
                data.mat[r, j]      <- pnum
                data.mat[r, j + 1L] <- pnum + 1L
                pnum <- pnum + 2L
            }
        }
    }
    title.row <- c(y.lbl, rep(1L, 8L))
    data.rows <- cbind(y.lbl, data.mat)
    xlbl.row  <- c(y.lbl, rep(x.lbl, 8L))
    mat <- rbind(title.row, data.rows, xlbl.row)
    list(mat     = mat,
         widths  = c(3, 8, 3, 8, 3, 8, 3, 8, 3),
         heights = c(3, rep(8, n.rows.data), 2),
         x.lbl   = x.lbl,
         y.lbl   = y.lbl)
}

.plotTraceMarginal <- function(cur.trace, codons, ylim) {
    par(mar = c(0.5, 0.2, 0.5, 1.2))
    plot(NULL, NULL, xlim = c(0, 1), ylim = ylim, axes = FALSE, xlab = "", ylab = "")
    for (k in seq_len(ncol(cur.trace))) {
        vals <- sort(cur.trace[, k], na.last = NA)
        n    <- length(vals)
        if (n < 2L) next
        lines(seq_len(n) / n, vals, type = "s",
              col = .codonColors[[codons[k]]], lwd = 0.8)
    }
    axis(4, las = 1, tck = 0.04, cex.axis = 0.55)
}

# Called from Plot Trace Object (plot for trace)
# NOT EXPOSED
#
#' Plot Codon Specific Parameter
#' @param trace An Rcpp trace object initialized with \code{initializeTraceObject}.
#'
#' @param mixture The mixture for which to plot values.
#'
#' @param type A string containing one of the following to graph: \code{Mutation, Selection, Alpha, LambdaPrime, MeanWaitingTime, VarWaitingTime}. 
#'
#' @param main The title of the plot.
#'
#' @param ROC.or.FONSE A logical value determining if the Parameter was ROC/FONSE or not.
#'
#' @param log.10.scale A logical value determining if figures should be plotted on the log.10.scale (default=F). Should not be applied to mutation and selection parameters estimated by ROC/FONSE.
#'
#' @param aa.names A character vector of amino acid names to include in the plot. Default is all amino acids from \code{aminoAcids()}.
#'
#' @return This function has no return value.
#'
#' @description Plots a codon-specific set of traces, specified with the \code{type} parameter.
#'
plotCodonSpecificParameters <- function(trace, mixture, type="Mutation", main="Mutation Parameter Traces", ROC.or.FONSE=TRUE, log.10.scale=F, aa.names = aminoAcids(), legacy.layout=FALSE)
{
  opar <- par(no.readonly = T)

  ### TODO change to groupList -> checks for ROC like model is not necessary!

  ## Check to ensure aa.names passed are valid
  if( is.null(aa.names) ) {
      aa.names <- aminoAcids()
  } else {
      aa.match <- (aa.names %in% aminoAcids())
      ## test to ensure there's no aa being called that don't exist in trace
      aa.mismatch <- aa.names[!aa.match]
      if(length(aa.mismatch) > 0){
          warning("Members ", aa.mismatch, "of aa.names argument absent from trace object and will be excluded.",
                  call. = TRUE, immediate. = FALSE, noBreaks. = FALSE,
                  domain = NULL)
          }
      aa.names <- aa.names[aa.match]
  }
  with.ref.codon <- ifelse(ROC.or.FONSE, TRUE, FALSE)

  ## Determine ylab and extraction parameters once (shared by both layout branches)
  if(type == "Mutation"){
    ylab <- expression(Delta~"M")
    paramType <- 0
    special <- FALSE
  }else if (type == "Selection"){
    ylab <- expression(Delta~eta)
    paramType <- 1
    special <- FALSE
  }else if (type == "Alpha"){
    ylab <- if (log.10.scale) expression("log"[10]*alpha) else expression(alpha)
    paramType <- 0
    special <- FALSE
  }else if (type == "Lambda"){
    ylab <- if (log.10.scale) expression("log"[10]*lambda) else expression(lambda)
    paramType <- 1
    special <- FALSE
  }else if (type == "MeanWaitingTime"){
    ylab <- if (log.10.scale) expression("log"[10]*alpha/lambda) else expression(alpha/lambda)
    special <- TRUE
  }else if (type == "VarWaitingTime"){
    ylab <- expression(alpha/lambda^"2")
    special <- TRUE
  }else if (type == "NSEProb"){
    ylab <- if (log.10.scale) expression("log"[10]*"Pr(NSE)") else expression("E[Pr(NSE)]")
    special <- TRUE
  }else if (type == "VarNSEProb"){
    ylab <- expression("Var[Pr(NSE)]")
    special <- TRUE
  }else if (type == "NSERate"){
    ylab <- if (log.10.scale) expression("log"[10]*"NSERate") else expression("NSERate")
    paramType <- 2
    special <- FALSE
  }else{
    stop("Parameter 'type' not recognized! Must be one of: 'Mutation', 'Selection', 'Alpha', 'Lambda', 'MeanWaitingTime', 'VarWaitingTime', 'NSEProb', 'NSERate'.")
  }

  ## Extract MCMC trace for one amino acid as a matrix (one column per codon)
  .extractTrace <- function(aa) {
    codons <- AAToCodon(aa, with.ref.codon)
    tr <- vector("list", length(codons))
    for (i in seq_along(codons)) {
      if (special) {
        tmpAlpha       <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], 0, with.ref.codon)
        tmpLambdaPrime <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], 1, with.ref.codon)
        if (type == "MeanWaitingTime") {
          tr[[i]] <- tmpAlpha / tmpLambdaPrime
        } else if (type == "VarWaitingTime") {
          tr[[i]] <- tmpAlpha / (tmpLambdaPrime * tmpLambdaPrime)
        } else if (type %in% c("NSEProb", "VarNSEProb")) {
          tmpNSERate <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], 2, with.ref.codon)
          tr[[i]] <- if (type == "NSEProb")
            tmpNSERate * (tmpAlpha / tmpLambdaPrime)
          else
            tmpNSERate * tmpNSERate * (tmpAlpha / (tmpLambdaPrime * tmpLambdaPrime))
        }
        if (log.10.scale) tr[[i]] <- log10(tr[[i]])
      } else {
        tr[[i]] <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], paramType, with.ref.codon)
        if (log.10.scale) tr[[i]] <- log10(tr[[i]])
      }
    }
    do.call("cbind", tr)
  }

  if (!legacy.layout) {
    ## --- compact layout branch ---

    ## Pre-compute valid AAs so n.slots is known before layout() is called
    valid.aas <- Filter(function(aa) {
      cods <- AAToCodon(aa, with.ref.codon)
      if (length(cods) == 0) return(FALSE)
      if (ROC.or.FONSE && aa %in% c("X", "M", "W")) return(FALSE)
      if (!ROC.or.FONSE && aa == "X") return(FALSE)
      TRUE
    }, aa.names)

    if (length(valid.aas) == 0) { par(opar); return(invisible(NULL)) }

    lo <- .buildTraceLayout(n.slots = length(valid.aas))
    layout(lo$mat, widths = lo$widths, heights = lo$heights, respect = FALSE)

    ## panel 1: title
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.6, main)
    text(0.5, 0.4, date(), cex = 0.6)

    for (aa in valid.aas) {
      codons    <- AAToCodon(aa, with.ref.codon)
      cur.trace <- .extractTrace(aa)
      if (length(cur.trace) == 0) next

      x    <- seq_len(nrow(cur.trace))
      xlim <- range(x)
      ylim <- range(cur.trace, na.rm = TRUE)

      ## trace panel (panels 2, 4, 6, ... in layout order)
      par(mar = c(0.5, 2.0, 1.5, 0.2))
      plot(NULL, NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes = FALSE)
      plot.order <- order(apply(cur.trace, 2, sd), decreasing = TRUE)
      for (i.codon in plot.order) {
        lines(x = x, y = cur.trace[, i.codon], col = .codonColors[[codons[i.codon]]])
      }
      axis(2, las = 1, cex.axis = 0.6)
      axis(1, tck = 0.02, labels = FALSE)
      colors <- unlist(.codonColors[codons])
      legend("topleft", legend = codons, col = colors,
             lty = rep(1, length(codons)), bty = "n", cex = 0.6)
      mtext(aa, side = 3, line = 0.2, cex = 0.9)

      ## marginal ECDF panel (panels 3, 5, 7, ... in layout order)
      .plotTraceMarginal(cur.trace, codons, ylim)
    }

    ## x-label panel (xpd=FALSE prevents axis-call residue from bleeding)
    par(mar = c(0, 0, 0, 0), xpd = FALSE)
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.5, "Samples")

    ## y-label panel
    par(mar = c(0, 0, 0, 0), xpd = FALSE)
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.5, ylab, srt = 90)

  } else {
    ## --- legacy layout branch (original behavior) ---
    if (ROC.or.FONSE) {
      nf <- layout(matrix(c(rep(1, 4), 2:21), nrow = 6, ncol = 4, byrow = TRUE),
                   rep(1, 4), c(2, 8, 8, 8, 8, 8), respect = FALSE)
    } else {
      nf <- layout(matrix(c(rep(1, 4), 2:25), nrow = 7, ncol = 4, byrow = TRUE),
                   rep(1, 4), c(2, 8, 8, 8, 8, 8, 8), respect = FALSE)
    }
    if (ROC.or.FONSE) {
      par(mar = c(0, 0, 0, 0))
    } else {
      par(mar = c(1, 1, 1, 1))
    }
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.6, main)
    text(0.5, 0.4, date(), cex = 0.6)
    par(mar = c(5.1, 4.1, 4.1, 2.1))

    for (aa in aa.names) {
      codons <- AAToCodon(aa, with.ref.codon)
      if (length(codons) == 0) next
      if (!ROC.or.FONSE && aa == "X") next
      if (ROC.or.FONSE && aa %in% c("X", "M", "W")) next
      cur.trace <- .extractTrace(aa)
      if (length(cur.trace) == 0) next

      x    <- seq_len(nrow(cur.trace))
      xlim <- range(x)
      ylim <- range(cur.trace, na.rm = TRUE)

      main.aa <- aa #TODO map to three letter code
      plot(NULL, NULL, xlim = xlim, ylim = ylim,
           xlab = "Samples", ylab = ylab, main = main.aa)
      plot.order <- order(apply(cur.trace, 2, sd), decreasing = TRUE)
      for (i.codon in plot.order) {
        lines(x = x, y = cur.trace[, i.codon], col = .codonColors[[codons[i.codon]]])
      }
      colors <- unlist(.codonColors[codons])
      legend("topleft", legend = codons, col = colors,
             lty = rep(1, length(codons)), bty = "n", cex = 0.75)
    }
  }

  par(opar)
}

# Called from Plot Trace Object (plot for trace)
# NOT EXPOSED
# 
#' Plot Acceptance ratios
#' @param trace An Rcpp trace object initialized with \code{initializeTraceObject}.
#'
#' @param main The title of the plot.
#'
#' @param aa.names A character vector of amino acid names to include in the plot. Default is all amino acids from \code{aminoAcids()}.
#'
#' @return This function has no return value.
#'
#' @description Plots acceptance ratios for codon-specific parameters. Will be by amino acid for ROC and FONSE models, but will be by codon for PA and PANSE models. Note assumes estimating parameters for all codons.

plotAcceptanceRatios <- function(trace,
                                 main="CSP Acceptance Ratio Traces",
                                 aa.names = aminoAcids())
{
  opar <- par(no.readonly = T) 
  
  ### Trace plot.
  acceptance.rate.traces <- trace$getCodonSpecificAcceptanceRateTrace()
  if (length(acceptance.rate.traces) == 61)
  {
    ROC.or.FONSE <- FALSE  
  } else {
    ROC.or.FONSE <- TRUE
  }



  if (ROC.or.FONSE)
  {
    nf <- layout(matrix(c(rep(1, 4), 2:21), nrow = 6, ncol = 4, byrow = TRUE),
               rep(1, 4), c(2, 8, 8, 8, 8, 8), respect = FALSE)  
  }else
  {    
    nf <- layout(matrix(c(rep(1, 4), 2:25), nrow = 7, ncol = 4, byrow = TRUE),
                    rep(1, 4), c(2, 8, 8, 8, 8, 8, 8), respect = FALSE) 
  }
  ### Plot title.
  if (ROC.or.FONSE){
    par(mar = c(0, 0, 0, 0))
  }else{
    par(mar = c(1,1,1,1))
  }
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, main)
  text(0.5, 0.4, date(), cex = 0.6)
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  
  # TODO change to groupList -> checks for ROC like model is not necessary!
  with.ref.codon <- ifelse(ROC.or.FONSE, TRUE, FALSE)
  
  for(aa in aa.names)
  { 
    if (!ROC.or.FONSE){
      if(aa == "X") next
    } else if(ROC.or.FONSE){
        if(aa == "X" || aa == "W" || aa == "M") next
    }
    
    codons <- AAToCodon(aa, with.ref.codon)
    if(length(codons) == 0) next

    if (!ROC.or.FONSE)
    {
      cur.trace <- vector("list", length(codons))
      for(i in 1:length(codons))
      { 
        cur.trace[[i]] <- trace$getCodonSpecificAcceptanceRateTraceForCodon(codons[i])
      }
    } else {
      cur.trace <- vector("list", 1)
      cur.trace[[1]] <- trace$getCodonSpecificAcceptanceRateTraceForAA(aa)
      }
    cur.trace <- do.call("cbind", cur.trace)
    if(length(cur.trace) == 0) next
    x <- 1:dim(cur.trace)[1]
    xlim <- range(x)
    ylim <- range(cur.trace, na.rm=T)
    
    main.aa <- aa #TODO map to three leter code
    plot(NULL, NULL, xlim = xlim, ylim = ylim,
         xlab = "Samples", ylab = "Accept. Rat.", main = main.aa)
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
# Currently can only be one of Sphi, Mphi, Aphi, and Sepsilon.
plotHyperParameterTrace <- function(trace, what = c("Sphi", "Mphi", "Aphi", "Sepsilon"))
{
#  opar <- par(no.readonly = T) 
#  par(oma=c(1,1,2,1), mgp=c(2,1,0), mar = c(3,4,2,1), mfrow=c(2, 1))
  xlab <- "Sample"
  
  if (what[1] == "Sphi")
  {
    sphi <- trace$getStdDevSynthesisRateTraces();
    numMixtures <- length(sphi)
    sphi <- do.call("cbind", sphi)
    
    ylimit <- range(sphi) + c(-0.1, 0.1)
    xlimit <- c(1, nrow(sphi))
    ylab <- expression("s"[phi])
    main <- expression("s"[phi]*"Trace")
    plot(NULL, NULL, type="l", xlab = xlab, ylab = ylab, xlim = xlimit, ylim = ylimit, main = main)
    
    for(i in 1:ncol(sphi))
    {
      lines(sphi[-1,i], col = .mixtureColors[i])
    }
    legend("topleft", legend = paste0("Mixture Element", 1:numMixtures), 
           col = .mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")
  }
  else if (what[1] == "Mphi")
  {
    sphi <- trace$getStdDevSynthesisRateTraces();
    numMixtures <- length(sphi)
    sphi <- do.call("cbind", sphi)
    mphi <- -(sphi * sphi) / 2;

    ylimit <- range(mphi) + c(-0.1, 0.1)
    xlimit <- c(1, nrow(mphi))
    ylab <- expression("m"[phi])
    main <- expression("m"[phi]*"Trace")
    plot(NULL, NULL, type="l", xlab = xlab, ylab = ylab, xlim  = xlimit, ylim = ylimit, main = main)
    
    for(i in 1:ncol(mphi))
    {
      lines(mphi[-1,i], col= .mixtureColors[i])
    }
    legend("topleft", legend = paste0("Mixture Element", 1:numMixtures), 
           col = .mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")    

  }
  else if (what[1] == "Aphi") 
  {
    aphi <- trace$getSynthesisOffsetTrace();
    aphi <- do.call("cbind", aphi)
    
    ylimit <- range(aphi) + c(-0.1, 0.1)
    xlimit <- c(1, nrow(aphi))
    ylab <- expression("A"[phi])
    main <- expression("A"[phi]~"Trace")
    plot(NULL, NULL, type="l", xlab = xlab, ylab = ylab, xlim  = xlimit, ylim = ylimit, main = main)
    
    num.obs.data <- ncol(aphi)
    for(i in 1:num.obs.data)
    {
      lines(aphi[-1,i], col = .mixtureColors[i])
    }
    legend("topleft", legend = paste0("Observed Data", 1:num.obs.data), 
           col = .mixtureColors[1:num.obs.data], lty = rep(1, num.obs.data), bty = "n")        
  }
  else if (what[1] == "Sepsilon")
  {
    sepsilon <- trace$getObservedSynthesisNoiseTrace();
    sepsilon <- do.call("cbind", sepsilon)

    ylimit <- range(sepsilon) + c(-0.1, 0.1)
    xlimit <- c(1, nrow(sepsilon))
    ylab <- expression("s"[epsilon])
    main <- expression("s"[epsilon]~"Trace")
    plot(NULL, NULL, type="l", xlab = xlab, ylab = ylab, xlim  = xlimit, ylim = ylimit, main = main)

    num.obs.data <- ncol(sepsilon)
    for(i in 1:num.obs.data)
    {
      lines(sepsilon[-1,i], col = .mixtureColors[i])
    }
    legend("topleft", legend = paste0("Observed Data", 1:num.obs.data), 
           col = .mixtureColors[1:num.obs.data], lty = rep(1, num.obs.data), bty = "n")  
  }
  #par(opar)
}


plotFONSEHyperParameterTrace <- function(trace, what = c("InitiationCost"))
{
#  opar <- par(no.readonly = T) 
#  par(oma=c(1,1,2,1), mgp=c(2,1,0), mar = c(3,4,2,1), mfrow=c(2, 1))
  xlab <- "Sample"
  
  if (what[1] == "InitiationCost")
  {
    a1 <- unlist(trace$getInitiationCostTrace())
    a1 <- a1[2:length(a1)]
    ylimit <- range(a1) + c(-0.1, 0.1)
    xlimit <- c(1, length(a1))
    ylab <- expression("a"[1])
    main <- expression("a"[1]*"Trace")
    plot(NULL, NULL, type="l", xlab = xlab, ylab = ylab, xlim = xlimit, ylim = ylimit, main = main)
    
    lines(a1, col = "black")
  }
  
  #par(opar)
}


plotPANSEHyperParameterTrace <- function(trace, what = c("PartitionFunction"))
{
#  opar <- par(no.readonly = T) 
#  par(oma=c(1,1,2,1), mgp=c(2,1,0), mar = c(3,4,2,1), mfrow=c(2, 1))
  xlab <- "Sample"
  
  if (what[1] == "PartitionFunction")
  {
    pf <- trace$getPartitionFunctionTraces();
    numMixtures <- length(pf)
    pf <- do.call("cbind", pf)
    
    ylimit <- range(pf) + c(-0.1, 0.1)
    xlimit <- c(1, nrow(pf))
    ylab <- expression("Partition Function")
    main <- expression("Partition Function Trace")
    plot(NULL, NULL, type="l", xlab = xlab, ylab = ylab, xlim = xlimit, ylim = ylimit, main = main)
    
    for(i in 1:ncol(pf))
    {
      lines(pf[-1,i], col = .mixtureColors[i])
    }
    legend("topleft", legend = paste0("Mixture Element", 1:numMixtures), 
           col = .mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")
  }
  
  #par(opar)
}





# NOT EXPOSED
plotMixtureProbability <- function(trace)
{
  samples <- length(trace$getMixtureProbabilitiesTraceForMixture(1))
  numMixtures <- trace$getNumberOfMixtures()
  
  plot(NULL, NULL, xlim = c(0, samples), ylim=c(0, 1), xlab = "Samples", ylab = "Mixture Probability", main = "Mixture Probability")
  for (i in 1:numMixtures)
  {
    lines(trace$getMixtureProbabilitiesTraceForMixture(i)[-1], col = .mixtureColors[i])    
  }
  legend("topleft", legend = paste0("Mixture Element", 1:numMixtures), 
         col = .mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")
}
