#' Plot Model Object
#' 
#' @param x An Rcpp model object initialized with \code{initializeModelObject}.
#' @param genome An Rcpp genome object initialized with \code{initializeGenomeObject}.
#' @param samples The number of samples in the trace
#' @param mixture The mixture for which to graph values.
#' @param simulated A boolean value that determines whether to use the simulated genome.
#' @param ... Optional, additional arguments.
#' For this function, a possible title for the plot in the form of a list if set with "main".
#'  
#' @return This function has no return value.
#'  
#' @description Plots traces from the model object such as synthesis rates for each gene.
#' Will work regardless of whether or not expression/synthesis rate levels are being
#' estimated. If you wish to plot observed/empirical values, these values MUST be set
#' using the initial.expression.values parameter found in initializeParameterObject.
#' Otherwise, the expression values plotted will just be SCUO values estimated upon
#' initialization of the Parameter object.

plot.Rcpp_ROCModel <- function(x, genome = NULL, samples = 100, mixture = 1,
                               simulated = FALSE, legacy.layout = FALSE, ...)
{
  model <- x
  opar <- par(no.readonly = T)

  input_list <- as.list(list(...))
  if("main" %in% names(input_list)){
    main <- input_list$main
    input_list$main <- NULL
  }else{
    main <- ""
  }

  num.genes <- length(genome)
  parameter <- model$getParameter()
  mixtureAssignment <- unlist(lapply(1:num.genes, function(geneIndex){
      parameter$getEstimatedMixtureAssignmentForGene(samples, geneIndex)}))
  genes.in.mixture <- which(mixtureAssignment == mixture)
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixture)
  num.genes <- length(genes.in.mixture)
  expressionValues <- unlist(lapply(genes.in.mixture, function(geneIndex){
      parameter$getSynthesisRatePosteriorMeanForGene(samples, geneIndex, FALSE)}))
  expressionValues <- log10(expressionValues)
  genome <- genome$getGenomeForGeneIndices(genes.in.mixture, simulated)
  names.aa <- aminoAcids()

  if(!legacy.layout) {
    # -- compact layout: interleaved AA + marginal columns, n-per-bin strip --
    # Pre-compute bins (shared across all AA panels).
    quantiles  <- quantile(expressionValues, probs = seq(0.05, 0.95, 0.05), na.rm = TRUE)
    n.bins     <- length(quantiles)
    xlimit.global <- range(expressionValues, na.rm = TRUE)
    # Bin counts computed from the first AA panel result; initialise here.
    bin.counts <- NULL
    bin.mids   <- NULL

    lay <- .buildModelLayout(20L)
    layout(lay$mat, lay$widths, lay$heights, respect = FALSE)

    # title
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.6, main)
    text(0.5, 0.4, date(), cex = 0.6)

    for(aa in names.aa) {
      if(aa == "M" || aa == "W" || aa == "X") next
      codon.probability <- calculateProbabilityVector(
          parameter, model, expressionValues, mixture, samples, aa, model.type = "ROC")
      result <- plotSinglePanel(parameter, model, genome, expressionValues,
                                samples, mixture, aa, codon.probability,
                                precomputed.quantiles = quantiles)
      xlimit <- result$xlimit
      if(is.null(bin.counts)) {
        bin.counts <- result$bin.counts
        bin.mids   <- result$bin.mids
      }
      box()
      text(mean(xlimit), 1, aa, cex = 1.5)
      if(aa %in% c("A", "F", "K", "Q", "V")) axis(2, las = 1)
      if(aa %in% c("T", "V", "Y", "Z"))       axis(1)
      if(aa %in% c("A", "C", "D", "E"))       axis(3)
      if(aa %in% c("E", "I", "P", "T"))       axis(4, las = 1)
      axis(1, tck = 0.02, labels = FALSE)
      axis(2, tck = 0.02, labels = FALSE)
      axis(3, tck = 0.02, labels = FALSE)
      axis(4, tck = 0.02, labels = FALSE)
      # marginal ECDF panel (paired with this AA panel)
      .plotCodonECDF(result$codonCounts, result$codons)
    }

    # phi histogram (right column, spans data rows)
    hist.values <- hist(expressionValues, plot = FALSE, nclass = 30)
    par(mar = c(2, 0.2, 0.5, 1.5))
    plot(hist.values, axes = FALSE, main = "", xlab = "", ylab = "")
    axis(1, cex.axis = 0.8); axis(4, las = 1, cex.axis = 0.8)

    # n-per-bin strip (type="h": vertical lines from 0, aligned to bin midpoints)
    par(mar = c(1.5, 2.5, 0.2, 0.5))
    plot(bin.mids, bin.counts, type = "h",
         xlim = xlimit.global, ylim = c(0, max(bin.counts, na.rm = TRUE) * 1.15),
         xlab = "", ylab = "", axes = FALSE, lwd = 3, col = "gray50", lend = 1)
    axis(1, cex.axis = 0.8)
    axis(2, las = 1, cex.axis = 0.7)
    mtext("n", side = 2, las = 1, line = 1.5, cex = 0.7)

    # x-label (xpd=FALSE prevents axis-call residue from bleeding into adjacent panels)
    par(mar = c(0, 0, 0, 0), xpd = FALSE)
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.5, expression("log"[10]~"(Protein Synthesis Rate"~phi~")"))

    # y-label
    par(mar = c(0, 0, 0, 0), xpd = FALSE)
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.5, "Proportion", srt = 90)

  } else {
    # -- legacy layout (original code, unchanged) --
    mat <- matrix(c(rep(1, 4), 2:21, rep(22, 4)),
                  nrow = 7, ncol = 4, byrow = TRUE)
    mat <- cbind(rep(23, 7), mat, rep(24, 7))
    nf <- layout(mat, c(3, rep(8, 4), 2), c(3, 8, 8, 8, 8, 8, 3), respect = FALSE)
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.6, main)
    text(0.5, 0.4, date(), cex = 0.6)

    for(aa in names.aa) {
      if(aa == "M" || aa == "W" || aa == "X") next
      codon.probability <- calculateProbabilityVector(
          parameter, model, expressionValues, mixture, samples, aa, model.type = "ROC")
      result <- plotSinglePanel(parameter, model, genome, expressionValues,
                                samples, mixture, aa, codon.probability)
      xlimit <- result$xlimit
      box()
      text(mean(xlimit), 1, aa, cex = 1.5)
      if(aa %in% c("A", "F", "K", "Q", "V")) axis(2, las = 1)
      if(aa %in% c("T", "V", "Y", "Z"))       axis(1)
      if(aa %in% c("A", "C", "D", "E"))       axis(3)
      if(aa %in% c("E", "I", "P", "T"))       axis(4, las = 1)
      axis(1, tck = 0.02, labels = FALSE)
      axis(2, tck = 0.02, labels = FALSE)
      axis(3, tck = 0.02, labels = FALSE)
      axis(4, tck = 0.02, labels = FALSE)
    }

    hist.values <- hist(expressionValues, plot = FALSE, nclass = 30)
    plot(hist.values, axes = FALSE, main = "", xlab = "", ylab = "")
    axis(1); axis(4, las = 1)
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.2, expression("log"[10]~"(Protein Synthesis Rate"~phi~")"))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.5, "Propotion", srt = 90)
  }

  par(opar)
}

#' Plot Model Object
#' 
#' @param x An Rcpp model object initialized with \code{initializeModelObject}.
#'
#' @param genome An Rcpp genome object initialized with \code{initializeGenomeObject}.
#'
#' @param samples The number of samples in the trace
#'
#' @param mixture The mixture for which to graph values.
#'
#' @param simulated A boolean value that determines whether to use the simulated genome.
#'
#' @param codon.window A boolean value that determines the codon window to use for calculating codon frequencies. If NULL (the default), use complete sequences.
#'
#' @param ... Optional, additional arguments.
#' For this function, a possible title for the plot in the form of a list if set with "main".
#'  
#' @return This function has no return value.
#'
#' @description Plots traces from the model object such as synthesis rates for each gene.
#' Will work regardless of whether or not expression/synthesis rate levels are being
#' estimated. If you wish to plot observed/empirical values, these values MUST be set
#' using the initial.expression.values parameter found in initializeParameterObject.
#' Otherwise, the expression values plotted will just be SCUO values estimated upon
#' initialization of the Parameter object.
#'
plot.Rcpp_FONSEModel <- function(x, genome, samples = 100, mixture = 1,
                               simulated = FALSE, codon.window = NULL,
                               legacy.layout = FALSE, ...)
{
  model <- x
  opar <- par(no.readonly = T)

  input_list <- as.list(list(...))
  if("main" %in% names(input_list)){
    main <- input_list$main
    input_list$main <- NULL
  }else{
    main <- ""
  }

  num.genes <- length(genome)
  parameter <- model$getParameter()
  mixtureAssignment <- unlist(lapply(1:num.genes, function(geneIndex){
      parameter$getEstimatedMixtureAssignmentForGene(samples, geneIndex)}))
  genes.in.mixture <- which(mixtureAssignment == mixture)
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixture)
  num.genes <- length(genes.in.mixture)
  expressionValues <- unlist(lapply(genes.in.mixture, function(geneIndex){
      parameter$getSynthesisRatePosteriorMeanForGene(samples, geneIndex, FALSE)}))
  expressionValues <- log10(expressionValues)
  genome <- genome$getGenomeForGeneIndices(genes.in.mixture, simulated)
  genes <- genome$getGenes(simulated)
  genome$clear()

  if (is.null(codon.window)) {
    codon.window <- seq(1, 100000)
  } else if (length(codon.window) == 2) {
    codon.window <- seq(codon.window[1], codon.window[2])
  }
  for (i in seq_along(genes)) {
    dna   <- genes[[i]]$seq
    start <- seq(1, nchar(dna), 3)
    stop  <- pmin(start + 2, nchar(dna))
    cods  <- substring(dna, start, stop)
    cods  <- cods[codon.window]
    cods  <- cods[!is.na(cods)]
    genes[[i]]$seq <- paste(cods, collapse = '')
    genome$addGene(genes[[i]], simulated)
  }
  names.aa <- aminoAcids()

  if (!legacy.layout) {
    # -- compact layout: interleaved AA + marginal columns, n-per-bin strip --
    quantiles     <- quantile(expressionValues, probs = seq(0.05, 0.95, 0.05), na.rm = TRUE)
    xlimit.global <- range(expressionValues, na.rm = TRUE)
    bin.counts    <- NULL
    bin.mids      <- NULL

    lay <- .buildModelLayout(20L)
    layout(lay$mat, lay$widths, lay$heights, respect = FALSE)

    # title
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.6, main)
    text(0.5, 0.4, date(), cex = 0.6)

    for(aa in names.aa) {
      if(aa == "M" || aa == "W" || aa == "X") next
      codon.probability <- calculateProbabilityVector(
          parameter, model, expressionValues, mixture, samples, aa,
          model.type = "FONSE", codon.window = codon.window)
      result <- plotSinglePanel(parameter, model, genome, expressionValues,
                                samples, mixture, aa, codon.probability,
                                precomputed.quantiles = quantiles)
      xlimit <- result$xlimit
      if(is.null(bin.counts)) {
        bin.counts <- result$bin.counts
        bin.mids   <- result$bin.mids
      }
      box()
      text(mean(xlimit), 1, aa, cex = 1.5)
      if(aa %in% c("A", "F", "K", "Q", "V")) axis(2, las = 1)
      if(aa %in% c("T", "V", "Y", "Z"))       axis(1)
      if(aa %in% c("A", "C", "D", "E"))       axis(3)
      if(aa %in% c("E", "I", "P", "T"))       axis(4, las = 1)
      axis(1, tck = 0.02, labels = FALSE)
      axis(2, tck = 0.02, labels = FALSE)
      axis(3, tck = 0.02, labels = FALSE)
      axis(4, tck = 0.02, labels = FALSE)
      # marginal ECDF panel (paired with this AA panel)
      .plotCodonECDF(result$codonCounts, result$codons)
    }

    # phi histogram (right column)
    hist.values <- hist(expressionValues, plot = FALSE, nclass = 30)
    par(mar = c(2, 0.2, 0.5, 1.5))
    plot(hist.values, axes = FALSE, main = "", xlab = "", ylab = "")
    axis(1, cex.axis = 0.8); axis(4, las = 1, cex.axis = 0.8)

    # n-per-bin strip
    par(mar = c(1.5, 2.5, 0.2, 0.5))
    plot(bin.mids, bin.counts, type = "h",
         xlim = xlimit.global, ylim = c(0, max(bin.counts, na.rm = TRUE) * 1.15),
         xlab = "", ylab = "", axes = FALSE, lwd = 3, col = "gray50", lend = 1)
    axis(1, cex.axis = 0.8)
    axis(2, las = 1, cex.axis = 0.7)
    mtext("n", side = 2, las = 1, line = 1.5, cex = 0.7)

    # x-label (xpd=FALSE prevents axis-call residue from bleeding into adjacent panels)
    par(mar = c(0, 0, 0, 0), xpd = FALSE)
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.5, expression("log"[10]~"(Protein Synthesis Rate"~phi~")"))

    # y-label
    par(mar = c(0, 0, 0, 0), xpd = FALSE)
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.5, "Proportion", srt = 90)

  } else {
    # -- legacy layout (original code, unchanged) --
    mat <- matrix(c(rep(1, 4), 2:21, rep(22, 4)),
                  nrow = 7, ncol = 4, byrow = TRUE)
    mat <- cbind(rep(23, 7), mat, rep(24, 7))
    nf <- layout(mat, c(3, rep(8, 4), 2), c(3, 8, 8, 8, 8, 8, 3), respect = FALSE)
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.6, main)
    text(0.5, 0.4, date(), cex = 0.6)

    for(aa in names.aa) {
      if(aa == "M" || aa == "W" || aa == "X") next
      codon.probability <- calculateProbabilityVector(
          parameter, model, expressionValues, mixture, samples, aa,
          model.type = "FONSE", codon.window = codon.window)
      result <- plotSinglePanel(parameter, model, genome, expressionValues,
                                samples, mixture, aa, codon.probability)
      xlimit <- result$xlimit
      box()
      main.aa <- aa #TODO map to three letter code
      text(mean(xlimit), 1, main.aa, cex = 1.5)
      if(aa %in% c("A", "F", "K", "Q", "V")) axis(2, las = 1)
      if(aa %in% c("T", "V", "Y", "Z"))       axis(1)
      if(aa %in% c("A", "C", "D", "E"))       axis(3)
      if(aa %in% c("E", "I", "P", "T"))       axis(4, las = 1)
      axis(1, tck = 0.02, labels = FALSE)
      axis(2, tck = 0.02, labels = FALSE)
      axis(3, tck = 0.02, labels = FALSE)
      axis(4, tck = 0.02, labels = FALSE)
    }

    hist.values <- hist(expressionValues, plot = FALSE, nclass = 30)
    plot(hist.values, axes = FALSE, main = "", xlab = "", ylab = "")
    axis(1); axis(4, las = 1)
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.2, expression("log"[10]~"(Protein Synthesis Rate"~phi~")"))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.5, "Propotion", srt = 90)
  }

  par(opar)
}

calculateProbabilityVector <- function(parameter,model,expressionValues,mixture,samples,aa,model.type="ROC",codon.window = c(1,300))
{
  codons <- AAToCodon(aa, T)
  
  # get codon specific parameter
  selection <- vector("numeric", length(codons))
  mutation <- vector("numeric", length(codons))
  for (i in 1:length(codons))
  {
    selection[i] <- parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 1, T,log_scale = F)
    mutation[i] <- parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 0, T,log_scale = F)
  }
  
  # calculate codon probabilities with respect to phi
  expression.range <- range(expressionValues)
  phis <- seq(from = expression.range[1], to = expression.range[2], by = 0.01)
  if (model.type == "ROC")
  {
    codonProbability <- lapply(10^phis,  
                               function(phi){
                                 model$CalculateProbabilitiesForCodons(mutation, selection, phi)
                               })
    codonProbability <- do.call("rbind",codonProbability)
  } else if (model.type == "FONSE")
  {
    codonProbability <- vector(mode = "list", length = length(codon.window))
    for (i in 1:length(codon.window))
    {
      codonProbability.tmp <- lapply(10^phis,
                               function(phi)
                                 {
                                 model$CalculateProbabilitiesForCodons(mutation, selection, phi,codon.window[i])
                               })
      codonProbability[[i]] <- do.call("rbind",codonProbability.tmp)
    }
    codonProbability <- Reduce("+",codonProbability)/length(codon.window)
  }
  return(codonProbability)
}



## Internal: build layout matrix for compact codon-freq vs phi plot.
## n.slots: number of AA panel slots to allocate (20 for ROC, covers 19 plotted AAs).
## Panel numbering (in drawing order):
##   1            = title
##   2 .. 2*n+1   = n AA+marg pairs (AA=even offset, marg=odd)
##   2*n+2        = phi histogram (right column)
##   2*n+3        = n-per-bin strip
##   2*n+4        = x-label
##   2*n+5        = y-label (left column)
.buildModelLayout <- function(n.slots = 20L) {
    n.rows.data <- ceiling(n.slots / 4L)
    base        <- 1L                     # title = panel 1
    phi.hist    <- base + 2L * n.slots + 1L
    n.strip     <- base + 2L * n.slots + 2L
    x.lbl       <- base + 2L * n.slots + 3L
    y.lbl       <- base + 2L * n.slots + 4L

    # 8-column data block (4 wide AA cols interleaved with 4 narrow marg cols)
    data.mat <- matrix(0L, nrow = n.rows.data, ncol = 8L)
    pnum <- 2L
    for(r in seq_len(n.rows.data)) {
        for(cp in seq_len(4L)) {
            j <- (cp - 1L) * 2L + 1L
            if(pnum <= base + 2L * n.slots - 1L) {
                data.mat[r, j]      <- pnum
                data.mat[r, j + 1L] <- pnum + 1L
                pnum <- pnum + 2L
            }
        }
    }

    # Assemble full matrix (n.rows.data+3 rows, 10 cols)
    title.row  <- c(y.lbl, rep(1L, 8L),          phi.hist)
    data.rows  <- cbind(y.lbl, data.mat,           phi.hist)
    nstrip.row <- c(y.lbl, rep(n.strip, 8L),       phi.hist)
    xlbl.row   <- c(y.lbl, rep(x.lbl,  8L),        0L)
    mat <- rbind(title.row, data.rows, nstrip.row, xlbl.row)

    list(
        mat      = mat,
        widths   = c(3, 8, 3, 8, 3, 8, 3, 8, 3, 2),
        heights  = c(3, rep(8, n.rows.data), 3, 2),
        phi.hist = phi.hist,
        n.strip  = n.strip,
        x.lbl    = x.lbl,
        y.lbl    = y.lbl
    )
}

## Internal: draw rotated ECDF of per-gene codon proportions.
## Each codon gets one step-function line; x = CDF (0->1), y = proportion.
## y-scale is fixed at (-0.05, 1.05) to match its paired AA panel.
.plotCodonECDF <- function(codonCounts, codons) {
    par(mar = c(0.5, 0.2, 0.5, 1.2))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(-0.05, 1.05),
         axes = FALSE, xlab = "", ylab = "")
    for(k in seq_along(codons)) {
        vals <- sort(codonCounts[, k], na.last = NA)
        n    <- length(vals)
        if(n < 2L) next
        lines(seq_len(n) / n, vals, type = "s",
              col = .codonColors[[ codons[k] ]], lwd = 0.8)
    }
    axis(4, las = 1, tck = 0.04, cex.axis = 0.55)
}

# NOT EXPOSED
# precomputed.quantiles: if provided (compact layout), use instead of recomputing.
# Returns a list(xlimit, codonCounts, codons, bin.mids, bin.counts) so the caller
# can draw the marginal ECDF and n-per-bin strip without re-reading the genome.
plotSinglePanel <- function(parameter, model, genome, expressionValues, samples,
                            mixture, aa, codon.probability,
                            precomputed.quantiles = NULL)
{
  codons <- AAToCodon(aa, T)
  #
  # get codon specific parameter
  selection <- vector("numeric", length(codons))
  mutation <- vector("numeric", length(codons))
  for (i in 1:length(codons))
  {
    selection[i] <- parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 1, T, log_scale = F)
    mutation[i] <- parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 0, T, log_scale = F)
  }
  #
  expression.range <- range(expressionValues)
  phis <- seq(from = expression.range[1], to = expression.range[2], by = 0.01)
  codonProbability <- codon.probability
  #get codon counts
  codons <- AAToCodon(aa, F)
  codonCounts <- vector("list", length(codons))
  for(i in 1:length(codons))
  {
    codonCounts[[i]] <- genome$getCodonCountsPerGene(codons[i])
  }
  codonCounts <- do.call("cbind", codonCounts)
  # codon proportions
  codonCounts <- codonCounts / rowSums(codonCounts)
  codonCounts[is.nan(codonCounts)] <- NA # necessary if AA does not appear in gene

  # make empty plot
  xlimit <- range(expressionValues, na.rm = T)
  par(mar = c(0.5, 2, 0.5, 0.2))
  plot(NULL, NULL, xlim=xlimit, ylim=c(-0.05,1.05),
       xlab = "", ylab="", axes = FALSE)
  # bin expression values of genes
  quantiles <- if(!is.null(precomputed.quantiles)) precomputed.quantiles else
                   quantile(expressionValues, probs = seq(0.05, 0.95, 0.05), na.rm = T)
  n.bins    <- length(quantiles)
  bin.counts <- integer(n.bins)
  bin.mids   <- numeric(n.bins)
  for(i in 1:n.bins)
  {
    if(i == 1){
      tmp.id <- expressionValues < quantiles[i]
    }else if(i == n.bins){
      tmp.id <- expressionValues > quantiles[i]
    }else{
      tmp.id <- expressionValues > quantiles[i] & expressionValues < quantiles[i + 1]
    }
    bin.counts[i] <- sum(tmp.id, na.rm = TRUE)
    bin.mids[i]   <- median(expressionValues[tmp.id], na.rm = TRUE)

    # plot quantiles
    means <- colMeans(codonCounts[tmp.id,], na.rm = T)
    std <- apply(codonCounts[tmp.id,], 2, sd, na.rm = T)
    for(k in 1:length(codons))
    {
      points(bin.mids[i], means[k],
             col=.codonColors[[ codons[k] ]] , pch=19, cex = 0.5)
      lines(rep(bin.mids[i], 2), c(means[k]-std[k], means[k]+std[k]),
            col=.codonColors[[ codons[k] ]], lwd=0.8)
    }
  }

  # draw model fit
  for(i in 1:length(codons))
  {
    lines(phis, codonProbability[, i], col=.codonColors[[ codons[i] ]])
  }
  colors <- unlist(.codonColors[codons])

  # add indicator to optimal codon
  optim.codon.index <- which(min(c(selection, 0)) == c(selection, 0))
  codons[optim.codon.index] <- paste0(codons[optim.codon.index], "*")
  legend("topleft", legend = codons, col=colors, bty = "n", lty=1, cex=0.75)

  invisible(list(xlimit = xlimit, codonCounts = codonCounts,
                 codons = codons, bin.mids = bin.mids, bin.counts = bin.counts))
}

