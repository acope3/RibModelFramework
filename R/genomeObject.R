#' Genome Initialization
#' 
#' \code{initializeGenomeObject} initializes the Rcpp Genome object
#' 
#' @param file A string containing the location of the input file.
#' @param fasta A boolean value which decides whether to initialize with a
#'  fasta file or an RFP value file. (TRUE for fasta, FALSE for RFP)
#' @param observed.expression.file A string containing the location of a file containing
#'  empirical expression rates, if needed.
#' @param append A boolean value that states whether the genes should be appended
#'  to the end of the genome
#' @param match.expression.by.id If TRUE (default) observed expression values will be assigned by matching sequence identifier.
#' If FALSE observed expression values will be assigned by order
#' 
#' @return This function returns the Genome Rcpp object created.
#' 
initializeGenomeObject <- function(file, observed.expression.file=NULL, append=FALSE, fasta=TRUE, match.expression.by.id=TRUE) {
  genome <- new(Genome)
  #genome <- new(Genome, 1, "ROC", TRUE) #CT ONLY
  if (fasta == TRUE) {
    genome$readFasta(file, append)
  } else {
    genome$readRFPFile(file)
  }
  if(!is.null(observed.expression.file)) {
    genome$readObservedPhiValues(observed.expression.file, match.expression.by.id)
  }
  return(genome)
}

#' Get Codon Counts For Each Amino Acid 
#' @param genome A genome object from which the counts of each
#' codon can be obtained.
#'  
#' @param aa A one character representation of an amino acid.
#' 
#' @return codonCounts Returns a matrix storing the codonCounts. 
#' 
#' @description \code{getCodonCountsForAA} returns a matrix filled with 
#' the number of times a codon is seen in each gene.
#' 
#' @details The returned matrix will have the row correspond to the
#' genes in the genome and the columns correspond to the codons for the 
#' given aa. The values will the number of times the codon is present in 
#' that gene.
#' 
getCodonCountsForAA <- function(genome, aa){
  # get codon count for aa
  codons <- AAToCodon(aa, F)
  codonCounts <- lapply(codons, function(codon){
    codonCounts <- genome$getCodonCountsPerGene(codon)
  })
  codonCounts <- do.call("cbind", codonCounts)
  return(codonCounts)
}



#' Length of Genome
#' 
#' \code{length} gives the length of a genome
#' 
#' @param x A genome object initialized with \code{\link{initializeGenomeObject}}.
#' 
#' @return returns the number of genes in a genome
length.Rcpp_Genome <- function(x) {
  return(x$getGenomeSize(F))
}

#' Summary of Genome
#' 
#' \code{summary} summarizes the description of a genome, such as number of genes and average gene length.
#' 
#' @param object A genome object initialized with \code{\link{initializeGenomeObject}}.
#' 
#' @param ... Optional, additional arguments to be passed to the main summary function 
#' that affect the summary produced.
#'
#' @return This function returns by default an object of class c("summaryDefault", table").
summary.Rcpp_Genome <- function(object, ...) {
  # TODO output stuff like:
  # - no. of genes
  # - avg. gene length
  # - avg. A,C,G,T content
  # - avg. AA composition
  # - ...
  summary(object, ...)
}

#' Gene Names of Genome
#' 
#' \code{getNames} prints the names of the genes within the genome specified.
#' 
#' @param genome A genome object initialized with \code{\link{initializeGenomeObject}}.
#' 
#' @param simulated A logical value denoting if the gene names to be listed are simulated or not.
#' The default value is FALSE.
#' 
#' @return gene.names Returns the names of the genes as a vector of strings.
#' 
getNames <- function(genome, simulated = FALSE)
{
  genes <- genome$getGenes(simulated)
  gene.names <- unlist(lapply(1:length(genes), function(i){return(genes[[i]]$id)}))
  return(gene.names)
}

getObservedSynthesisRateSet <- function(genome, simulated = FALSE)
{
  genes <- genome$getGenes(simulated)
  expression <- lapply(1:length(genes), function(i){return(genes[[i]]$getObservedSynthesisRateValues())})
  
  return(do.call(rbind, expression))
}


