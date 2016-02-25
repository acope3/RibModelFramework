#' Genome Initialization
#' 
#' \code{initializeGenomeObject} initializes the Rcpp Genome object
#' 
#' @param file A string containing the location of the input file.
#' @param fasta A boolean value which decides whether to initialize with a
#'  fasta file or an RFP value file. (TRUE for fasta, FALSE for RFP)
#' @param expression.file A string containing the location of a file containing
#'  empirical expression rates, if needed.
#' @param append A boolean value that states whether the genes should be appended
#'  to the end of the genome
#' 
#' @return This function returns the Genome Rcpp object created.
#' 
initializeGenomeObject <- function(file, fasta=TRUE, expression.file=NULL, append=FALSE) {
  genome <- new(Genome)
  #genome <- new(Genome, 1, "ROC", TRUE) #CT ONLY
  if (fasta == TRUE) {
    genome$readFasta(file, append)
  } else {
    genome$readRFPFile(file)
  }
  if(!is.null(expression.file)) {
    genome$readObservedPhiValues(expression.file, FALSE)
  }
  return(genome)
}

#' Length of Genome
#' 
#' \code{length} gives the length of a genome
#' 
#' @param genome A genome object initialized with \code{\link{initializeGenomeObject}}.
#' 
#' @return returns the number of genes in a genome
length.Rcpp_Genome <- function(x) {
  return(x$getGenomeSize(F))
}

summary.Rcpp_Genome <- function(genome) {
  # TODO output stuff like:
  # - no. of genes
  # - avg. gene length
  # - avg. A,C,G,T content
  # - avg. AA composition
  # - ...
  summary(object, ...)
}

getNames <- function(genome, simulated = FALSE)
{
  genes <- genome$getGenes(simulated)
  gene.names <- unlist(lapply(1:length(genes), function(i){return(genes[[i]]$id)}))
  return(gene.names)
}
