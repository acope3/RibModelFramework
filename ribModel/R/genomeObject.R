#' Genome initialization
#' 
#' \code{initializeGenomeObject} initializes the RCpp Genome object
#' 
#' @param file A string containing the location of the input file.
#' @param fasta A boolean value which decides whether to initialize with a
#'  fasta file or an RFP value file. (TRUE for fasta, FALSE for RFP)
#' @param expression.file A string containing the location of a file containing
#'  empirical expression rates, if needed.
#' @param append A boolean value that states whether the genes should be appended
#'  to the end of the genome
#'
#' @return This function has no return value.
#' 
initializeGenomeObject <- function(file, fasta=TRUE, expression.file=NULL, append=FALSE)
{
  genome <- new(Genome)
  if (fasta == TRUE)
  {
    genome$readFasta(file, append)
  }
  else
  {
    genome$readRFPFile(file)
  }
  if(!is.null(expression.file))
  {
    genome$readObservedPhiValues(expression.file, FALSE)
  }
  return(genome)
}

length.Rcpp_Genome <- function(genome)
{
  return(genome$getGenomeSize())
}

summary.Rcpp_Genome <- function(genome)
{
  # TODO output stuff like:
  # - no. of genes
  # - avg. gene length
  # - avg. A,C,G,T content
  # - avg. AA composition
  # - ...
}