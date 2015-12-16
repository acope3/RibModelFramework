
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

length.Rcpp_Genome <- function(x)
{
  return(x$getGenomeSize())
}

summary.Rcpp_Genome <- function(object, ...)
{
  # TODO output stuff like:
  # - no. of genes
  # - avg. gene length
  # - avg. A,C,G,T content
  # - avg. AA composition
  # - ...
  summary(object, ...)
}