
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