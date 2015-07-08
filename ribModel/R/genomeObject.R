
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
    #TODO implement reading expressin csv file and attaching expression values to genes
  }
  return(genome)
}