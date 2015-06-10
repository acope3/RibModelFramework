
initializeGenomeObject <- function(fasta.file, expression.file=NULL, append=FALSE)
{
  genome <- new(Genome)
  genome$readFasta(fasta.file, append)
  if(!is.null(expression.file))
  {
    #TODO implement reading expressin csv file and attaching expression values to genes
  }
  return(genome)
}