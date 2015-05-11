source("http://bioconductor.org/biocLite.R")
library(Biostrings)
library(annotate)



genome <- seqinr::read.fasta("~/CodonUsageBias/organisms/bos_taurus/data/bt_mgc_cds_nt.fasta")

blastset <- vector("character", length=5)
i <- 1
for(j in 10:15)
{
  blastset[i] <- paste(seqinr::getSequence.SeqFastadna(genome[[j]]), collapse = '')
  i <- i + 1
}

getGeneAnnotationUsingBLAST <- function(blastset, hitlist=10, timeout=50)
{
  returnSet <- vector("list", length = length(blastset))
  k <- 1
  for(seq in blastset)
  {
    cat("blasting sequence ", k, " of ", length(blastset), "\n")
    resultSet <- NA
    try(resultSet <- blastSequences(x = seq, hitListSize=hitlist, timeout = timeout, as = "data.frame"))
    returnSet[[k]] <- as.character(resultSet$Hit_def)
    k <- k + 1
  }
  return(returnSet)
}



