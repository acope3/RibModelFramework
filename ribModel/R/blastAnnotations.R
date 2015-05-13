source("http://bioconductor.org/biocLite.R")
library(Biostrings)
library(annotate)



genome <- seqinr::read.fasta("~/CodonUsageBias/organisms/bos_taurus/data/bt_mgc_cds_nt.fasta")
estm.phi <- read.csv("~/CodonUsageBias/organisms/bos_taurus/results/Btaurus_thin50_use500/without_xobs_singlechain.phi", as.is=T)
phi.name <- estm.phi[, 1]
estm.phi <- estm.phi[, 2]
names(estm.phi) <- phi.name
estm.phi <- sort(estm.phi, decreasing = T)
phi.name <- names(estm.phi)
numGenes <- round(length(estm.phi)*0.01)


blastset <- vector("character", length=5)
for(j in 1:numGenes)
{
  seq <- seqinr::getSequence.SeqFastadna(genome[[ phi.name[j] ]])
  blastset[j] <- paste(seq, collapse = '')
}

geneDescriptions <- getGeneAnnotationUsingBLAST(blastset, timeout=100)

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



