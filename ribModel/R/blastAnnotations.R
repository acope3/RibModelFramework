#source("http://bioconductor.org/biocLite.R")
library(Biostrings)
library(annotate)

getGeneAnnotationUsingBLAST <- function(blastset, cut.to=NULL, hitlist=10, timeout=50)
{
  returnSet <- vector("list", length = length(blastset))
  k <- 1
  for(seq in blastset)
  {
    cat("blasting sequence ", k, " of ", length(blastset), "\n")
    resultSet <- NA
    if(!is.null(cut.to))
    {
      if(nchar(seq) > cut.to) seq <- substr(seq, start = 1, stop = cut.to)
    }
    try(resultSet <- blastSequences(x = seq, hitListSize=hitlist, timeout = timeout, as = "data.frame"), program = "megablast")
    try(returnSet[[k]] <- as.character(resultSet$Hit_def))
    
    k <- k + 1
  }
  return(returnSet)
}


genome <- seqinr::read.fasta("~/CodonUsageBias/organisms/bos_taurus/data/bt_mgc_cds_nt.fasta")
estm.phi <- read.csv("~/CodonUsageBias/organisms/bos_taurus/results/Btaurus_thin50_use500/without_xobs_singlechain.phi", as.is=T)
phi.name <- estm.phi[, 1]
estm.phi <- estm.phi[, 2]
names(estm.phi) <- phi.name
estm.phi <- sort(estm.phi, decreasing = T)
phi.name <- names(estm.phi)
numGenes <- round(length(estm.phi)*0.01)


blastset <- vector("character", length=numGenes)
blastnames <- vector("character", length=numGenes)
for(j in 1:numGenes)
{
  seq <- seqinr::getSequence.SeqFastadna(genome[[ phi.name[j] ]])
  blastnames[j] <- phi.name[j]
  blastset[j] <- paste(seq, collapse = '')
}

geneDescriptions <- getGeneAnnotationUsingBLAST(blastset, cut.to=300, timeout=600)
sink("outfile.txt")
for(i in 1:length(blastset))
{
  cat(blastnames[i], "\n")
  cur <- geneDescriptions[[i]]
  for(k in 1:length(cur))
  {
    cat(cur[k], "\n")
  }
  cat("\n")
}
sink()





