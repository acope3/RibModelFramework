table <- read.table("../ribModel/data/GSM1406455_untreated_3.percodon.txt", header=TRUE, sep="\t")
gene.names <- as.character(unique(table[,1]))

ORF <- vector("character")
RFP_Counts <- numeric()
Codon_Counts <- numeric()
codon <- numeric()
rfpVals <- numeric()


for(gene.name in gene.names)
{
  idx <- table[,1] == gene.name
  codon.names <- as.character(unique(table[idx,3]))
  codon.names <- sort(codon.names)
  rfpVals <- vector("numeric", length=length(codon.names))
  start <- Sys.time()
  i <- 1
  for(codon.name in codon.names) {
    rfpidx <- (as.character(table[idx,3]) == codon.name)
    rfpVals[i] <- sum(table[idx,][rfpidx,5])
    i <- i + 1
  }
  names(rfpVals) <- codon.names
  Sys.time() - start
  
  seq.vector <- table[idx, 3]
  rfp.count <- table[idx, 5]
  codon.count <- plyr::count(df=as.data.frame(seq.vector))$freq
  #idx2 <- rfp.count != 0
  #rfp.codons <- seq.vector[idx2]
  #rfp.codon.counts <- rfp.count[idx2]
  #codon.sumation <- unlist(lapply(X = rfp.codons, FUN= function(codon){sum(rfp.codon.counts[rfp.codons == codon])}))
  
  ORF <- c(ORF, rep(gene.name, length(codon.names)))
  RFP_Counts <- c(RFP_Counts, rfpVals)
  Codon_Counts <- c(Codon_Counts, codon.count)
  codon <- c(codon, codon.names)
}


newtable <- data.frame(cbind(ORF, RFP_Counts, Codon_Counts, codon), row.names = NULL)
write.table(newtable, file = "GSM1406455_untreated_3.percodon.csv", sep = ",", row.names = FALSE, 
            quote = FALSE)
