#NOTE!!!
#Expecting file format as follows:
#ORF	Position	Codon	Short	Long
#to be changed to 
#ORF	RFPCount	CodonCount	Codon


table <- read.table("../ribModel/data/GSM1406455_untreated_3.percodon.txt", header=TRUE, sep="\t")
gene.names <- as.character(unique(table[,1])) #Pull all unique gene names from the file


#Set up 4 arrays which will be appended with the information each contains for each gene
ORF <- vector("character") #list of eache gene; gene name dupicated by num of unique codons found in the gene
RFP_Counts <- numeric() #count of "long" (RFP counts) found in the gene
Codon_Counts <- numeric() #codon counts per gene
codon <- numeric() 
rfpVals <- numeric()


#loop through all the gene names
for(gene.name in gene.names)
{
  idx <- table[,1] == gene.name #get which lines of the data frame contain info on the current gene

  #For the current gene, get the list of codons found in the gene (no duplicates) and sort them
  codon.names <- as.character(unique(table[idx,3]))
  codon.names <- sort(codon.names)

  rfpVals <- vector("numeric", length=length(codon.names)) #reset the vector to correct size
  
  #Sum the rfp counts for each codon in the current gene
  i <- 1
  for(codon.name in codon.names) {
    rfpidx <- (as.character(table[idx,3]) == codon.name)
    rfpVals[i] <- sum(table[idx,][rfpidx,5])
    i <- i + 1
  }
  names(rfpVals) <- codon.names

  
  seq.vector <- table[idx, 3]
  rfp.count <- table[idx, 5]
  codon.count <- plyr::count(df=as.data.frame(seq.vector))$freq

  #append the relevant information on the end of the arrays
  ORF <- c(ORF, rep(gene.name, length(codon.names)))
  RFP_Counts <- c(RFP_Counts, rfpVals)
  Codon_Counts <- c(Codon_Counts, codon.count)
  codon <- c(codon, codon.names)
}


newtable <- data.frame(cbind(ORF, RFP_Counts, Codon_Counts, codon), row.names = NULL)
write.table(newtable, file = "GSM1406455_untreated_3.percodon.csv", sep = ",", row.names = FALSE, 
            quote = FALSE)
