for(mixture in 1:parameter$numMixtures)
{
  names.aa <- aminoAcids()
  Amino_Acid <- c()
  Value <- c()
  Codon <- c()
  Std_Deviation <- c()
  for(aa in names.aa)
  {
    if(aa == "M" || aa == "W" || aa == "X") next
    codons <- AAToCodon(aa, T)
    for(i in 1:length(codons))
    {
      Amino_Acid <- c(Amino_Acid, aa)
      Codon <- c(Codon, codons[i])
      Value <- c(Value,parameter$getMutationPosteriorMeanForCodon(mixture, samples*0.1, codons[i]))
      Std_Deviation <- c(Std_Deviation, parameter$getMutationVarianceForCodon(mixture, samples*0.1, codons[i], TRUE))
    }
  }
  data <- data.frame(Amino_Acid,Codon,Value, Std_Deviation)
  filename <- paste("mutationMixture", mixture, ".csv", sep="")
  write.csv(data, file = filename, row.names = FALSE, quote=FALSE)
}

