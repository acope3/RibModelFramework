

library(seqinr)

genome <- read.fasta(file = "../data/realGenomes/Skluyveri.fasta")

main.genes <- (1:5321)[-(1039:1646)]
gene.idx <- sample(x = main.genes, size = 2000, replace = FALSE)
gene.idx <- c(gene.idx, 1039:1646)
gene.idx <- sort(gene.idx)

genome <- genome[gene.idx]

phi.vals <- read.table(file = "../data/realGenomes/Skluyveri_phi.csv", sep=",", header = T)[gene.idx, ]
write.table(file = "../data/realGenomes/Skluyveri_reduced_phi.csv", quote = F, sep = ",", row.names = F, x = phi.vals)
write.fasta(sequences = genome, names = getName(genome), file.out = "../data/realGenomes/Skluyveri_reduced.fasta")

