
mut <- read.table(file="../ribModel/data/genome_2000.mutation.tsv", header=T, sep="\t")
sel <- read.table(file="../ribModel/data/genome_2000.selection.tsv", header=T, sep="\t")

for(aa in as.character(unique(mut[, 1])))
{
  idx <- which(sel[, 1] == aa)
  for(id in idx)
  {
    mut[id, 3:6] <- mut[id, 3:6] - mut[idx[length(idx)], 3:6]
    sel[id, 3:6] <- sel[id, 3:6] - sel[idx[length(idx)], 3:6]
  }
}

mut <- mut[mut[, 3] != 0, ]
sel <- sel[sel[, 3] != 0, ]

write.table(file = "../ribModel/data/genome_2000.mutation.csv", quote = F, sep = ",", x = mut)
write.table(file = "../ribModel/data/genome_2000.selection.csv", quote = F, sep = ",", x = sel)
