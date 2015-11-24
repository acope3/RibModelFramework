t <- read.table("../data/FONSE/selection2.csv", sep=",")
ref <- t
for (aa in unique(t[,1])) {
  ref[(t[,1] == aa),3] <- t[(t[,1] == aa),3] - t[max(which(t[,1] == aa)),3]
}

ref <- ref[-which(ref[,3] == 0), ]
write.table(ref, file = "../data/FONSE/selection2ref.csv", sep=",", row.names = F, col.names = F, quote = F)
