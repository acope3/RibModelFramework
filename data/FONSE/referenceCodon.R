t <- read.table("selection1.csv", sep=",")
AA.list <- unique(t[,1])
for (aa in AA.list) {
  idx = (t[,1] == aa)
  m <- max(which(idx == TRUE))
  
}