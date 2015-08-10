
x <- read.csv("../ribModel/data/SimulatedRFPData.csv")
y <- unique(x[,1])
y <- as.character(y)
names <- sample(y, 1500)

locations <- NULL
for (i in names){
  locations <- c(locations, which(i==x[,1]))
}
z <- x[locations, ]

dim(x[i==x[,1],])
write.csv(z,"1500SimulatedRFPData.csv")
