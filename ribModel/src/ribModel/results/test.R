test.phi <- read.table("phiPosterior.csv", sep=",")[,2]
test.phi.names <- as.character(read.table("phiPosterior.csv", sep=",")[,1])
true.phi <- read.table("../SimulatedGenome_allUnique_phi.csv", sep=",", header=T)[, 2]


idx <- 1:500
plot(log10(true.phi[idx]), log10(test.phi[idx]), xlim = range(log10(true.phi)), ylim = range(log10(test.phi)))
points(log10(true.phi[-idx]), log10(test.phi[-idx]), col = "red")
legend("topleft", legend = c("Category 0", "Category 1"), col = c("black", "red"), pch = c(1,1))
abline(0,1, col="blue", lwd=2)
cor(log10(true.phi), log10(test.phi))

mutation <- read.table("mutationPosterior_Cat0.csv", sep=",")[,2]
mutation <- mutation[mutation != 0]
selection <- read.table("selectionPosterior_Cat0.csv", sep=",")[,2]
selection <- selection[selection != 0]
trueCSP <- read.table("../SimulatedGenome_CSP1.csv", sep=",", header=T)
dm.idx <- grepl(pattern = "^[A-Z].[ACGT]{3}.log", x = trueCSP[,1])
trueMutation <- trueCSP[dm.idx, 2]
trueSelection <- trueCSP[!dm.idx, 2]

plot(trueMutation, mutation)
points(trueSelection, selection, col="red")
abline(0, 1, col="blue", lwd=2)
cor(trueMutation, mutation)
cor(trueSelection, selection)

likTrace <- unlist(c(read.table("liklihoodTrace.csv", sep=",")))
plot(likTrace, type = "l")

which(log10(test.phi) > 0.5)
gene <- 303
expressionTrace <- read.table("expressionLevelTrace.csv", sep=",")
plot(log10(expressionTrace[, gene]), type = "l")
test.phiTrace0 <- read.table("phiTrace_nmix_0.csv", sep=",")
test.phiTrace1 <- read.table("phiTrace_nmix_1.csv", sep=",")
plot(log10(test.phiTrace1[, gene]), type = "l")
plot(log10(test.phiTrace0[, gene]), type = "l")


mixAssignment <- read.table("mixAssignment.csv", sep=",")[,2]
idx <- round(mixAssignment) == 1

test.scuo <- read.table("scuo.csv", sep=",")[,2]
plot(log10(test.scuo), log10(true.phi))

