test.phi0 <- read.table("expressionPosterior_Cat0.csv", sep=",")[,2]
test.phi1 <- read.table("expressionPosterior_Cat1.csv", sep=",")[,2]
test.phi.names <- as.character(read.table("expressionPosterior_Cat0.csv", sep=",")[,1])
true.phi <- read.table("../SimulatedGenome_allUnique_phi.csv", sep=",", header=T)[, 2]


idx <- c(rep(TRUE, times = 500), rep(FALSE, times = 500))
mixAssignment <- read.table("mixAssignment.csv", sep=",")[,2]
idx <- round(mixAssignment) == 1

plot(log10(true.phi[idx]), log10(test.phi1[idx]), xlim = range(log10(true.phi)), ylim = range(log10(test.phi1), na.rm = T))
points(log10(true.phi[!idx]), log10(test.phi0[!idx]), col = "red")
legend("topleft", legend = c("Category 0", "Category 1"), col = c("black", "red"), pch = c(1,1))
abline(0,1, col="blue", lwd=2)
cor(c(log10(true.phi[!idx]), log10(true.phi[idx])), c(log10(test.phi0[!idx]), log10(test.phi1[idx])))


mutation <- read.table("mutationPosterior_Cat0.csv", sep=",")[,2]
mutation <- mutation[mutation != 0]
selection <- read.table("selectionPosterior_Cat1.csv", sep=",")[,2]
selection <- selection[selection != 0]
trueCSP <- read.table("../SimulatedGenome_CSP0.csv", sep=",", header=T)
dm.idx <- grepl(pattern = "^[A-Z].[ACGT]{3}.log", x = trueCSP[,1])
trueMutation <- trueCSP[dm.idx, 2]
trueSelection <- trueCSP[!dm.idx, 2]

plot(trueMutation, mutation, xlab = "TRUE VALUES", ylab = "ESTM. VALUES", xlim = range(c(trueMutation, trueSelection)), 
     ylim = range(c(mutation, selection)) )
points(trueSelection, selection, col="red")
abline(0, 1, col="blue", lwd=2)
cor(trueMutation, mutation)
cor(trueSelection, selection)

likTrace <- unlist(c(read.table("liklihoodTrace.csv", sep=",")))
plot(likTrace[10:500], type = "l")

expressionTrace <- read.table("expressionLevelTrace.csv", sep=",")
which(log10(test.phi) > 2)
gene <- 756
plot(log10(expressionTrace[, gene]), type = "l")
test.phiTrace0 <- read.table("phiTrace_nmix_0.csv", sep=",")
test.phiTrace1 <- read.table("phiTrace_nmix_0.csv", sep=",")
ylim <- range(c(log10(test.phiTrace1[, gene]), log10(test.phiTrace0[, gene]), log10(true.phi[gene])))
ylim <- range(c(log10(test.phiTrace0[, gene])), log10(true.phi[gene]))
plot(log10(test.phiTrace0[, gene]), type = "l", ylim = ylim, ylab = "estm. phi", xlab="sample")
lines(log10(test.phiTrace1[, gene]), type = "l", col="red")
abline(h = log10(true.phi[gene]), lty = 2)
abline(h = ifelse(idx[gene], log10(test.phi1[gene]), log10(test.phi0[gene])), lty = 2, col = "blue")
lines(log10(expressionTrace[, gene]), type = "l", col="blue")

mutationcat0 <- read.table("mutationParamTrace_0.csv", sep=",", header=T)
mutationcat1 <- read.table("mutationParamTrace_1.csv", sep=",", header=T)
plot(mutationcat0[,1], type = "l")
lines(mutationcat1[,1], type = "l", col="red")

selectincat0 <- read.table("selectionParamTrace_0.csv", sep=",", header=T)
selectioncat1 <- read.table("selectionParamTrace_1.csv", sep=",", header=T)
plot(selectincat0[,10], type = "l")
lines(selectincat1[,1], type = "l", col="red")


test.scuo <- read.table("scuo.csv", sep=",")[,2]
plot(log10(true.phi), log10(test.scuo))
cor(log10(true.phi), log10(test.scuo))


