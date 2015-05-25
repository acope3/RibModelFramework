test.phi <- read.table("phiPosterior.csv", sep=",")[,2]
test.phi.names <- as.character(read.table("phiPosterior.csv", sep=",")[,1])
kluyveri.phi <- read.table("../Skluyveri_ChrA_ChrCleft_phi_est.csv", sep=",", header=T)[, 2]

idx <- grepl(pattern = "^SAKL0C", x = test.phi.names)
plot(log10(kluyveri.phi[!idx]), log10(test.phi[!idx]), xlim = range(log10(kluyveri.phi)), ylim = range(log10(test.phi)))
points(log10(kluyveri.phi[idx]), log10(test.phi[idx]), col = "red")
abline(0,1, col="blue", lwd=2)
cor(log10(kluyveri.phi), log10(test.phi))



test.phi.names[log(test.phi) > 4]
which(log(test.phi) > 4)

if(has.emp.phi)
{
  seq.string.names <- test.phi
  names(seq.string.names) <- do.call("rbind", strsplit(test.phi.names, split = "_"))[, 1]
  emp <- read.empirical.data("~/CodonUsageBias/organisms/yeast/data/LKluyveri/expression/Skluyveri_GSM552569.csv", seq.string.names, 1, th=0)
  emp <- emp$empirical[names(emp$empirical) %in% test.phi.names]
  emp <- emp[order(names(emp))]
}
idx <- log10(test.phi) > -3
idx <- grepl(pattern = "^SAKL0C", x = test.phi.names)
plot(log10(emp[!idx]), log10(test.phi[!idx]), xlim = range(log10(emp)), ylim = range(log10(test.phi)))
points(log10(emp[idx]), log10(test.phi[idx]), col="red")
cor(log10(emp[!idx]), log10(test.phi[!idx]))^2
cor(log10(emp[idx]), log10(test.phi[idx]))^2

test.phiTrace0 <- read.table("phiTrace_nmix_0.csv", sep=",")
test.phiTrace1 <- read.table("phiTrace_nmix_1.csv", sep=",")
mixAssignment <- read.table("mixAssignment.csv", sep=",")[,2]
idx <- round(mixAssignment) == 1
which(log10(test.phi) > 3)
gene <- 208
if(idx[gene]){
  plot(log10(test.phiTrace1[, gene]), type = "l")
}else{
  plot(log10(test.phiTrace0[, gene]), type = "l")
}

expressionTrace <- read.table("expressionLevelTrace.csv", sep=",")
plot(log10(expressionTrace[, 208]), type = "l")


plot(log10(kluyveri.phi[!idx]), log10(test.phi[!idx]), xlim = range(log10(kluyveri.phi)), ylim = range(log10(test.phi)))
points(log10(kluyveri.phi[idx]), log10(test.phi[idx]), col = "red")

hist(log10(test.phi), nclass = 100)

test.scuo <- read.table("scuo.csv", sep=",")[,2]
plot(log10(test.scuo), log10(kluyveri.phi))

csp.kluyveri <- read.table("../Skluyveri_CSP_ChrA.csv", sep=",")

likTrace <- unlist(c(read.table("liklihoodTrace.csv", sep=",")))
plot(likTrace, type = "l")


