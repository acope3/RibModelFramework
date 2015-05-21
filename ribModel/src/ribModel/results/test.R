test.phi <- read.table("phiPosterior.csv", sep=",")[,2]
test.phi.names <- as.character(read.table("phiPosterior.csv", sep=",")[,1])
kluyveri.phi <- read.table("../Skluyveri_ChrA_phi_est.csv", sep=",", header=T)[, 2]
plot(log10(kluyveri.phi), log10(test.phi))
cor(log10(kluyveri.phi), log10(test.phi))

test.phi.names[log(test.phi) > 4]
which(log(test.phi) > 4)

if(has.emp.phi)
{
  seq.string.names <- test.phi
  names(seq.string.names) <- test.phi.names
  emp <- read.empirical.data("~/CodonUsageBias/organisms/yeast/data/LKluyveri/expression/Skluyveri_GSM552569.csv", seq.string.names, 1, th=0)
  emp <- emp$empirical[names(emp$empirical) %in% test.phi.names]
  emp <- emp[order(names(emp))]
}
idx <- log10(test.phi) > -3
idx <- grepl(pattern = "^SAKL0C", x = test.phi.names)
plot(log10(emp[!idx]), log10(test.phi[!idx]))
points(log10(emp[idx]), log10(test.phi[idx]), col="red")
cor(log10(emp[!idx]), log10(test.phi[!idx]))^2
cor(log10(emp[idx]), log10(test.phi[idx]))^2

test.phiTrace <- read.table("phiTrace_nmix_0.csv", sep=",")
which(log(test.phi) < -50)
plot(log(test.phiTrace[, 173]), type = "l")


mixAssignment <- read.table("mixAssignment.csv", sep=",")[,2]
mixAssignment[272]



test.scuo <- read.table("scuo.csv", sep=",")[,2]
plot(log10(test.scuo), log10(kluyveri.phi))

csp.kluyveri <- read.table("../Skluyveri_CSP_ChrA.csv", sep=",")



