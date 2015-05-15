test.phi <- read.table("test.phi", sep=",")[,2]
test.phi.names <- as.character(read.table("test.phi", sep=",")[,1])
kluyveri.phi <- read.table("../Skluyveri_ChrA_phi_est.csv", sep=",", header=T)[, 2]
plot(log(test.phi), log(kluyveri.phi))