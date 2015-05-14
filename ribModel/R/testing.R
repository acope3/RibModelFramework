test.phi <- read.table("/home/clandere/CodonUsageBias/organisms/yeast/results/test.phi", sep=",")[,2]
test.phi.names <- as.character(read.table("/home/clandere/CodonUsageBias/organisms/yeast/results/test.phi", sep=",")[,1])

kluyveri.phi <- read.table("/home/clandere/CodonUsageBias/organisms/yeast/results/Skluyvery_study/Skluyveri_chromosomeA/without_xobs_singlechain_main.phi", sep=",", header=T)[, 2]
test.scuo <- read.table("/home/clandere/CodonUsageBias/organisms/yeast/results/test.scuo", sep=",")[,2]

test.phiTrace <- read.table("/home/clandere/CodonUsageBias/organisms/yeast/results/test.phiTrace", sep=",")



plot(test.phiTrace[, 5], type = "l")
#bla <- colMeans(test.phiTrace[500:1000, 1:448])

test.lik <- as.numeric(scan("/home/clandere/CodonUsageBias/organisms/yeast/results/test.lik"))
test.sphi <- as.numeric(scan("/home/clandere/CodonUsageBias/organisms/yeast/results/test.sphi"))
#plot(test.lik[10:90], type="l")
plot(test.lik, type="l")
plot(50:300, test.lik[50:300], type="l")

plot(test.sphi, type="l")

#plot(log(test.phi[test.phi > 0.0067]), log(kluyveri.phi[test.phi > 0.0067]))
#plot(log(bla[bla > 0.0067]), log(kluyveri.phi[bla > 0.0067]))
plot(log(test.scuo), log(kluyveri.phi))
plot(log(test.phi), log(kluyveri.phi))

cor(log(test.phi), log(kluyveri.phi))

mean(test.phi)
mean(kluyveri.phi)
mean(test.lik[150:300])

library(cubfits)
seq.string <- readGenome("/home/clandere/CodonUsageBias/organisms/yeast/data/LKluyveri/Skluyveri_chromosomeA.fasta")
aa.names <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "Z")
codon.counts <- gen.scuo(seq.string, aa.names)
scuo <- calc_scuo_values(codon.counts)$SCUO

aa.names <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "Z")
y <- gen.y(seq.string, aa.names)




