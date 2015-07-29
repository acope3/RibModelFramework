library(ribModel)
rm(list=ls())


#read genome
genome <- initializeGenomeObject(file = "../ribModel/data/rfp.counts.by.codon.and.gene.GSE63789.wt.csv", FALSE)


#initialize parameter object
sphi_init <- 2
numMixtures <- 1
mixDef <- "allUnique"
geneAssignment <- c(rep(1, genome$getGenomeSize()))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, model= "RFP", split.serine = TRUE, mixture.definition = mixDef)


#init from "true" values

#DO I NEED THE PREVIOUS PHI VALUES???????
parameter$initMutationSelectionCategories("FILENAME HERE", 1, "Alpha")
parameter$initMutationSelectionCategories("FILENAME HERE", 1, "LambdaPrime")
