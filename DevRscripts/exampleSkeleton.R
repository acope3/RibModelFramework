
gene.files <- system("ls -1 rokas/*.fasta", intern=TRUE)
queue = "medium*"
for (gene.index in sequence(length(gene.files))) {
    cat("source('selacNewton.R')\nyeast.gene<-read.dna('",gene.files[gene.index],"', format='fasta')\nyeast.gene<-as.list(as.matrix(cbind(yeast.gene))[1:7,])\nchars<-DNAbinToCodonNumeric(yeast.gene)\ntree<-read.tree('rokas/Rokas_bestTree.rooted.2.tre')\ntree.pruned <- drop.tip(tree,'Calb')\nresult<-EstimateParametersCodon(codon.data=chars, phy=tree.pruned)\nsave(result,file='gene",gene.index,".Rdata')",file=paste("gene",gene.index,".R",sep=""),sep="")
    pbsCommands=paste('#!/bin/bash','#$ -cwd',sep="\n")
    pbsCommands=paste(pbsCommands,'\n#$ -q ',queue, sep="")
    pbsCommands=paste(pbsCommands,'#$ -M lefty913@gmail.com', '#$ -m beas', '#$ -S /bin/bash',sep="\n")
    pbsCommands=paste(pbsCommands,"\n","/data/apps/R/3.1.0/bin/R CMD BATCH gene", gene.index, ".R",sep="")
    pbsFile <- paste("gene", gene.index, ".sh", sep="")
    cat(pbsCommands,file=pbsFile,append=FALSE)
    system(paste("chmod u+x ", pbsFile, sep=""))
    system(paste("/opt/grid/bin/lx26-amd64/qsub ", pbsFile, sep=""))
}

