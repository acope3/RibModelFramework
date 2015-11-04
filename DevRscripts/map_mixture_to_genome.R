
mixture.probabilities <- read.table(file="1_kluyveri_full_mixProbs_50.csv", sep=",")
gene.names <- as.character(mixture.probabilities[,1])
prob.for.mix.one <- as.numeric(mixture.probabilities[,2])
prob.for.mix.two <- as.numeric(mixture.probabilities[,3])
mixture.probabilities <- data.frame(gene.name=gene.names, mix1=prob.for.mix.one, mix2=prob.for.mix.two)


genome.annotation <- read.table(file="../../../organisms/yeast/data/LKluyveri/Sakl_genome_annotation_minimal.csv", sep="\t", as.is=T)
genome.annotation <- genome.annotation[genome.annotation[,2] == "CDS", ]


kluyveri.chr.split.genome <- seqinr::read.fasta(file="../../../organisms/yeast/data/wholeGenomes/Lkluyveri_wholeGenome_chr.fasta")

plot.GC.By.Chr.with.Genes <- function(genome=NULL, annotation=NULL, mix.prob=NULL, gc.content=TRUE)
{
  nchr <- length(genome)
  num.panels <- nchr
  add.gene <- !is.null(annotation)
  if(add.gene)
  {
    num.panels <- num.panels + nchr
    add.certainty <- !is.null(mix.prob)
    if(add.certainty){num.panels <- num.panels + nchr}
  }
  par(mar=c(2, 4, 2, 2), oma=c(3,1,3,1))
  layout(matrix(1:num.panels, nrow=num.panels))  
  
  # collect gc content information for every chromosome over a window
  maxlength <- max(seqinr::getLength.list(genome))
  x.range <- c(-100, maxlength+100)
  for(curChr in genome)
  {
    n.nucleotides <- length(curChr)
    start <- 1
    window.size <- 5000
    window.move <- 1000
    gc.result <- vector("numeric")
    position.result <- vector("numeric")
    i <- 1
    while(T)
    {
      end <- start+window.size
      if(end > n.nucleotides) end <- n.nucleotides
      
      # get GC-content or GC-skew
      if(gc.content)
      {
        gc.result[i] <- sum(seqinr::count(seq = curChr[start:end], wordsize = 1, alphabet = c("g", "c"))) / window.size #gc content
      }else{
        g.count <- seqinr::count(seq = curChr[start:end], wordsize = 1, alphabet = c("g"))
        c.count <- seqinr::count(seq = curChr[start:end], wordsize = 1, alphabet = c("c"))
        gc.result[i] <- (g.count-c.count) / (g.count+c.count) #gc skew
      }
      position.result[i] <- (start + end) / 2
      start <- start + window.move
      if(end == n.nucleotides) break
      i <- i + 1    
    }
  
    
    chr.name <- seqinr::getName.SeqFastadna(curChr)
    plot(NULL, NULL, xlim=x.range, ylim=range(gc.result), xlab="Position", ylab="GC content", axes=F, main=chr.name)
    lines(position.result, gc.result)
    axis(side = 2, at = pretty(range(gc.result)), las = 1)
    abline(h=0.404, col="blue", lty=2)
    abline(h=0.529, col="blue", lty=2)
    
    if(add.gene)
    {
      idx <- grepl(pattern = chr.name, x = annotation[, 1])
      current.chr.subset <- annotation[idx,]
      current.mix.prob.subset <- mix.prob[idx,]
      plot(NULL, NULL, xlim=x.range, ylim=c(0.05, 0.2), xlab="Position", ylab="Strand", axes=F)
      for(j in 1:length(current.chr.subset[,1]))
      {
        if(as.character(current.chr.subset[j, 7]) == "1"){
          y.value <- c(0.1, 0.1)
        }else{
          y.value <- c(0.15, 0.15)
        }
        
        colour <- ribModel:::.mixtureColors[which(current.mix.prob.subset[j, -1] == max(current.mix.prob.subset[j, -1]))]
        lines(current.chr.subset[j, 3:4], y.value+runif(1, -0.01, 0.01), col=colour) 
      }
      axis(side = 2, at = c(0.1, 0.15), labels = c("+", "-"), las = 1)
      
      if(add.certainty)
      {
        plot(NULL, NULL, xlim=x.range, ylim=c(0.5, 1), xlab="Position", ylab="Certainty", axes=F)
        for(j in 1:length(current.chr.subset[,1]))
        {
          mix.index <- which(current.mix.prob.subset[j, -1] == max(current.mix.prob.subset[j, -1]))
          if(length(mix.index) > 1)
          {
            warning(paste("Gene", current.chr.subset[j,8], "on chromosome", current.chr.subset[,1], "can not be assign to a categorty\n",
                          "Gene will be ignored\n"))
            mix.index <- 1
          }
          colour <- ribModel:::.mixtureColors[mix.index]
          mix.certainty <- current.mix.prob.subset[j, mix.index + 1]
          lines(mean(as.numeric(current.chr.subset[j, 3:4])), mix.certainty, col=colour, type="h") 
        }
        axis(side = 2, at = c(0.5, 1), labels = c("0.5", "1.0"), las = 1)
      }
    }
    axis(side = 1, at = pretty(position.result))
  }
}

pdf(file="1_kluyveri_full_mixassignment_genomemap_50.pdf", height=26, width=10)
plot.GC.By.Chr.with.Genes(genome=kluyveri.chr.split.genome, annotation=genome.annotation, mix.prob=mixture.probabilities)
dev.off()