
initializeParameterObject <- function(genome, sphi, numMixtures, geneAssignment, expressionValues = NULL, model = "ROC",
                                      split.serine = TRUE, mixture.definition = "allUnique", mixture.definition.matrix = NULL,
                                      restart.file = NULL)
{
  if(model == "ROC")
  {
    if(is.null(restart.file))
    {
      parameter <- initializeROCParameterObject(genome, sphi, numMixtures, geneAssignment, expressionValues,
                                   split.serine, mixture.definition, mixture.definition.matrix)    
    }else{
      parameter <- new(ROCParameter, restart.file)
    }
  }else if(model == "FONSE"){
    if(is.null(restart.file))
    {
      parameter <- initializeFONSEParameterObject(genome, sphi, numMixtures, geneAssignment, expressionValues,
                                    split.serine, mixture.definition, mixture.definition.matrix)
    }
  }else if(model == "RFP"){
    if(is.null(restart.file))
    {
      parameter <- initializeRFPParameterObject(genome, sphi, numMixtures, geneAssignment, expressionValues,
                                                split.serine, mixture.definition, mixture.definition.matrix) 
    }
    else{
      parameter <- new(RFPParameter, restart.file)
    }
  }else{
    stop("Unknown model.")
  }
  return(parameter)
}

initializeROCParameterObject <- function(genome, sphi, numMixtures, geneAssignment, expressionValues = NULL,
                                          split.serine = TRUE, mixture.definition = "allUnique", mixture.definition.matrix = NULL)
{
  # test input integrity
  if(genome$getGenomeSize() != length(geneAssignment)) 
  {
    stop("Gene assignment length does not match genome size. Not every Gene has a mixture assignment")
  }
  # create parameter object
  if(is.null(mixture.definition.matrix))
  { # keyword constructor
    parameter <- new(ROCParameter, sphi, numMixtures, geneAssignment, split.serine, mixture.definition)
  }else{
    #matrix constructor
    mixture.definition <- c(mixture.definition.matrix[, 1], mixture.definition.matrix[, 2])
    parameter <- new(ROCParameter, sphi, numMixtures, geneAssignment, mixture.definition, split.serine)
  }
  # initialize expression values
  if(is.null(expressionValues)){
    parameter$initializeSynthesisRateByGenome(genome, sphi)
  }else{
    parameter$initializeSynthesisRateByList(expressionValues)
  }
  
  numMutationCategory <- parameter$numMutationCategories
  numSelectionCategory <- parameter$numSelectionCategories
  
  phi <- parameter$getCurrentSynthesisRateForMixture(1) # phi values are all the same initially
  names.aa <- aminoAcids()
  for(aa in names.aa)
  {
    if(aa == "M" || aa == "W" || aa == "X") next
    
    codonCounts <- getCodonCountsForAA(aa)
    numCodons <- dim(codonCounts)[2] - 1
    #-----------------------------------------
    # TODO WORKS CURRENTLY ONLY FOR ALLUNIQUE!
    #-----------------------------------------
    covmat <- vector("list", numMixtures)
    for(mixElement in 1:numMixtures)
    {    
      idx <- geneAssignment == mixElement
      csp <- getCSPbyLogit(codonCounts[idx, ], phi[idx])
      
      parameter$initMutation(csp$coef.mat[1,], mixElement, aa)
      parameter$initSelection(csp$coef.mat[2,], mixElement, aa)
      # split matrix into sup matrices (dM and dEta)
      covmat[[mixElement]] <- split.matrix(t(csp$R) %*% csp$R, numCodons, numCodons)  # we expect the covariance matrix, but get the decomposition.
    }
    compl.covMat <- matrix(0, ncol = numMixtures * numCodons * 2, nrow = numMixtures * numCodons * 2)
    matrix.positions <- sub.matrices(compl.covMat, numCodons, numCodons)
    
    compl.seq <- seq(1, dim(compl.covMat)[1], numCodons)
    mut.seq <- compl.seq[1:(length(compl.seq)/2)]
    i <- 1
    for(pos in mut.seq)
    { 
      compl.covMat[matrix.positions == matrix.positions[pos, pos]] <- unlist(covmat[[i]][1])
      i <- i + 1
      i <- ifelse(i > numMutationCategory, 1, i)
    }
    sel.seq <- compl.seq[(length(compl.seq)/2 + 1):length(compl.seq)]
    i <- 1
    for(pos in sel.seq)
    { 
      compl.covMat[matrix.positions == matrix.positions[pos, pos]] <- unlist(covmat[[i]][4])
      i <- i + 1
      i <- ifelse(i > numMutationCategory, 1, i)
    }
    
    ofdiag.seq <- mut.seq + numCodons*numMutationCategory
    for(i in 1:length(mut.seq))
    {
      compl.covMat[matrix.positions == matrix.positions[mut.seq[i], ofdiag.seq[i]]] <- unlist(covmat[[i]][2])
      compl.covMat[matrix.positions == matrix.positions[ofdiag.seq[i], mut.seq[i]]] <- unlist(covmat[[i]][3])
    }
    #for testing
    compl.covMat <- diag((numMutationCategory + numSelectionCategory) * numCodons) *0.05
    #compl.covMat / max(compl.covMat)
    parameter$initCovarianceMatrix(compl.covMat, aa)
    
  }
  
  return(parameter)
}


initializeRFPParameterObject <- function(genome, sphi, numMixtures, geneAssignment, expressionValues = NULL,
                                         split.serine = TRUE, mixture.definition = "allUnique", mixture.definition.matrix = NULL)
{
  # test input integrity
  if(genome$getGenomeSize() != length(geneAssignment)) 
  {
    stop("Gene assignment length does not match genome size. Not every Gene has a mixture assignment")
  }
  # create parameter object
  if(is.null(mixture.definition.matrix))
  { # keyword constructor
    parameter <- new(RFPParameter, sphi, numMixtures, geneAssignment, split.serine, mixture.definition)
  }else{
    #matrix constructor
    mixture.definition <- c(mixture.definition.matrix[, 1], mixture.definition.matrix[, 2])
    parameter <- new(RFPParameter, sphi, numMixtures, geneAssignment, mixture.definition, split.serine)
  }
  # initialize expression values
  if(is.null(expressionValues)){
    parameter$initializeSynthesisRateByGenome(genome, sphi)
  }else{
    parameter$initializeSynthesisRateByList(expressionValues)
  }
  
  return (parameter)
}



initializeFONSEParameterObject <- function(genome, sphi, numMixtures, geneAssignment, expressionValues = NULL,
                                         split.serine = TRUE, mixture.definition = "allUnique", mixture.definition.matrix = NULL)
{
  # test input integrity
  if(genome$getGenomeSize() != length(geneAssignment)) 
  {
    stop("Gene assignment length does not match genome size. Not every Gene has a mixture assignment")
  }
  # create parameter object
  if(is.null(mixture.definition.matrix))
  { # keyword constructor
    parameter <- new(FONSEParameter, sphi, numMixtures, geneAssignment, split.serine, mixture.definition)
  }else{
    #matrix constructor
    mixture.definition <- c(mixture.definition.matrix[, 1], mixture.definition.matrix[, 2])
    parameter <- new(FONSEParameter, sphi, numMixtures, geneAssignment, mixture.definition, split.serine)
  }
  # initialize expression values
  if(is.null(expressionValues)){
    parameter$initializeSynthesisRateByGenome(genome, sphi)
  }else{
    parameter$initializeSynthesisRateByList(expressionValues)
  }
  
  numMutationCategory <- parameter$numMutationCategories
  numSelectionCategory <- parameter$numSelectionCategories
  
  phi <- parameter$getCurrentSynthesisRateForMixture(1) # phi values are all the same initially
  names.aa <- aminoAcids()
  for(aa in names.aa)
  {
    if(aa == "M" || aa == "W" || aa == "X") next
    
    codonCounts <- getCodonCountsForAA(aa)
    numCodons <- dim(codonCounts)[2] - 1
    #-----------------------------------------
    # TODO WORKS CURRENTLY ONLY FOR ALLUNIQUE!
    #-----------------------------------------
    covmat <- vector("list", numMixtures)
    for(mixElement in 1:numMixtures)
    {    
      idx <- geneAssignment == mixElement
      csp <- getCSPbyLogit(codonCounts[idx, ], phi[idx])
      
      parameter$initMutation(csp$coef.mat[1,], mixElement, aa)
      parameter$initSelection(csp$coef.mat[2,], mixElement, aa)
      # split matrix into sup matrices (dM and dEta)
      covmat[[mixElement]] <- split.matrix(t(csp$R) %*% csp$R, numCodons, numCodons)  # we expect the covariance matrix, but get the decomposition.
    }
    compl.covMat <- matrix(0, ncol = numMixtures * numCodons * 2, nrow = numMixtures * numCodons * 2)
    matrix.positions <- sub.matrices(compl.covMat, numCodons, numCodons)
    
    compl.seq <- seq(1, dim(compl.covMat)[1], numCodons)
    mut.seq <- compl.seq[1:(length(compl.seq)/2)]
    i <- 1
    for(pos in mut.seq)
    { 
      compl.covMat[matrix.positions == matrix.positions[pos, pos]] <- unlist(covmat[[i]][1])
      i <- i + 1
      i <- ifelse(i > numMutationCategory, 1, i)
    }
    sel.seq <- compl.seq[(length(compl.seq)/2 + 1):length(compl.seq)]
    i <- 1
    for(pos in sel.seq)
    { 
      compl.covMat[matrix.positions == matrix.positions[pos, pos]] <- unlist(covmat[[i]][4])
      i <- i + 1
      i <- ifelse(i > numMutationCategory, 1, i)
    }
    
    ofdiag.seq <- mut.seq + numCodons*numMutationCategory
    for(i in 1:length(mut.seq))
    {
      compl.covMat[matrix.positions == matrix.positions[mut.seq[i], ofdiag.seq[i]]] <- unlist(covmat[[i]][2])
      compl.covMat[matrix.positions == matrix.positions[ofdiag.seq[i], mut.seq[i]]] <- unlist(covmat[[i]][3])
    }
    #for testing
    compl.covMat <- diag((numMutationCategory + numSelectionCategory) * numCodons) *0.05
    #compl.covMat / max(compl.covMat)
    parameter$initCovarianceMatrix(compl.covMat, aa)
    
  }
  
  return(parameter)
}

writeParameterToCSV <- function(parameter, filename, CSP, mixture)
{
  UseMethod("writeParameterToCSV", parameter)
}

writeParameterToCSV.Rcpp_ROCParameter <- function(parameter, filename=NULL, CSP=NULL, mixture=1)
{
  names.aa <- aminoAcids()
  Amino_Acid <- c()
  Value <- c()
  Codon <- c()
  Std_Deviation <- c()
  for(aa in names.aa)
  {
    if(aa == "M" || aa == "W" || aa == "X") next
    codons <- AAToCodon(aa, T)
    for(i in 1:length(codons))
    {
      Amino_Acid <- c(Amino_Acid, aa)
      Codon <- c(Codon, codons[i])
      if(CSP == "Mutation")
      {
        Value <- c(Value,parameter$getMutationPosteriorMeanForCodon(mixture, samples*0.1, codons[i]))
        Std_Deviation <- c(Std_Deviation, sqrt(parameter$getMutationVarianceForCodon(mixture, samples*0.1, codons[i], TRUE)))
      }
      else if(CSP == "Selection")
      {
        Value <- c(Value,parameter$getSelectionPosteriorMeanForCodon(mixture, samples*0.1, codons[i]))
        Std_Deviation <- c(Std_Deviation, sqrt(parameter$getSelectionVarianceForCodon(mixture, samples*0.1, codons[i], TRUE)))
      }else 
      {
        stop("Unknown Parameter type given")
      }
    }
  }
  data <- data.frame(Amino_Acid,Codon,Value, Std_Deviation)
  if(is.null(filename))
  {
    print(data)
  }else 
  {
    write.csv(data, file = filename, row.names = FALSE, quote=FALSE)
  }
}

getCodonCountsForAA <- function(aa)
{
  # get codon count for aa
  codons <- AAToCodon(aa, F)
  codonCounts <- lapply(codons, function(codon){
    codonCounts <- genome$getCodonCountsPerGene(codon)
  })
  codonCounts <- do.call("cbind", codonCounts)
  return(codonCounts)
}

getCSPbyLogit <- function(codonCounts, phi, coefstart = NULL, x.arg = FALSE, y.arg = FALSE, qr.arg = FALSE)
{
  #avoid cases with 0 aa count
  idx <- rowSums(codonCounts) != 0
  
  ### Obtain new beta (M, S_1) from vglm.
  ret <- VGAM::vglm(codonCounts[idx, ] ~ phi[idx],
                    VGAM::multinomial, coefstart = coefstart,
                    x.arg = x.arg, y.arg = y.arg, qr.arg = qr.arg)
  coefficients <- ret@coefficients
  ## convert delta.t to delta.eta
  coefficients <- -coefficients
  
  ret <- list(coefficients = coefficients,
              coef.mat = matrix(coefficients, nrow = 2, byrow = TRUE),
              R = ret@R)
}

sub.matrices <- function(M, r, c)
{
  rg <- (row(M) - 1) %/% r + 1
  cg <- (col(M) - 1) %/% c + 1
  rci <- (rg - 1) * max(cg) + cg
  return(rci)
}
split.matrix <- function(M, r, c)
{
  rci <- sub.matrices(M, r, c)
  N <- prod(dim(M)) / r / c
  cv <- lapply(1:N, function(x) M[rci==x])
  
  return(lapply(1:N, function(i) matrix(cv[[i]], nrow = r)))
} 