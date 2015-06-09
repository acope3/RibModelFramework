
initializeParameterObject <- function(genome, sphi, numMixtures, geneAssignment, expressionValues = NULL,
                                      split.serine = TRUE, mixture.definition = "allUnique")
{
  # create parameter object
  parameter <- new(ROCParameter, genome$getGenomeSize(), sphi, numMixtures, geneAssignment, split.serine, mixture.definition)
  # initialize expression values
  if(is.null(expressionValues)){
    parameter$initializeExpressionByGenome(genome, sphi)
  }else{
    parameter$initializeExpressionByList(expressionValues)
  }
  
  
  phi <- parameter$getCurrentExpressionForMixture(1) # phi values are all the same initially
  names.aa <- aminoAcids()
  for(aa in names.aa)
  {
    if(aa == "M" || aa == "W" || aa == "X") next
    
    # TODO WORKS CURRENTLY ONLY FOR ALLUNIQUE!
    covmat <- NULL
    for(mixElement in 1:numMixtures)
    {    
      idx <- geneAssignment == mixElement
      codonCounts <- getCodonCountsForAA(aa)
      csp <- getCSPbyLogit(codonCounts[idx, ], phi[idx])
      
      parameter$initMutation(csp$coef.mat[1,], mixElement, aa)
      parameter$initSelection(csp$coef.mat[2,], mixElement, aa)
      covmat <- t(csp$R) %*% csp$R  # we expect the covariance matrix, but get the decomposition.
    }
    parameter$initCovarianceMatrix(covmat, aa)
  }
  
  return(parameter)
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
