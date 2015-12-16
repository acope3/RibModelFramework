#' Initialize Parameter 
#' 
#' @param genome A genome object which the parameter use.
#' 
#' @param sphi Initial sphi values that corrosponds with the sphi for
#' each mixture. sphi is a vector whose length should equal numMixtures.
#' 
#' @param numMixtures The number of mixtures the parameter 
#' should use with the genome. This should be a positive number.
#' 
#' @param geneAssignment A vector holding that corrosponds to each
#' gene in the genome. The vector size should equal the number of 
#' genes in the genome. The assignment is to which mixture the 
#' gene initially starts in. Valid values for the vector go from
#' 1 to numMixtures.
#' 
#' @param expressionValues A vector containing the starting phi 
#' values. This is an optional argument and is set to null if
#' no value is given. If it is provided, the length of the 
#' vector should equal the number of genes in the genome.
#' 
#' @param model The string name of the model and parameter type 
#' to create. Valid options are "ROC", "RFP", or "FONSE". The default
#' value for model is "ROC".
#' 
#' @param split.serine Whether serine should be considered as 
#' one or two amino acids when running the model. TRUE and FALSE
#' are the only valid values. The default value for split.serine is
#' TRUE.
#' 
#' @param mixture.definition A string describing how each mixture should
#' be treated with respect to mutation and selection. Valid values consist
#' of "allUnique", "mutationShared", and "selectionShared". The default value
#' for mixture.definition is "allUnique". See details for more information.
#' 
#' @param mixture.definition.matrix A matrix representation of how
#' the mutation and selection categories corrospond to the mixtures.
#' The default value for mixture.definition.matrix is NULL. If provided,
#' the model will use the matrix to initialize the mutation and selection
#' categories instead of the definition listed directly above. See details
#' for more information.
#' 
#' @param  restart.file File name containing information to reinitialize a 
#' previous parameter object. If given, all other arguments will be ignored.
#' The default value for restart.file is NULL.
#' 
#' @param mutation_prior_sd Controling the standard deviation of the normal prior on the mutation parameters
#' 
#' @return parameter Returns an initialized parameter object.
#' 
#' @description \code{initializeParameterObject} will call the appropriate followup
#' call to writeXXXParameterObject based off of the value of model and restart.file.
#' 
#' @details \code{initializeParameterObject} checks the values of the arguments 
#' given to insure the values are valid. Additionally, if a restart file is given,
#' no follow up function calls are made - a new call is made instead which calls
#' the C++ constructor that only takes a file name.
#' 
#' The mixture definition and mixture definition matrix describe how the mutation
#' and selection categories are set up with respect to the number of mixtures. For
#' example, if mixture.definition = "allUnique" and numMixtures = 3, a matrix
#' representation would be as follows:
#' 
#' 1 1
#' 
#' 2 2
#' 
#' 3 3
#' 
#' where each row represents a mixture, the first column represents the mutation
#' category, and the second column represents the selection category. Another 
#' example would be mixture.definition = "selectionShared" and numMixtures = 4.
#' 
#' 1 1
#' 
#' 2 1
#' 
#' 3 1
#' 
#' 4 1
#' 
#' In this case, the selection category is the same for every mixture. If a matrix
#' is given, and it is valid, then the mutation/selection relationship will be
#' defined by the given matrix as opposed to the keyword. A matrix should only
#' be given in cases where the keywords would not create the desired valid matrix.
#' 
#' 
initializeParameterObject <- function(genome, sphi, numMixtures, geneAssignment, expressionValues = NULL, model = "ROC",
                                      split.serine = TRUE, mixture.definition = "allUnique", mixture.definition.matrix = NULL,
                                      restart.file = NULL, mutation_prior_sd = 0.35)
{
  # check input integrity
  if(length(sphi) != numMixtures){
    stop("Not all mixtures have an Sphi value assigned!\n")
  }
  
  if(length(genome) != length(geneAssignment)){
    stop("Not all Genes have a mixture assignment!\n")
  }
  
  if(max(geneAssignment) > numMixtures){
    stop("Gene is assigned to non existing mixture!\n")
  }  
  #TODO: should we check integraty of other values, such as numMixtures being
  #positive?
  
  
  
  if(model == "ROC"){
    if(is.null(restart.file)){
      parameter <- initializeROCParameterObject(genome, sphi, numMixtures, 
                            geneAssignment, expressionValues, split.serine, 
                            mixture.definition, mixture.definition.matrix, 
                            mutation_prior_sd)    
    }
    else{
      parameter <- new(ROCParameter, restart.file)
    }
  }
  else if(model == "FONSE"){
    if(is.null(restart.file)){
      parameter <- initializeFONSEParameterObject(genome, sphi, numMixtures, 
                            geneAssignment, expressionValues, split.serine, 
                            mixture.definition, mixture.definition.matrix)
    }
    else{
      parameter <- new(FONSEParameter, restart.file)
    }
  }
  else if(model == "RFP"){
    if(is.null(restart.file)){
      parameter <- initializeRFPParameterObject(genome, sphi, numMixtures, 
                            geneAssignment, expressionValues, split.serine, 
                            mixture.definition, mixture.definition.matrix) 
    }
    else{
      parameter <- new(RFPParameter, restart.file)
    }
  }
  else{
    stop("Unknown model.")
  }
  
  return(parameter)
}



initializeROCParameterObject <- function(genome, sphi, numMixtures, geneAssignment, expressionValues = NULL,
                                         split.serine = TRUE, mixture.definition = "allUnique", mixture.definition.matrix = NULL,
                                         mutation_prior_sd = 0.35)
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
    parameter$initializeSynthesisRateByGenome(genome, mean(sphi))
  }else{
    parameter$initializeSynthesisRateByList(expressionValues)
  }
  
  parameter$mutation_prior_sd <- (mutation_prior_sd)
  
  numMutationCategory <- parameter$numMutationCategories
  numSelectionCategory <- parameter$numSelectionCategories
  
  phi <- parameter$getCurrentSynthesisRateForMixture(1) # phi values are all the same initially
  names.aa <- aminoAcids()
  for(aa in names.aa)
  {
    if(aa == "M" || aa == "W" || aa == "X") next
    
    codonCounts <- getCodonCountsForAA(aa, genome)
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
      covmat[[mixElement]] <- splitMatrix(t(csp$R) %*% csp$R, numCodons, numCodons)  # we expect the covariance matrix, but get the decomposition.
    }
    compl.covMat <- matrix(0, ncol = numMixtures * numCodons * 2, nrow = numMixtures * numCodons * 2)
    matrix.positions <- subMatrices(compl.covMat, numCodons, numCodons)
    
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
    compl.covMat <- diag((numMutationCategory + numSelectionCategory) * numCodons) * 0.05
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
    
    codonCounts <- getCodonCountsForAA(aa, genome)
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
      covmat[[mixElement]] <- splitMatrix(t(csp$R) %*% csp$R, numCodons, numCodons)  # we expect the covariance matrix, but get the decomposition.
    }
    compl.covMat <- matrix(0, ncol = numMixtures * numCodons * 2, nrow = numMixtures * numCodons * 2)
    matrix.positions <- subMatrices(compl.covMat, numCodons, numCodons)
    
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

writeParameterToCSV <- function(parameter, filename, CSP, mixture, samples)
{
  UseMethod("writeParameterToCSV", parameter)
}

writeParameterToCSV.Rcpp_ROCParameter <- function(parameter, filename=NULL, CSP=NULL, mixture = 1, samples = 10)
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
        Value <- c(Value,parameter$getMutationPosteriorMeanForCodon(mixture, samples, codons[i]))
        Std_Deviation <- c(Std_Deviation, sqrt(parameter$getMutationVarianceForCodon(mixture, samples, codons[i], TRUE)))
      }
      else if(CSP == "Selection")
      {
        Value <- c(Value,parameter$getSelectionPosteriorMeanForCodon(mixture, samples, codons[i]))
        Std_Deviation <- c(Std_Deviation, sqrt(parameter$getSelectionVarianceForCodon(mixture, samples, codons[i], TRUE)))
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

getCodonCountsForAA <- function(aa, genome)
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

subMatrices <- function(M, r, c)
{
  rg <- (row(M) - 1) %/% r + 1
  cg <- (col(M) - 1) %/% c + 1
  rci <- (rg - 1) * max(cg) + cg
  return(rci)
}
splitMatrix <- function(M, r, c)
{
  rci <- subMatrices(M, r, c)
  N <- prod(dim(M)) / r / c
  cv <- lapply(1:N, function(x) M[rci==x])
  
  return(lapply(1:N, function(i) matrix(cv[[i]], nrow = r)))
} 


#' Write Parameter Object to a File
#' 
#' @param parameter A parameter object that corrosponds to
#' one of the model types, such as "ROC", or "FONSE".
#' 
#' @param file A filename that where the data will be stored.
#' The file should end with the extension "Rdat".
#' 
#' @return This function has no return value.
#' 
#' @description \code{writeParameterObject} will call the appropriate followup
#' call to writeParameterObject based off of the parameter type 
#' given.
#' 
#' @details For example, if a ROCParameter is passed, the the writeParameterObject
#' for the ROCParameter will be called. This allows us to not have an if-else
#' block in the code - making use of the how R handles these situations.
#' 
#' 
writeParameterObject <- function(parameter, file)
{
  UseMethod("writeParameterObject", parameter)
}


extractBaseInfo <- function(parameter)
{
  trace <- parameter$getTraceObject()
  sPhiTraces <- trace$getSphiTraces()
  sphiAcceptRatTrace <- trace$getSphiAcceptanceRatioTrace()
  synthRateTrace <- trace$getSynthesisRateTrace()
  synthAcceptRatTrace <- trace$getSynthesisRateAcceptanceRatioTrace()
  mixAssignTrace <- trace$getMixutreAssignmentTrace()
  mixProbTrace <- trace$getMixtureProbabilitiesTrace()
  cspAcceptRatTrace <- trace$getCspAcceptanceRatioTrace()
  numMix <- parameter$numMixtures
  numMut <- parameter$numMutationCategories
  numSel <- parameter$numSelectionCategories
  categories <- parameter$getCategories()
  curMixAssignment <- parameter$getMixtureAssignment()
  lastIteration <- parameter$getLastIteration()
  
  varList <- list(sPhiTraces = sPhiTraces, 
                    sphiAcceptRatTrace = sphiAcceptRatTrace,
                    synthRateTrace = synthRateTrace,
                    synthAcceptRatTrace = synthAcceptRatTrace,
                    mixAssignTrace = mixAssignTrace,
                    mixProbTrace = mixProbTrace,
                    cspAcceptRatTrace = cspAcceptRatTrace,
                    numMix = numMix,
                    numMut = numMut,
                    numSel = numSel,
                    categories = categories
                    )
  return(varList)
}


writeParameterObject.Rcpp_ROCParameter <- function(parameter, file)
{
  paramBase <- extractBaseInfo(parameter)
  
  currentMutation <- parameter$currentMutationParameter
  currentSelection <- parameter$currentSelectionParameter
  proposedMutation <- parameter$proposedMutationParameter
  proposedSelection <- parameter$proposedSelectionParameter
  
  trace <- parameter$getTraceObject()
  
  mutationTrace <- trace$getMutationParameterTrace()
  selectionTrace <- trace$getSelectionParameterTrace()
  aphiAcceptRatTrace <- trace$getAphiAcceptanceRatioTrace()
  aphiTrace <- trace$getAphiTraces()
  sepisolonTrace <- trace$getSepsilonTraces()
  
  save(list = c("paramBase", "currentMutation", "currentSelection",
                "proposedMutation", "proposedSelection",  
                "mutationTrace", "selectionTrace", 
                "aphiAcceptRatTrace", "aphiTrace", "sepisolonTrace"),
       file=file)
}


writeParameterObject.Rcpp_RFPParameter <- function(parameter, file)
{
  paramBase <- extractBaseInfo(parameter)
  
  currentAlpha <- parameter$currentAlphaParameter
  currentLambdaPrime <- parameter$currentLambdaPrimeParameter
  proposedAlpha <- parameter$proposedAlphaParameter
  proposedLambdaPrime <- parameter$proposedLambdaPrimeParameter
  
  
  trace <- parameter$getTraceObject()
  alphaTrace <- trace$getAlphaParameterTrace()
  lambdaPrimeTrace <- trace$getLambdaPrimeParameterTrace()

  save(list = c("paramBase", "currentAlpha", "currentLambdaPrime", "proposedAlpha",
                "proposedLambdaPrime", "alphaTrace", "lambdaPrimeTrace"),
       file=file)
}


writeParameterObject.Rcpp_FONSEParameter <- function(parameter, file)
{
  trace <- parameter$getTraceObject()
  sPhiTraces <- trace$getSphiTraces()
  sphiAcceptRatTrace <- trace$getSphiAcceptanceRatioTrace()
  synthRateTrace <- trace$getSynthesisRateTrace()
  synthAcceptRatTrace <- trace$getSynthesisRateAcceptanceRatioTrace()
  mixAssignTrace <- trace$getMixutreAssignmentTrace()
  mixProbTrace <- trace$getMixtureProbabilitiesTrace()
  cspAcceptRatTrace <- trace$getCspAcceptanceRatioTrace()
  mutationTrace <- trace$getMutationParameterTrace()
  selectionTrace <- trace$getSelectionParameterTrace()
  save(list = c("sPhiTraces", "sphiAcceptRatTrace", "synthRateTrace", "synthAcceptRatTrace", 
                "mixAssignTrace", "mixProbTrace", "cspAcceptRatTrace", "mutationTrace", "selectionTrace",
                "loglikeTrace"),
       file=file)
}


loadParameterObject <- function(parameter, file, model)
{
  UseMethod("loadParameterObject", parameter)
}


setBaseInfo <- function(parameter, file, model)
{
  tempEnv <- new.env();
  load(file = file, envir = tempEnv)
  parameter$setCategories(tempEnv$paramBase$categories)
  parameter$setCategoriesForTrace()
  parameter$numMixtures <- tempEnv$paramBase$numMix
  parameter$numMutationCategories <- tempEnv$paramBase$numMut
  parameter$numSelectionCategories <- tempEnv$paramBase$numSel
  parameter$setMixtureAssignment(tempEnv$paramBase$curMixAssignment)
  parameter$setLastIteration(tempEnv$paramBase$lastIteration)
  
  trace <- parameter$getTraceObject()
  trace$setSphiTraces(tempEnv$paramBase$sPhiTraces)
  trace$setSphiAcceptanceRatioTrace(tempEnv$paramBase$sphiAcceptRatTrace)
  trace$setSynthesisRateTrace(tempEnv$paramBase$synthRateTrace)
  trace$setSynthesisRateAcceptanceRatioTrace(tempEnv$paramBase$synthAcceptRatTrace)
  trace$setMixtureAssignmentTrace(tempEnv$paramBase$mixAssignTrace)
  trace$setMixtureProbabilitiesTrace(tempEnv$paramBase$mixProbTrace)
  trace$setCspAcceptanceRatioTrace(tempEnv$paramBase$cspAcceptRatTrace)
  if (model == "ROC"){
    parameter$setROCTrace(trace)
  }
  else if (model == "RFP"){
    parameter$setRFPTrace(trace)
  }
  else if (model == "FONSE"){
    parameter$setFONSETrace(trace)
  }
  return(parameter)
}


loadParameterObject.Rcpp_ROCParameter <- function(parameter, file, model)
{
  tempEnv <- new.env();
  load(file = file, envir = tempEnv)
  setBaseInfo(parameter, file, model)
  parameter$currentMutationParameter <- tempEnv$currentMutation
  parameter$currentSelectionParameter <- tempEnv$currentSelection
  parameter$proposedMutationParameter <- tempEnv$proposedMutation
  parameter$proposedSelectionParameter <- tempEnv$proposedSelection
  trace <- parameter$getTraceObject()
  trace$setAphiTrace(tempEnv$aphiTrace)
  trace$setAphiAcceptanceRatioTrace(tempEnv$aphiAcceptRatTrace)
  trace$setSepsilonTrace(tempEnv$sepisolonTrace)
  trace$setMutationParameterTrace(tempEnv$mutationTrace)
  trace$setSelectionParameterTrace(tempEnv$selectionTrace)
  
  
  parameter$setROCTrace(trace)
  return(parameter) #Because of concerns with R passing arguments, we return explicitly.
}


loadParameterObject.Rcpp_RFPParameter <- function(parameter, file, model)
{
  tempEnv <- new.env();
  load(file = file, envir = tempEnv)
  setBaseInfo(parameter, file, model)
  parameter$currentAlphaParameter <- tempEnv$currentAlpha
  parameter$proposedAlphaParameter <- tempEnv$proposedAlpha
  parameter$currentLambdaPrimeParameter <- tempEnv$currentLambdaPrime
  parameter$proposedLambdaPrimeParameter <- tempEnv$proposedLambdaPrime
  
  trace <- parameter$getTraceObject()
  trace$setAlphaParameterTrace(tempEnv$alphaTrace)
  trace$setLambdaPrimeParameterTrace(tempEnv$lambdaPrimeTrace)
  
  parameter$setRFPTrace(trace)
  return(parameter) #R seems to produce copies, not pointers.
}