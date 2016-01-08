#' Initialize Parameter 
#' 
#' @param genome An object of type Genome necessary for the initialization of the Parameter object. The default value 
#' is NULL.
#' 
#' @param sphi Initial values for sphi. Expected is a vector of length numMixtures.
#' The default value is NULL.
#' 
#' @param numMixtures The number of mixtures elements for the underlying mixture distribution (numMixtures > 0).
#' The default value is 1 
#' 
#' @param geneAssignment A vector holding the initial mixture assignment for each gene. 
#' The vector length has to equal the number of 
#' genes in the genome. Valid values for the vector range from
#' 1 to numMixtures. It is possible but not advised to leave a mixture element empty. The default Value is NULL.
#' 
#' @param expressionValues (Optional) A vector with intial phi values.
#' The length of the vector has to equal the number of genes in the Genome object.
#' The default value is NULL.
#' 
#' @param model Specifies the model used. Valid options are "ROC", "RFP", or "FONSE". The default
#' model is "ROC". 
#' ROC is described in Gilchrist et al. 2015
#' RFP and FONSE are currently unpublished 
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
#' @param mutation_prior_sd Controling the standard deviation of the normal 
#' prior on the mutation parameters
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
#' If expressionValues is given, then the phi values will be initialized with 
#' them. If not, we calculate starting phi values by doing an SCUO caculation.
#' 
initializeParameterObject <- function(genome = NULL, sphi = NULL, numMixtures = 1, 
                                    geneAssignment = NULL, expressionValues = NULL,
                                    model = "ROC", split.serine = TRUE, 
                                    mixture.definition = "allUnique", 
                                    mixture.definition.matrix = NULL,
                                    restart.file = NULL, mutation_prior_sd = 0.35){
  # check input integrity
  if(is.null(restart.file)){
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
  }

  
  
  
  if(model == "ROC"){
    if(is.null(restart.file)){
      parameter <- initializeROCParameterObject(genome, sphi, numMixtures, 
                            geneAssignment, expressionValues, split.serine, 
                            mixture.definition, mixture.definition.matrix, 
                            mutation_prior_sd)    
    }else{
      parameter <- new(ROCParameter, restart.file)
    }
  }else if(model == "FONSE"){
    if(is.null(restart.file)){
      parameter <- initializeFONSEParameterObject(genome, sphi, numMixtures, 
                            geneAssignment, expressionValues, split.serine, 
                            mixture.definition, mixture.definition.matrix)
    }else{
      parameter <- new(FONSEParameter, restart.file)
    }
  }else if(model == "RFP"){
    if(is.null(restart.file)){
      parameter <- initializeRFPParameterObject(genome, sphi, numMixtures, 
                            geneAssignment, expressionValues, split.serine, 
                            mixture.definition, mixture.definition.matrix) 
    }else{
      parameter <- new(RFPParameter, restart.file)
    }
  }else{
    stop("Unknown model.")
  }
  
  return(parameter)
}


#Called from initializeParameterObject. 
initializeROCParameterObject <- function(genome, sphi, numMixtures, geneAssignment,
                      expressionValues = NULL, split.serine = TRUE,
                      mixture.definition = "allUnique", 
                      mixture.definition.matrix = NULL, mutation_prior_sd = 0.35){

  if(is.null(mixture.definition.matrix)){ 
    # keyword constructor
    parameter <- new(ROCParameter, sphi, numMixtures, geneAssignment, 
                     split.serine, mixture.definition)
  }else{
    #matrix constructor
    mixture.definition <- c(mixture.definition.matrix[, 1], 
                            mixture.definition.matrix[, 2])
    parameter <- new(ROCParameter, sphi, numMixtures, geneAssignment, 
                     mixture.definition, split.serine)
  }
  
  
  # initialize expression values
  if(is.null(expressionValues)){
    parameter$initializeSynthesisRateByGenome(genome, mean(sphi))
  }else{
    parameter$initializeSynthesisRateByList(expressionValues)
  }
  
  parameter$mutation_prior_sd <- (mutation_prior_sd)
  parameter <- initializeCovarianceMatricies(parameter, genome, numMixtures)
  
  
  return(parameter)
}


#Called from initializeParameterObject.
initializeRFPParameterObject <- function(genome, sphi, numMixtures, geneAssignment, 
                          expressionValues = NULL, split.serine = TRUE, 
                          mixture.definition = "allUnique", 
                          mixture.definition.matrix = NULL){

  if(is.null(mixture.definition.matrix))
  { # keyword constructor
    parameter <- new(RFPParameter, sphi, numMixtures, geneAssignment, 
                     split.serine, mixture.definition)
  }else{
    #matrix constructor
    mixture.definition <- c(mixture.definition.matrix[, 1], 
                            mixture.definition.matrix[, 2])
    parameter <- new(RFPParameter, sphi, numMixtures, geneAssignment, 
                     mixture.definition, split.serine)
  }
  
  
  # initialize expression values
  if(is.null(expressionValues)){
    parameter$initializeSynthesisRateByGenome(genome, sphi)
  }else{
    parameter$initializeSynthesisRateByList(expressionValues)
  }
  
  return (parameter)
}

#Called from initializeParameterObject.
initializeFONSEParameterObject <- function(genome, sphi, numMixtures, 
                        geneAssignment, expressionValues = NULL, split.serine = TRUE,
                        mixture.definition = "allUnique", 
                        mixture.definition.matrix = NULL){

  
  # create parameter object
  if(is.null(mixture.definition.matrix))
  { # keyword constructor
    parameter <- new(FONSEParameter, sphi, numMixtures, geneAssignment, 
                     split.serine, mixture.definition)
  }else{
    #matrix constructor
    mixture.definition <- c(mixture.definition.matrix[, 1], 
                            mixture.definition.matrix[, 2])
    parameter <- new(FONSEParameter, sphi, numMixtures, geneAssignment, 
                     mixture.definition, split.serine)
  }
  
  
  # initialize expression values
  if(is.null(expressionValues)){
    parameter$initializeSynthesisRateByGenome(genome, sphi)
  }else{
    parameter$initializeSynthesisRateByList(expressionValues)
  }
  
  parameter <- initializeCovarianceMatricies(parameter, genome, numMixtures)
  
  return(parameter)
}



#' Write Parameter To CSV File 
#' 
#' @param parameter A parameter object that corrosponds to
#' one of the model types. Valid values are "ROC", "RFP", and
#' "FONSE".
#' 
#' @param filename A filename where the data will be written to.
#' This file should end with a "csv" extension.
#' 
#' @param CSP Tells what codon specific parameter should be written to the file.
#' This will vary between models.
#' 
#' @param mixture Tells which mixture the data should be retrieved from to write.
#' 
#' @param samples The number of samples that should be used when calculating
#' the posteriors.
#' 
#' @return This function has no return value.
#' 
#' @description \code{writeParameterToCSV} will obtain the codon specific
#' parameter data for a given parameter and mixture and write this data
#' to a csv file.
#' 
#' @details \code{writeParameterToCSV} will make the necessary calls
#' to writeXXXParameterToCSV based off of the parameter given.
#' 
writeParameterToCSV <- function(parameter, filename, CSP, mixture, samples){
  UseMethod("writeParameterToCSV", parameter)
}


#Called from writeParameterToCSV
writeParameterToCSV.Rcpp_ROCParameter <- function(parameter, filename=NULL, 
                                            CSP=NULL, mixture = 1, samples = 10){
  names.aa <- aminoAcids()
  Amino_Acid <- c()
  Value <- c()
  Codon <- c()
  Std_Deviation <- c()
  
  
  for(aa in names.aa){
    if(aa == "M" || aa == "W" || aa == "X") next
    codons <- AAToCodon(aa, T)
    
    for(i in 1:length(codons)){
      Amino_Acid <- c(Amino_Acid, aa)
      Codon <- c(Codon, codons[i])
      
      
      if(CSP == "Mutation"){
        Value <- c(Value,parameter$getMutationPosteriorMeanForCodon(mixture, samples, codons[i]))
        Std_Deviation <- c(Std_Deviation, sqrt(parameter$getMutationVarianceForCodon(mixture, samples, codons[i], TRUE)))
      }
      else if(CSP == "Selection"){
        Value <- c(Value,parameter$getSelectionPosteriorMeanForCodon(mixture, samples, codons[i]))
        Std_Deviation <- c(Std_Deviation, sqrt(parameter$getSelectionVarianceForCodon(mixture, samples, codons[i], TRUE)))
      }
      else {
        stop("Unknown Parameter type given")
      }
    }
  }
  
  
  data <- data.frame(Amino_Acid,Codon,Value, Std_Deviation)
  if(is.null(filename))
  {
    print(data)
  }else {
    write.csv(data, file = filename, row.names = FALSE, quote=FALSE)
  }
}

#TODO: implement writeCSV for RFP and FONSE




#' Get Codon Counts For Each Amino Acid 
#' 
#' @param aa A one character representation of an amino acid.
#' 
#' @param genome A genome object from which the counts of each
#' codon can be obtained.
#' 
#' @return codonCounts Returns a matrix storing the codonCounts. 
#' 
#' @description \code{getCodonCountsForAA} returns a matrix filled with 
#' the number of times a codon is seen in each gene.
#' 
#' @details The returned matrix will have the row corrospond to the
#' genes in the genome and the columns corrospond to the codons for the 
#' given aa. The values will the number of times the codon is present in 
#' that gene.
#' 
getCodonCountsForAA <- function(aa, genome){
  # get codon count for aa
  codons <- AAToCodon(aa, F)
  codonCounts <- lapply(codons, function(codon){
    codonCounts <- genome$getCodonCountsPerGene(codon)
  })
  codonCounts <- do.call("cbind", codonCounts)
  return(codonCounts)
}


# Uses a multinomial logistic regression to estimate the codon specific parameters for every category.
# Delta M is the intercept - and Delta eta is the slope of the regression.
# The package VGAM is used to perform the regression.
getCSPbyLogit <- function(codonCounts, phi, coefstart = NULL, x.arg = FALSE, 
                          y.arg = FALSE, qr.arg = FALSE){
  #avoid cases with 0 aa count
  idx <- rowSums(codonCounts) != 0
  
  # performs the regression and returns Delta M and Delta eta as well as other information no used here
  ret <- VGAM::vglm(codonCounts[idx, ] ~ phi[idx],
                    VGAM::multinomial, coefstart = coefstart,
                    x.arg = x.arg, y.arg = y.arg, qr.arg = qr.arg)
  coefficients <- ret@coefficients
  
  ## convert delta.t to delta.eta
  coefficients <- -coefficients
  
  ret <- list(coefficients = coefficients,
              coef.mat = matrix(coefficients, nrow = 2, byrow = TRUE),
              R = ret@R)
  return(ret)
}


#TODO: Need comments explaining what is going on
subMatrices <- function(M, r, c){
  rg <- (row(M) - 1) %/% r + 1
  cg <- (col(M) - 1) %/% c + 1
  rci <- (rg - 1) * max(cg) + cg
  return(rci)
}



#TODO: Need comments explaining what is going on
splitMatrix <- function(M, r, c){
  rci <- subMatrices(M, r, c)
  N <- prod(dim(M)) / r / c
  cv <- lapply(1:N, function(x) M[rci==x])
  
  return(lapply(1:N, function(i) matrix(cv[[i]], nrow = r)))
} 



# Also initializes the mutaiton and selection parameter
initializeCovarianceMatricies <- function(parameter, genome, numMixtures) {
  numMutationCategory <- parameter$numMutationCategories
  numSelectionCategory <- parameter$numSelectionCategories
  
  phi <- parameter$getCurrentSynthesisRateForMixture(1) # phi values are all the same initially
  names.aa <- aminoAcids()
  
  for(aa in names.aa){
    if(aa == "M" || aa == "W" || aa == "X") next
    #should go away when CT is up and running
    
    codonCounts <- getCodonCountsForAA(aa, genome)
    numCodons <- dim(codonCounts)[2] - 1
    #-----------------------------------------
    # TODO WORKS CURRENTLY ONLY FOR ALLUNIQUE!
    #-----------------------------------------
    covmat <- vector("list", numMixtures)
    for(mixElement in 1:numMixtures){    
      idx <- geneAssignment == mixElement
      csp <- getCSPbyLogit(codonCounts[idx, ], phi[idx])
      
      parameter$initMutation(csp$coef.mat[1,], mixElement, aa)
      parameter$initSelection(csp$coef.mat[2,], mixElement, aa)
    }
    
    # One covariance matrix for all mixtures.
    # Currently only variances used.
    compl.covMat <- diag((numMutationCategory + numSelectionCategory) * numCodons) * 0.05
    parameter$initCovarianceMatrix(compl.covMat, aa)
  }
  
  
  
  #for(aa in names.aa){
  # if(aa == "M" || aa == "W" || aa == "X") next
  #should go away when CT is up and running
  
  #codonCounts <- getCodonCountsForAA(aa, genome)
  #numCodons <- dim(codonCounts)[2] - 1
  #-----------------------------------------
  # TODO WORKS CURRENTLY ONLY FOR ALLUNIQUE!
  #-----------------------------------------
  # covmat <- vector("list", numMixtures)
  #for(mixElement in 1:numMixtures){    
  # idx <- geneAssignment == mixElement
  #csp <- getCSPbyLogit(codonCounts[idx, ], phi[idx])
  
  #   parameter$initMutation(csp$coef.mat[1,], mixElement, aa)
  #  parameter$initSelection(csp$coef.mat[2,], mixElement, aa)
  # split matrix into sup matrices (dM and dEta)
  # covmat[[mixElement]] <- splitMatrix(t(csp$R) %*% csp$R, numCodons, numCodons)  # we expect the covariance matrix, but get the decomposition.
  #  }
  # compl.covMat <- matrix(0, ncol = numMixtures * numCodons * 2, nrow = numMixtures * numCodons * 2)
  #matrix.positions <- subMatrices(compl.covMat, numCodons, numCodons)
  
  #compl.seq <- seq(1, dim(compl.covMat)[1], numCodons)
  #mut.seq <- compl.seq[1:(length(compl.seq)/2)]
  #i <- 1
  #for(pos in mut.seq){ 
  # compl.covMat[matrix.positions == matrix.positions[pos, pos]] <- unlist(covmat[[i]][1])
  #  i <- i + 1
  # i <- ifelse(i > numMutationCategory, 1, i)
  #  }
  # sel.seq <- compl.seq[(length(compl.seq)/2 + 1):length(compl.seq)]
  #  i <- 1
  # for(pos in sel.seq){ 
  #  compl.covMat[matrix.positions == matrix.positions[pos, pos]] <- unlist(covmat[[i]][4])
  # i <- i + 1
  #i <- ifelse(i > numMutationCategory, 1, i)
  #}
  
  #ofdiag.seq <- mut.seq + numCodons*numMutationCategory
  #for(i in 1:length(mut.seq)){
  #  compl.covMat[matrix.positions == matrix.positions[mut.seq[i], ofdiag.seq[i]]] <- unlist(covmat[[i]][2])
  # compl.covMat[matrix.positions == matrix.positions[ofdiag.seq[i], mut.seq[i]]] <- unlist(covmat[[i]][3])
  #}
  #for testing - in actuallity this is used, it is currently overwriting 
  #previous steps.
  #compl.covMat <- diag((numMutationCategory + numSelectionCategory) * numCodons) * 0.05
  #compl.covMat / max(compl.covMat)
  #parameter$initCovarianceMatrix(compl.covMat, aa)
  #}
    
  return(parameter)
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
#' call to writeXXXParameterObject based off of the parameter type 
#' given.
#' 
#' @details For example, if a ROCParameter is passed, the the writeParameterObject
#' for the ROCParameter will be called. This allows us to not have an if-else
#' block in the code - making use of the how R handles these situations.
#' 
writeParameterObject <- function(parameter, file)
{
  UseMethod("writeParameterObject", parameter)
}


# extracts traces and parameter information from the base class Parameter
extractBaseInfo <- function(parameter){
  trace <- parameter$getTraceObject()
  stdDevSynthesisRateTraces <- trace$getStdDevSynthesisRateTraces()
  stdDevSynthesisRateAcceptRatTrace <- trace$getStdDevSynthesisRateAcceptanceRatioTrace()
  synthRateTrace <- trace$getSynthesisRateTrace()
  synthAcceptRatTrace <- trace$getSynthesisRateAcceptanceRatioTrace()
  mixAssignTrace <- trace$getMixutreAssignmentTrace()
  mixProbTrace <- trace$getMixtureProbabilitiesTrace()
  codonSpecificAcceptRatTrace <- trace$getCodonSpecificAcceptanceRatioTrace()
  numMix <- parameter$numMixtures
  numMut <- parameter$numMutationCategories
  numSel <- parameter$numSelectionCategories
  categories <- parameter$getCategories()
  curMixAssignment <- parameter$getMixtureAssignment()
  lastIteration <- parameter$getLastIteration()
  
  varList <- list(stdDevSynthesisRateTraces = stdDevSynthesisRateTraces, 
                    stdDevSynthesisRateAcceptRatTrace = stdDevSynthesisRateAcceptRatTrace,
                    synthRateTrace = synthRateTrace,
                    synthAcceptRatTrace = synthAcceptRatTrace,
                    mixAssignTrace = mixAssignTrace,
                    mixProbTrace = mixProbTrace,
                    codonSpecificAcceptRatTrace = codonSpecificAcceptRatTrace,
                    numMix = numMix,
                    numMut = numMut,
                    numSel = numSel,
                    categories = categories,
                    curMixAssignment = curMixAssignment,
                    lastIteration = lastIteration
                    )
  return(varList)
}


#called from "writeParameterObject."
writeParameterObject.Rcpp_ROCParameter <- function(parameter, file){
  paramBase <- extractBaseInfo(parameter)
  
  currentMutation <- parameter$currentMutationParameter
  currentSelection <- parameter$currentSelectionParameter
  proposedMutation <- parameter$proposedMutationParameter
  proposedSelection <- parameter$proposedSelectionParameter
  model = "ROC"
  mutationPrior <- parameter$getMutationPriorStandardDeviation()
  
  trace <- parameter$getTraceObject()
  
  mutationTrace <- trace$getCodonSpecificParameterTrace(0)
  selectionTrace <- trace$getCodonSpecificParameterTrace(1)
  synthesisOffsetAcceptRatTrace <- trace$getSynthesisOffsetAcceptanceRatioTrace()
  synthesisOffsetTrace <- trace$getSynthesisOffsetTrace()
  observedSynthesisNoiseTrace <- trace$getObservedSynthesisNoiseTrace()
  
  save(list = c("paramBase", "currentMutation", "currentSelection",
                "proposedMutation", "proposedSelection", "model",  
                "mutationPrior", "mutationTrace", "selectionTrace", 
                "synthesisOffsetAcceptRatTrace", "synthesisOffsetTrace", 
                "observedSynthesisNoiseTrace"),
       file=file)
}


#called from "writeParameterObject."
writeParameterObject.Rcpp_RFPParameter <- function(parameter, file){
  paramBase <- extractBaseInfo(parameter)
  
  currentAlpha <- parameter$currentAlphaParameter
  currentLambdaPrime <- parameter$currentLambdaPrimeParameter
  proposedAlpha <- parameter$proposedAlphaParameter
  proposedLambdaPrime <- parameter$proposedLambdaPrimeParameter
  model = "RFP"
  
  
  trace <- parameter$getTraceObject()
  alphaTrace <- trace$getCodonSpecificParameterTrace(0)
  lambdaPrimeTrace <- trace$getCodonSpecificParameterTrace(1)

  save(list = c("paramBase", "currentAlpha", "currentLambdaPrime", "proposedAlpha",
                "proposedLambdaPrime", "model", "alphaTrace", "lambdaPrimeTrace"),
       file=file)
}


#called from "writeParameterObject."
writeParameterObject.Rcpp_FONSEParameter <- function(parameter, file)
{
  #TODO: implement
}




#' Load Parameter Object
#' 
#' @param parameter A parameter object that corrosponds to
#' one of the model types, such as "ROC", or "FONSE".
#' 
#' @param file A filename that where the data will be stored.
#' 
#' @param model Type of the model. Should corrospond to the parameter type.
#' 
#' @return This function has no return value.
#' 
#' @description \code{loadParameterObject} will call the appropriate followup
#' call to loadXXXParameterObject based off of the parameter type 
#' given.
#' 
#' @details For example, if a ROCParameter is passed, the the loadParameterObject
#' for the ROCParameter will be called. This allows us to not have an if-else
#' block in the code - making use of the how R handles these situations.
#' 
loadParameterObject <- function(file)
{
  #A temporary env is set up to stop R errors.
  tempEnv <- new.env();
  load(file = file, envir = tempEnv)
  model <- tempEnv$model
  
  
  if (model == "ROC"){
    parameter <- new(ROCParameter)
    parameter <- loadROCParameterObject(parameter, file)
  }else if (model == "RFP") {
    parameter <- new(RFPParameter)
    parameter <- loadRFPParameterObject(parameter, file)
  }else if (model == "FONSE") {
    parameter <- new(FONSEParameter)
    parameter <- loadFONSEParameterObject(parameter, file)
  }else{
    stop("File data corrupted")
  }
  return(parameter)
}


#Sets all the common variables in the parameter objects.
setBaseInfo <- function(parameter, file)
{
  tempEnv <- new.env();
  load(file = file, envir = tempEnv)
  parameter$setCategories(tempEnv$paramBase$categories)
  #parameter$setCategoriesForTrace() #TODO: FIX 
  parameter$numMixtures <- tempEnv$paramBase$numMix
  parameter$numMutationCategories <- tempEnv$paramBase$numMut
  parameter$numSelectionCategories <- tempEnv$paramBase$numSel
  parameter$setMixtureAssignment(tempEnv$paramBase$curMixAssignment)
  parameter$setLastIteration(tempEnv$paramBase$lastIteration)
  
  trace <- parameter$getTraceObject()
  trace$setStdDevSynthesisRateTraces(tempEnv$paramBase$stdDevSynthesisRateTraces)
  trace$setStdDevSynthesisRateAcceptanceRatioTrace(tempEnv$paramBase$stdDevSynthesisRateAcceptRatTrace)
  trace$setSynthesisRateTrace(tempEnv$paramBase$synthRateTrace)
  trace$setSynthesisRateAcceptanceRatioTrace(tempEnv$paramBase$synthAcceptRatTrace)
  trace$setMixtureAssignmentTrace(tempEnv$paramBase$mixAssignTrace)
  trace$setMixtureProbabilitiesTrace(tempEnv$paramBase$mixProbTrace)
  trace$setCodonSpecificAcceptanceRatioTrace(tempEnv$paramBase$codonSpecificAcceptRatTrace)
  
  parameter$setTraceObject(trace)
}


#Called from "loadParameterObject."
loadROCParameterObject <- function(parameter, file)
{
  tempEnv <- new.env();
  load(file = file, envir = tempEnv)
  
  
  setBaseInfo(parameter, file)
  parameter$currentMutationParameter <- tempEnv$currentMutation
  parameter$currentSelectionParameter <- tempEnv$currentSelection
  parameter$proposedMutationParameter <- tempEnv$proposedMutation
  parameter$proposedSelectionParameter <- tempEnv$proposedSelection
  trace <- parameter$getTraceObject()
  trace$setSynthesisOffsetTrace(tempEnv$synthesisOffsetTrace)
  trace$setSynthesisOffsetAcceptanceRatioTrace(tempEnv$synthesisOffsetAcceptRatTrace)
  trace$setObservedSynthesisNoiseTrace(tempEnv$observedSynthesisNoiseTrace)
  trace$setCodonSpecificParameterTrace(tempEnv$mutationTrace, 0)
  trace$setCodonSpecificParameterTrace(tempEnv$selectionTrace, 1)
  
  
  parameter$setTraceObject(trace)
  return(parameter) 
}


#Called from "loadParameterObject."
loadRFPParameterObject <- function(parameter, file)
{
  tempEnv <- new.env();
  load(file = file, envir = tempEnv)
  
  setBaseInfo(parameter, file)
  parameter$currentAlphaParameter <- tempEnv$currentAlpha
  parameter$proposedAlphaParameter <- tempEnv$proposedAlpha
  parameter$currentLambdaPrimeParameter <- tempEnv$currentLambdaPrime
  parameter$proposedLambdaPrimeParameter <- tempEnv$proposedLambdaPrime
  
  trace <- parameter$getTraceObject()
  trace$setCodonSpecificParameterTrace(tempEnv$alphaTrace, 0)
  trace$setCodonSpecificParameterTrace(tempEnv$lambdaPrimeTrace, 1)
  
  parameter$setTraceObject(trace)
  return(parameter) 
}


#Called from "loadParameterObject."
loadFONSEParameterObject <- function(parameter, file)
{
 #TODO: 
}