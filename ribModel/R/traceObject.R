# NOT EXPOSED
convergence.test.Rcpp_ROCTrace <- function(object, nsamples = 10, frac1 = 0.1, frac2 = 0.5, thin = 1, plot = FALSE, what, mixture = 1, ...)
{
  # TODO: extend to work with multiple chains once we have that capability.
  current.trace <- 0
  
  if(what[1] == "Mutation" || what[1] == "Selection")
  {
    names.aa <- aminoAcids()
    numCodons <- 0
    for(aa in names.aa) {numCodons <- numCodons + length(codons)}
    
    index <- 1
    cur.trace <- vector("list", numCodons)
    for(aa in names.aa)
    {
      codons <- AAToCodon(aa, T)
      for(i in 1:length(codons))
      {
        if(what[1] == "Mutation"){
          cur.trace[[index]] <- object$getMutationParameterTraceByMixtureElementForCodon(mixture, codons[i])
        }else{
          cur.trace[[index]] <- object$getSelectionParameterTraceByMixtureElementForCodon(mixture, codons[i])
        }
        index <- index + 1
      }
    }
    current.trace <- do.call("rbind", cur.trace)
  }
 
  if(what[1] == "MixtureProbability")
  {
    numMixtures <- object$getNumberOfMixtures()
    cur.trace <- vector("list", numMixtures)
    for(i in 1:numMixtures)
    {
      cur.trace[[i]] <- object$getMixtureProbabilitiesTraceForMixture(i)
    }
    current.trace <- do.call("rbind", cur.trace)
  }
  if(what[1] == "Sphi")
  {
    current.trace <- object$getSPhiTrace()
  }
  if(what[1] == "Mphi") 
  {
    sphi <- object$getSPhiTrace();
    mphi <- -(sphi * sphi) / 2;
    current.trace <- mphi
  }
  if(what[1] == "Aphi")
  {
    # TODO need way to determine number of Aphi traces
  }
  if(what[1] == "Sepsilon") 
  {
    # TODO need way to determine number of Sepsilon traces
  }
  if(what[1] == "ExpectedPhi")
  {
    current.trace <- trace$getExpectedPhiTrace()
  }
  if(what[1] == "Expression")
  {
    # TODO need way to determine number of expression traces
  } 

  trace.length <- length(current.trace)
  start <- max(0, trace.length - nsamples)
  
  mcmcobj <- coda::mcmc(data=current.trace, start=start, thin=thin)
  diag <- coda::geweke.diag(mcmcobj, frac1=frac1, frac2=frac2)
  if(plot){ 
    coda::geweke.plot(diag, frac1=frac1, frac2=frac2, ...)
  }else{
    return(diag)
  }
}