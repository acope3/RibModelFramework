
convergence.test <- function(trace, what, mixture = 1, n.samples = 10, frac1 = 0.1, frac2 = 0.5, plot = FALSE)
{
  UseMethod("convergence.test", trace)
}
convergence.test.Rcpp_ROCTrace <- function(trace, what, mixture, n.samples, frac1, frac2, plot)
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
          cur.trace[[index]] <- trace$getMutationParameterTraceByMixtureElementForCodon(mixture, codons[i])
        }else{
          cur.trace[[index]] <- trace$getSelectionParameterTraceByMixtureElementForCodon(mixture, codons[i])
        }
        index <- index + 1
      }
    }
    current.trace <- do.call("rbind", cur.trace)
  }
 
  if(what[1] == "MixtureProbability")
  {
    numMixtures <- trace$getNumberOfMixtures()
    cur.trace <- vector("list", numMixtures)
    for(i in 1:numMixtures)
    {
      cur.trace[[i]] <- trace$getMixtureProbabilitiesTraceForMixture(i)
    }
    current.trace <- do.call("rbind", cur.trace)
  }
  if(what[1] == "Sphi")
  {
    current.trace <- trace$getSPhiTrace()
  }
  if(what[1] == "Mphi") 
  {
    sphi <- trace$getSPhiTrace();
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
  
  mcmcobj <- mcmc(data=current.trace, start=start, thin=thin)
  diag <- geweke.diag(mcmcobj, frac1=frac1, frac2=frac2)
  if(plot){ 
    geweke.plot(diag, frac1=frac1, frac2=frac2)
  }else{
    return(diag)
  }
}