
initializeModelObject <- function(parameter, model = "ROC")
{  
  if(model == "ROC")
  {
    c.model <- new(ROCModel)
  }else if(model == "FONSE"){
    c.model = new(FONSEModel)
  }else if (model == "RFP"){
    c.model <- new(RFPModel)
  }else{
    stop("Unknown model.")
  }
  c.model$setParameter(parameter)
  return(c.model)
}