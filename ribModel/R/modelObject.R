
initializeModelObject <- function(parameter, model = "ROC")
{  
  if(model == "ROC")
  {
    c.model <- new(ROCModel)
  }else if(model == "NSE"){
    cat("MODEL NOT IMPLEMENTED")
  }else if (model == "RFP"){
    c.model <- new(RFPModel)
  }else{
    stop("Unkown model.")
  }
  c.model$setParameter(parameter)
  return(c.model)
}