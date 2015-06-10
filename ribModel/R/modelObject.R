
initializeModelObject <- function(model = "ROC")
{  
  if(model == "ROC")
  {
    c.model <- new(ROCModel)
  }else if(model == "NSE"){
    cat("MODEL NOT IMPLEMENTED")
  }else{
    stop("Unkown model.")
  }
  return(c.model)
}