#'  Model Initialization
#'  
#'  \code{initializeModelObject} initializes the Rcpp Model Object
#'  
#'  @param parameter An Rcpp parameter object
#'  @param model A string containing the model to run (currently ROC, FONSE, or 
#'    RFP)
#'  @param with.phi A boolean that determines whether or not to use empirical
#'    phi values (expression rates) for the calculations. Currently, this is
#'    only implemented for the ROC model.
#'    
#'  @return This function returns the model object created.
initializeModelObject <- function(parameter, model = "ROC", with.phi = FALSE) {  
  if(model == "ROC") {
    c.model <- new(ROCModel, with.phi)
  } else if(model == "FONSE") {
    c.model = new(FONSEModel)
  } else if (model == "RFP") {
    c.model <- new(RFPModel)
  } else {
    stop("Unknown model.")
  }
  c.model$setParameter(parameter)
  return(c.model)
}