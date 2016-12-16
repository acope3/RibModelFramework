#'  Model Initialization
#'  
#'  \code{initializeModelObject} initializes the Rcpp Model Object
#'  
#' @param parameter An Rcpp parameter object
#' @param model A string containing the model to run (ROC, FONSE, or RFP)
#' @param with.phi (ROC only) A boolean that determines whether or not to use empirical
#'    phi values (expression rates) for the calculations.
#' @param fix.observation.noise (ROC only) Allows to fix the noise in the observed expression dataset to the initial condition.
#'	The initial condition for the observed expression noise can be set in the parameter object.
#'    
#' @return This function returns the model object created.
initializeModelObject <- function(parameter, model = "ROC", with.phi = FALSE, fix.observation.noise = FALSE) {  
  if(model == "ROC") {
    c.model <- new(ROCModel, with.phi, fix.observation.noise)
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
