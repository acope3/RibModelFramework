#'  Model Initialization 
#'
#' @param parameter An object created with \code{initializeParameterObject}. 
#'  
#' @param model A string containing the model to run (ROC, FONSE, or RFP), has to match parameter object. 
#'  
#' @param with.phi (ROC only) A boolean that determines whether or not to include empirical
#'    phi values (expression rates) for the calculations. 
#'    
#' @param fix.observation.noise (ROC only) Allows to fix the noise in the observed expression dataset to the initial condition.
#'	The initial condition for the observed expression noise can be set in the parameter object. 
#'  
#' @return This function returns the model object created. 
#'  
#' @description initializes the model object. 
#' 
#' @details initializeModelObject initializes a model. The type of model is determined based on the string passed to the \code{model} argument.
#'  The Parameter object has to match the model that is initialized. E.g. to initialize a ROC model, 
#'  it is required that a ROC parameter object is passed to the function.
#'        
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
