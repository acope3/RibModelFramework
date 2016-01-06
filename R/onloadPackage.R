NAMESPACE <- environment()
.onLoad <- function(libname, pkgname){
  #library.dynam("ribModel", pkgname, libname)

  Rcpp::loadRcppModules()
  
#  invisible()
} # End of .onLoad().

#.onUnload <- function(libpath){
#  library.dynam.unload("ribModel", libpath)
#  invisible()
#} # End of .onUnload().

#.onAttach <- function(libname, pkgname){
#  invisible()
#} # End of .onAttach().
