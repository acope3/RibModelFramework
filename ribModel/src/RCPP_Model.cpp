#ifndef STANDALONE
#include "include/Model.h"
#include <Rcpp.h>
using namespace Rcpp;
RCPP_EXPOSED_CLASS(ROCParameter)
RCPP_EXPOSED_CLASS(Parameter)
RCPP_EXPOSED_CLASS(Genome)
RCPP_MODULE(Model_mod)
{
	class_<Model>("Model")
		;

  class_<ROCModel>( "ROCModel" )
    .derives<Model>("Model")
		.constructor<ROCParameter>()
    .method("CalculateProbabilitiesForCodons", &ROCModel::CalculateProbabilitiesForCodons, "Calculated codon probabilities. Input is one element shorter than output")
  ;
}
#endif
