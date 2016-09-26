#ifndef STANDALONE
#include "include/ROC/ROCModel.h"
#include "include/RFP/RFPModel.h"
#include "include/FONSE/FONSEModel.h"
#include <Rcpp.h>
using namespace Rcpp;

/*
void roc_finalizer(ROCModel* m)
{
	delete m;
}
*/

RCPP_EXPOSED_CLASS(ROCParameter)
RCPP_EXPOSED_CLASS(RFPParameter)
RCPP_EXPOSED_CLASS(FONSEParameter)
RCPP_EXPOSED_CLASS(Parameter)
RCPP_EXPOSED_CLASS(Genome)
RCPP_MODULE(Model_mod)
{
	class_<Model>("Model")
		;

	class_<ROCModel>( "ROCModel" )
		.derives<Model>("Model")
		.constructor<bool>()
		//.finalizer(&roc_finalizer)
		.method("CalculateProbabilitiesForCodons", &ROCModel::CalculateProbabilitiesForCodons,
		        "Calculated codon probabilities. Input is one element shorter than output")
  		.method("getParameter", &ROCModel::getParameter)
  		.method("setParameter", &ROCModel::setParameter)
  		.method("simulateGenome", &ROCModel::simulateGenome)
		;
	
	class_<RFPModel>("RFPModel")
		.derives<Model>("Model")
		.constructor()
  		.method("getParameter", &RFPModel::getParameter)
		.method("setParameter", &RFPModel::setParameter)
		.method("simulateGenome", &RFPModel::simulateGenome) //TODO: Debug this. Does NOT work in R (unknown crash).
		;

	class_<FONSEModel>("FONSEModel")
		.derives<Model>("Model")
		.constructor()
  		.method("getParameter", &FONSEModel::getParameter)
		.method("setParameter", &FONSEModel::setParameter)
		.method("simulateGenome", &FONSEModel::simulateGenome)
		;
}
#endif
