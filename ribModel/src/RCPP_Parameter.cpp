#ifndef STANDALONE
#include "include/base/Parameter.h"
#include "include/ROC/ROCParameter.h"
#include "include/RFP/RFPParameter.h"
#include <Rcpp.h>
using namespace Rcpp;

RCPP_EXPOSED_CLASS(ROCTrace)
RCPP_EXPOSED_CLASS(RFPTrace)
RCPP_EXPOSED_CLASS(Genome)
RCPP_EXPOSED_CLASS(CovarianceMatrix)

RCPP_MODULE(Parameter_mod)
{
	class_<Parameter>(	"Parameter" )
		.method("readPhiValues", &Parameter::readPhiValues)
		.method("setMixtureAssignmentForGene", &Parameter::setMixtureAssignmentForGene)
		.method("getMutationCategoryForMixture", &Parameter::getMutationCategoryForMixture)
		.method("getSelectionCategoryForMixture", &Parameter::getSelectionCategoryForMixture)
		.method("getSynthesisRateCategoryForMixture", &Parameter::getSynthesisRateCategoryForMixture)
		.method("initCovarianceMatrix", &Parameter::initCovarianceMatrix)
		.method("getCovarianceMatrixForAA", &Parameter::getCovarianceMatrixForAA)
		.method("initializeSynthesisRateByGenome", &Parameter::initializeSynthesisRateByGenome)
		.method("initializeSynthesisRateByList", &Parameter::initializeSynthesisRateByList)
		.method("initializeSynthesisRateByRandom", &Parameter::initializeSynthesisRateByRandom)
		.method("getEstimatedMixtureAssignmentForGene", &Parameter::getEstimatedMixtureAssignmentForGene, "returns the mixture assignment for a given gene")
		.method("getEstimatedMixtureAssignmentProbabilitiesForGene", &Parameter::getEstimatedMixtureAssignmentProbabilitiesForGene, "returns the probabilities assignment for a given gene")
		.method("getSynthesisRatePosteriorMeanByMixtureElementForGene", &Parameter::getSynthesisRatePosteriorMeanByMixtureElementForGene)
		.method("getSynthesisRateVarianceByMixtureElementForGene", &Parameter::getSynthesisRateVarianceByMixtureElementForGene)
		.method("getCurrentSynthesisRateForMixture", &Parameter::getCurrentSynthesisRateForMixture)
		

		.property("numMutationCategories", &Parameter::getNumMutationCategories)
		.property("numSelectionCategories", &Parameter::getNumSelectionCategories)
		.property("numMixtures", &Parameter::getNumMixtureElements)
		;


	class_<ROCParameter>( "ROCParameter" )
		.derives<Parameter>("Parameter")
		.constructor <std::string>()
		.constructor <double, std::vector<unsigned>, std::vector<unsigned>, bool>()
		.constructor <double, unsigned, std::vector<unsigned>, bool, std::string>()
		.method("initMutationSelectionCategories", &ROCParameter::initMutationSelectionCategoriesR)
		.method("initSelection", &ROCParameter::initSelection)
		.method("initMutation", &ROCParameter::initMutation)
		.method("getTraceObject", &ROCParameter::getTraceObject)

		//R wrapper functions
		.method("calculateSelectionCoefficients", &ROCParameter::calculateSelectionCoefficientsR)
		

		// Posterior functions
		.method("getSphiPosteriorMean", &ROCParameter::getSphiPosteriorMean)

		//R wrapper functions
		.method("getMutationPosteriorMeanForCodon", &ROCParameter::getMutationPosteriorMeanForCodon)
		.method("getSelectionPosteriorMeanForCodon", &ROCParameter::getSelectionPosteriorMeanForCodon)

		// Variance functions
		.method("getSphiVariance", &ROCParameter::getSphiVariance)

		//R wrapper functions
		.method("getMutationVarianceForCodon", &ROCParameter::getMutationVarianceForCodon)
		.method("getSelectionVarianceForCodon", &ROCParameter::getSelectionVarianceForCodon)

		;

	class_<RFPParameter>( "RFPParameter" )
		.derives<Parameter>("Parameter")
		.constructor <double, std::vector<unsigned>, std::vector<unsigned>, bool>()
		.constructor <double, unsigned, std::vector<unsigned>, bool, std::string>()
		
		.method("getTraceObject", &RFPParameter::getTraceObject)
		.method("getParameterForCategory", &RFPParameter::getParameterForCategoryR) //need to implement the R wrapper
		;
}
#endif
