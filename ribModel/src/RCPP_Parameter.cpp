#ifndef STANDALONE
#include "include/Parameter.h"
#include "include/ROCParameter.h"
#include <Rcpp.h>
using namespace Rcpp;

RCPP_EXPOSED_CLASS(ROCTrace)
RCPP_EXPOSED_CLASS(Genome)
RCPP_EXPOSED_CLASS(CovarianceMatrix)

RCPP_MODULE(Parameter_mod)
{
	class_<Parameter>(	"Parameter" )
		.property("numMutationCategories", &Parameter::getNumMutationCategories)
		.property("numSelectionCategories", &Parameter::getNumSelectionCategories)
		.property("numMixtures", &Parameter::getNumMixtureElements)
		;
	class_<ROCParameter>( "ROCParameter" )
		.derives<Parameter>("Parameter")
		.constructor("empty constructor")
		.constructor <double, std::vector<unsigned>, std::vector<unsigned>, bool>()
		.constructor <double, unsigned, std::vector<unsigned>, bool, std::string>()
		.method("initMutationSelectionCategories", &ROCParameter::initMutationSelectionCategoriesR)
		.method("readPhiValues", &ROCParameter::readPhiValues)
		.method("setMixtureAssignmentForGene", &ROCParameter::setMixtureAssignmentForGene)
		.method("getMutationCategoryForMixture", &ROCParameter::getMutationCategoryForMixture)
		.method("getSelectionCategoryForMixture", &ROCParameter::getSelectionCategoryForMixture)
		.method("getSynthesisRateCategoryForMixture", &ROCParameter::getSynthesisRateCategoryForMixture)

		.method("initSelection", &ROCParameter::initSelection)
		.method("initMutation", &ROCParameter::initMutation)
		.method("initCovarianceMatrix", &ROCParameter::initCovarianceMatrix)
		.method("getCovarianceMatrixForAA", &ROCParameter::getCovarianceMatrixForAA)
		.method("getTraceObject", &ROCParameter::getTraceObject)

		.method("initSelection", &ROCParameter::initSelection)
		.method("initMutation", &ROCParameter::initMutation)
		.method("initCovarianceMatrix", &ROCParameter::initCovarianceMatrix)
		.method("getCovarianceMatrixForAA", &ROCParameter::getCovarianceMatrixForAA)
		.method("getTraceObject", &ROCParameter::getTraceObject)

		//R wrapper functions
		.method("initializeSynthesisRateByGenome", &ROCParameter::initializeSynthesisRateByGenome)
		.method("initializeSynthesisRateByList", &ROCParameter::initializeSynthesisRateByList)
		.method("initializeSynthesisRateByRandom", &ROCParameter::initializeSynthesisRateByRandom)
		.method("getEstimatedMixtureAssignmentForGene", &ROCParameter::getEstimatedMixtureAssignmentForGene, "returns the mixture assignment for a given gene")
		.method("getEstimatedMixtureAssignmentProbabilitiesForGene", &ROCParameter::getEstimatedMixtureAssignmentProbabilitiesForGene, "returns the probabilities assignment for a given gene")
		.method("calculateSelectionCoefficients", &ROCParameter::calculateSelectionCoefficientsR)
		// Posterior functions
		.method("getSphiPosteriorMean", &ROCParameter::getSphiPosteriorMean)
		//.method("getMixtureAssignmentPosteriorMean", &ROCParameter::getMixtureAssignmentPosteriorMeanR)

		//R wrapper functions
		.method("getSynthesisRatePosteriorMeanByMixtureElementForGene", &ROCParameter::getSynthesisRatePosteriorMeanByMixtureElementForGene)
		.method("getMutationPosteriorMeanForCodon", &ROCParameter::getMutationPosteriorMeanForCodon)
		.method("getSelectionPosteriorMeanForCodon", &ROCParameter::getSelectionPosteriorMeanForCodon)

		// Variance functions
		.method("getSphiVariance", &ROCParameter::getSphiVariance)

		//R wrapper functions
		.method("getMutationVarianceForCodon", &ROCParameter::getMutationVarianceForCodon)
		.method("getSelectionVarianceForCodon", &ROCParameter::getSelectionVarianceForCodon)
		.method("getSynthesisRateVarianceByMixtureElementForGene", &ROCParameter::getSynthesisRateVarianceByMixtureElementForGene)
		.method("getCurrentSynthesisRateForMixture", &ROCParameter::getCurrentSynthesisRateForMixture)

		;
}
#endif
