#ifndef STANDALONE
#include "include/base/Parameter.h"
#include "include/ROC/ROCParameter.h"
#include "include/RFP/RFPParameter.h"
#include "include/FONSE/FONSEParameter.h"
#include <Rcpp.h>
using namespace Rcpp;

RCPP_EXPOSED_CLASS(ROCTrace)
RCPP_EXPOSED_CLASS(RFPTrace)
RCPP_EXPOSED_CLASS(FONSETrace)
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

		.method("initializeSynthesisRateByGenome", &Parameter::initializeSynthesisRateByGenome)
		.method("initializeSynthesisRateByList", &Parameter::initializeSynthesisRateByList)
		.method("initializeSynthesisRateByRandom", &Parameter::initializeSynthesisRateByRandom)
		.method("getEstimatedMixtureAssignmentForGene", &Parameter::getEstimatedMixtureAssignmentForGene, "returns the mixture assignment for a given gene")
		.method("getEstimatedMixtureAssignmentProbabilitiesForGene", &Parameter::getEstimatedMixtureAssignmentProbabilitiesForGene, "returns the probabilities assignment for a given gene")
		.method("getSynthesisRatePosteriorMeanByMixtureElementForGene", &Parameter::getSynthesisRatePosteriorMeanByMixtureElementForGene)
		.method("getSynthesisRateVarianceByMixtureElementForGene", &Parameter::getSynthesisRateVarianceByMixtureElementForGene)
		.method("getCurrentSynthesisRateForMixture", &Parameter::getCurrentSynthesisRateForMixture)
		.method("getMixtureAssignmentForGene", &Parameter::getMixtureAssignmentForGeneR)
		.method("getLastIteration", &Parameter::getLastIteration)
		.method("getGroupList", &Parameter::getGroupList)	
		.method("getMixtureAssignment", &Parameter::getMixtureAssignmentR)
		.method("getSynthesisRate", &Parameter::getSynthesisRateR)
        .method("getCategories", &Parameter::getCategories)
        .method("setCategories", &Parameter::setCategories)
    


		.property("numMutationCategories", &Parameter::getNumMutationCategories, &Parameter::setNumMutationCategories)
		.property("numSelectionCategories", &Parameter::getNumSelectionCategories, &Parameter::setNumSelectionCategories)
		.property("numMixtures", &Parameter::getNumMixtureElements, &Parameter::setNumMixtureElements)
		;


	class_<ROCParameter>( "ROCParameter" )
		.derives<Parameter>("Parameter")
        .constructor()
		.constructor <std::string>()
		.constructor <std::vector<double>, std::vector<unsigned>, std::vector<unsigned>, bool>()
		.constructor <std::vector<double>, unsigned, std::vector<unsigned>, bool, std::string>()


		.method("initSelection", &ROCParameter::initSelection)
		.method("initMutation", &ROCParameter::initMutation)
		.method("getTraceObject", &ROCParameter::getTraceObject)
		.method("initCovarianceMatrix", &ROCParameter::initCovarianceMatrix)
		.method("getCovarianceMatrixForAA", &ROCParameter::getCovarianceMatrixForAA)
		.method("initMutationCategories", &ROCParameter::initMutationCategories)
		.method("initSelectionCategories", &ROCParameter::initSelectionCategories)

		.method("calculateSelectionCoefficients", &ROCParameter::calculateSelectionCoefficientsR)
		

		.method("getSphiPosteriorMean", &ROCParameter::getSphiPosteriorMean)

		.method("getMutationPosteriorMeanForCodon", &ROCParameter::getMutationPosteriorMeanForCodon)
		.method("getSelectionPosteriorMeanForCodon", &ROCParameter::getSelectionPosteriorMeanForCodon)

		.method("getSphiVariance", &ROCParameter::getSphiVariance)
		.method("getAphiVariance", &ROCParameter::getAphiVariance)

		.method("getMutationPriorStandardDeviation", &ROCParameter::getMutationPriorStandardDeviation)
		.method("setMutationPriorStandardDeviation", &ROCParameter::setMutationPriorStandardDeviation)

		.method("getMutationVarianceForCodon", &ROCParameter::getMutationVarianceForCodon)
		.method("getSelectionVarianceForCodon", &ROCParameter::getSelectionVarianceForCodon)
        .method("setROCTrace", &ROCParameter::setTraceObject)
		.method("setCategoriesForTrace", &ROCParameter::setCategoriesForTrace)
		;

	class_<RFPParameter>( "RFPParameter" )
		.derives<Parameter>("Parameter")
        .constructor()
		.constructor <std::string>()
		.constructor <std::vector<double>, std::vector<unsigned>, std::vector<unsigned>, bool>()
		.constructor <std::vector<double>, unsigned, std::vector<unsigned>, bool, std::string>()

		.method("initAlpha", &RFPParameter::initAlphaR)
		.method("initLambdaPrime", &RFPParameter::initLambdaPrimeR)
		.method("getTraceObject", &RFPParameter::getTraceObject)
		.method("getParameterForCategory", &RFPParameter::getParameterForCategoryR)
		.method("initMutationSelectionCategories", &RFPParameter::initMutationSelectionCategoriesR)
		.method("getAlphaPosteriorMeanForCodon", &RFPParameter::getAlphaPosteriorMeanForCodon)
		.method("getLambdaPrimePosteriorMeanForCodon", &RFPParameter::getLambdaPrimePosteriorMeanForCodon)
		.method("getAlphaVarianceForCodon", &RFPParameter::getAlphaVarianceForCodon)
		.method("getLambdaPrimeVarianceForCodon", &RFPParameter::getLambdaPrimeVarianceForCodon)
		.method("getTmpTrace", &RFPParameter::getTmp)
        .method("getProposedAlphaParameter", &RFPParameter::getProposedAlphaParameter)
        .method("getProposedLambdaPrimeParameter", &RFPParameter::getProposedLambdaPrimeParameter)
        .method("getCurrentAlphaParameter", &RFPParameter::getCurrentAlphaParameter)
        .method("getCurrentLambdaPrimeParameter", &RFPParameter::getCurrentLambdaPrimeParameter)
        .method("setCurrentAlphaParameter", &RFPParameter::setCurrentAlphaParameter)
        .method("setProposedAlphaParameter", &RFPParameter::setProposedAlphaParameter)
        .method("setCurrentLambdaPrimeParameter", &RFPParameter::setCurrentLambdaPrimeParameter)
        .method("setProposedLambdaPrimeParameter", &RFPParameter::setProposedLambdaPrimeParameter)
        .method("setRFPTrace", &RFPParameter::setTraceObject)
        .method("setCategoriesForTrace", &RFPParameter::setCategoriesForTrace)
		;

	class_<FONSEParameter>("FONSEParameter")
		.derives<Parameter>("Parameter")
		.constructor <std::string>()
		.constructor <std::vector<double>, std::vector<unsigned>, std::vector<unsigned>, bool>()
		.constructor <std::vector<double>, unsigned, std::vector<unsigned>, bool, std::string>()
		.method("initMutationSelectionCategories", &FONSEParameter::initMutationSelectionCategoriesR)
		.method("initSelection", &FONSEParameter::initSelection)
		.method("initMutation", &FONSEParameter::initMutation)
		.method("getTraceObject", &FONSEParameter::getTraceObject)
        .method("setFONSETrace", &FONSEParameter::setTraceObject)

		//R wrapper functions
		.method("calculateSelectionCoefficients", &FONSEParameter::calculateSelectionCoefficientsR)


		// Posterior functions
		.method("getSphiPosteriorMean", &FONSEParameter::getSphiPosteriorMean)

		//R wrapper functions
		.method("getMutationPosteriorMeanForCodon", &FONSEParameter::getMutationPosteriorMeanForCodon)
		.method("getSelectionPosteriorMeanForCodon", &FONSEParameter::getSelectionPosteriorMeanForCodon)

		// Variance functions
		.method("getSphiVariance", &FONSEParameter::getSphiVariance)

		//R wrapper functions
		.method("getMutationVarianceForCodon", &FONSEParameter::getMutationVarianceForCodon)
		.method("getSelectionVarianceForCodon", &FONSEParameter::getSelectionVarianceForCodon)
		.method("initCovarianceMatrix", &FONSEParameter::initCovarianceMatrix)
		.method("getCovarianceMatrixForAA", &FONSEParameter::getCovarianceMatrixForAA)
		;
}
#endif
