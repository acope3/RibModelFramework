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
	class_<Parameter>("Parameter")

		//Initialization and Restart Functions:
		.method("initializeSynthesisRateByGenome", &Parameter::initializeSynthesisRateByGenome)
		.method("initializeSynthesisRateByList", &Parameter::initializeSynthesisRateByList)
		.method("initializeSynthesisRateByRandom", &Parameter::initializeSynthesisRateByRandom)
		//checkIndex is not listed/exposed since it is only called from the other R functions
		.method("readPhiValues", &Parameter::readPhiValues) //Not a R wrapper



		//Mixture Definition Matrix and Category Functions:
		.method("getMutationCategoryForMixture", &Parameter::getMutationCategoryForMixture)
		.method("getSelectionCategoryForMixture", &Parameter::getSelectionCategoryForMixture)
		.method("getSynthesisRateCategoryForMixture", &Parameter::getSynthesisRateCategoryForMixture)
		.method("getCategories", &Parameter::getCategories)
        .method("setCategories", &Parameter::setCategories)
		//setNumMutationCategories and setNumSelectionCategories are taken care of in the
		//properties section below.



		//Group List Functions:
		.method("getGroupList", &Parameter::getGroupList)



		//Synthesis Rate Functions:
		.method("getSynthesisRate", &Parameter::getSynthesisRateR)
		.method("getCurrentSynthesisRateForMixture", &Parameter::getCurrentSynthesisRateForMixture)



		//Iteration Functions:
		.method("getLastIteration", &Parameter::getLastIteration) //Not a R wrapper
		.method("setLastIteration", &Parameter::setLastIteration) //Not a R wrapper



		//Posterior, Variance, and Estimates Functions:
		.method("getSynthesisRatePosteriorMeanByMixtureElementForGene", &Parameter::getSynthesisRatePosteriorMeanByMixtureElementForGene)
		.method("getSynthesisRateVarianceByMixtureElementForGene", &Parameter::getSynthesisRateVarianceByMixtureElementForGene)
		.method("getEstimatedMixtureAssignmentForGene", &Parameter::getEstimatedMixtureAssignmentForGene, "returns the mixture assignment for a given gene")
		.method("getEstimatedMixtureAssignmentProbabilitiesForGene", &Parameter::getEstimatedMixtureAssignmentProbabilitiesForGene, "returns the probabilities assignment for a given gene")



		//Other Functions:
		.method("getMixtureAssignment", &Parameter::getMixtureAssignmentR)
		.method("setMixtureAssignment", &Parameter::setMixtureAssignmentR)
		.method("getMixtureAssignmentForGene", &Parameter::getMixtureAssignmentForGeneR)
		.method("setMixtureAssignmentForGene", &Parameter::setMixtureAssignmentForGene)
		//setNumMixtureElements it taken care in the properties section below



		//Used for getters and setters
		.property("numMutationCategories", &Parameter::getNumMutationCategories, &Parameter::setNumMutationCategories)
		.property("numSelectionCategories", &Parameter::getNumSelectionCategories, &Parameter::setNumSelectionCategories)
		.property("numMixtures", &Parameter::getNumMixtureElements, &Parameter::setNumMixtureElements)
		;


	class_<ROCParameter>( "ROCParameter" )
		.derives<Parameter>("Parameter")


		//Constructors & Destructors:
        .constructor()
		.constructor <std::string>()
		.constructor <std::vector<double>, std::vector<unsigned>, std::vector<unsigned>, bool>()
		.constructor <std::vector<double>, unsigned, std::vector<unsigned>, bool, std::string>()


		//Initialization, Restart, Index Checking:
		.method("initCovarianceMatrix", &ROCParameter::initCovarianceMatrix)
		.method("initMutationCategories", &ROCParameter::initMutationCategories) //Not an R wrapper
		.method("initSelectionCategories", &ROCParameter::initSelectionCategories) //Not an R wrapper
		.method("getCovarianceMatrixForAA", &ROCParameter::getCovarianceMatrixForAA) //Not an R wrapper
		.method("initSelection", &ROCParameter::initSelection)
		.method("initMutation", &ROCParameter::initMutation)



		//Trace Functions:
		.method("getTraceObject", &ROCParameter::getTraceObject) //TODO: only used in R?
		.method("setROCTrace", &ROCParameter::setTraceObject)
		.method("setCategoriesForTrace", &ROCParameter::setCategoriesForTrace)




		//CSP Functions:
		//Listed in the properties section below. NOTE: these getter/setters are ONLY
		//used in R




		//Posterior, Variance, and Estimates Functions:
		.method("getSphiPosteriorMean", &ROCParameter::getSphiPosteriorMean)
		//TODO: should we have a posterior mean for synthesis rate?
		.method("getMutationPosteriorMeanForCodon", &ROCParameter::getMutationPosteriorMeanForCodon)
		.method("getSelectionPosteriorMeanForCodon", &ROCParameter::getSelectionPosteriorMeanForCodon)
		.method("getSphiVariance", &ROCParameter::getSphiVariance)
		//TODO: should we have a variance for synthesis rate?
		.method("getAphiVariance", &ROCParameter::getAphiVariance) //TODO: this is not wrapped! May not run correctly
		.method("getMutationVarianceForCodon", &ROCParameter::getMutationVarianceForCodon)
		.method("getSelectionVarianceForCodon", &ROCParameter::getSelectionVarianceForCodon)



		//Other Functions:
		.method("calculateSelectionCoefficients", &ROCParameter::calculateSelectionCoefficientsR)


		.property("proposedMutationParameter", &ROCParameter::getProposedMutationParameter, &ROCParameter::setProposedMutationParameter) //R Specific
		.property("proposedSelectionParameter", &ROCParameter::getProposedSelectionParameter, &ROCParameter::setProposedSelectionParameter) //R Specific
		.property("currentMutationParameter", &ROCParameter::getCurrentMutationParameter, &ROCParameter::setCurrentMutationParameter) //R Specific
		.property("currentSelectionParameter", &ROCParameter::getCurrentSelectionParameter, &ROCParameter::setCurrentSelectionParameter) //R Specific
		.property("mutation_prior_sd", &ROCParameter::getMutationPriorStandardDeviation, &ROCParameter::setMutationPriorStandardDeviation)
		;


	class_<RFPParameter>("RFPParameter")
		.derives<Parameter>("Parameter")



		//Constructors & Destructors:
        .constructor()
		.constructor <std::string>()
		.constructor <std::vector<double>, std::vector<unsigned>, std::vector<unsigned>, bool>()
		.constructor <std::vector<double>, unsigned, std::vector<unsigned>, bool, std::string>()



		//Initialization, Restart, Index Checking:
		.method("initAlpha", &RFPParameter::initAlphaR)
		.method("initLambdaPrime", &RFPParameter::initLambdaPrimeR)
		.method("initMutationSelectionCategories", &RFPParameter::initMutationSelectionCategoriesR)


		//Trace Functions:
		.method("getTraceObject", &RFPParameter::getTraceObject) //TODO: is this only used in R?
		.method("setRFPTrace", &RFPParameter::setTraceObject)
        .method("setCategoriesForTrace", &RFPParameter::setCategoriesForTrace)



		//CSP Functions:
		//Listed in the properties section below. NOTE: these getter/setters are ONLY
		//used in R



		//Posterior, Variance, and Estimates Functions:
		.method("getAlphaPosteriorMeanForCodon", &RFPParameter::getAlphaPosteriorMeanForCodon)
		.method("getLambdaPrimePosteriorMeanForCodon", &RFPParameter::getLambdaPrimePosteriorMeanForCodon)
		.method("getAlphaVarianceForCodon", &RFPParameter::getAlphaVarianceForCodon)
		.method("getLambdaPrimeVarianceForCodon", &RFPParameter::getLambdaPrimeVarianceForCodon)



		//Other Functions:
		.method("getParameterForCategory", &RFPParameter::getParameterForCategoryR)



		.property("proposedAlphaParameter", &RFPParameter::getProposedAlphaParameter, &RFPParameter::setProposedAlphaParameter) //R Specific
		.property("proposedLambdaPrimeParameter", &RFPParameter::getProposedLambdaPrimeParameter, &RFPParameter::setProposedLambdaPrimeParameter) //R Specific
		.property("currentAlphaParameter", &RFPParameter::getCurrentAlphaParameter, &RFPParameter::setCurrentAlphaParameter) //R Specific
		.property("currentLambdaPrimeParameter", &RFPParameter::getCurrentLambdaPrimeParameter, &RFPParameter::setCurrentLambdaPrimeParameter) //R Specific
		;


	class_<FONSEParameter>("FONSEParameter")
		.derives<Parameter>("Parameter")



		//Constructors & Destructors:
		.constructor <std::string>()
		.constructor <std::vector<double>, std::vector<unsigned>, std::vector<unsigned>, bool>()
		.constructor <std::vector<double>, unsigned, std::vector<unsigned>, bool, std::string>()



		//Initialization, Restart, Index Checking:
		.method("initCovarianceMatrix", &FONSEParameter::initCovarianceMatrix)
		.method("getCovarianceMatrixForAA", &FONSEParameter::getCovarianceMatrixForAA) //Not an R wrapper
		.method("initMutation", &FONSEParameter::initMutation)
		.method("initSelection", &FONSEParameter::initSelection)




		//Trace Functions:
		.method("getTraceObject", &FONSEParameter::getTraceObject) //TODO: only used in R?
        .method("setFONSETrace", &FONSEParameter::setTraceObject)
        .method("setCategoriesForTrace", &FONSEParameter::setCategoriesForTrace)



		//CSP Functions:




		//Posterior, Variance, and Estimates Functions:
		.method("getSphiPosteriorMean", &FONSEParameter::getSphiPosteriorMean)
		.method("getMutationPosteriorMeanForCodon", &FONSEParameter::getMutationPosteriorMeanForCodon)
		.method("getSelectionPosteriorMeanForCodon", &FONSEParameter::getSelectionPosteriorMeanForCodon)
		.method("getSphiVariance", &FONSEParameter::getSphiVariance)
		.method("getMutationVarianceForCodon", &FONSEParameter::getMutationVarianceForCodon)
		.method("getSelectionVarianceForCodon", &FONSEParameter::getSelectionVarianceForCodon)


		//Other Functions:
		.method("calculateSelectionCoefficients", &FONSEParameter::calculateSelectionCoefficientsR)
		;
}
#endif
