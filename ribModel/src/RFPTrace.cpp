#include "include/RFP/RFPTrace.h"
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

RFPTrace::RFPTrace() : Trace()
{
	//CTOR
}

void RFPTrace::initAllTraces(unsigned samples, unsigned num_genes, unsigned numMutationCategories, unsigned numSelectionCategories,
		unsigned numParam, unsigned numMixtures, std::vector<mixtureDefinition> &_categories, unsigned maxGrouping)
{
	initBaseTraces(samples, num_genes, numMutationCategories, numMixtures, _categories, maxGrouping);
	initRFPTraces(samples, numMutationCategories, numSelectionCategories, numParam);
}


void RFPTrace::initRFPTraces(unsigned samples, unsigned numMutationCategories, unsigned numSelectionCategories, unsigned numParam)
{
	initAlphaParameterTrace(samples, numMutationCategories, numParam);
	initLambdaPrimeParameterTrace(samples, numSelectionCategories, numParam);
}


void RFPTrace::initAlphaParameterTrace(unsigned samples, unsigned numMutationCategories, unsigned numParam)
{
	alphaParameterTrace.resize(numMutationCategories);
	for (unsigned category = 0; category < numMutationCategories; category++)
	{
		alphaParameterTrace[category].resize(numParam);
		for (unsigned i = 0; i < numParam; i++)
		{
			alphaParameterTrace[category][i].resize(samples);
			std::vector <double> temp(samples, 0.0);
			alphaParameterTrace[category][i] = temp;
		}
	}
}


void RFPTrace::initLambdaPrimeParameterTrace(unsigned samples, unsigned numSelectionCategories, unsigned numParam)
{
	lambdaPrimeParameterTrace.resize(numSelectionCategories);
	for (unsigned category = 0; category < numSelectionCategories; category++)
	{
		lambdaPrimeParameterTrace[category].resize(numParam);
		for (unsigned i = 0; i < numParam; i++)
		{
			lambdaPrimeParameterTrace[category][i].resize(samples);
			std::vector <double> temp(samples, 0.0);
			lambdaPrimeParameterTrace[category][i] = temp;
		}
	}
}


std::vector<double> RFPTrace::getAlphaParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon)
{
	unsigned codonIndex = SequenceSummary::codonToIndex(codon);
	unsigned category = getAlphaCategory(mixtureElement);
	return alphaParameterTrace[category][codonIndex];
}


std::vector<double> RFPTrace::getLambdaPrimeParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon)
{
	unsigned codonIndex = SequenceSummary::codonToIndex(codon);
	unsigned category = getLambdaPrimeCategory(mixtureElement);
	return lambdaPrimeParameterTrace[category][codonIndex];
}


std::vector<std::vector<std::vector<double>>> RFPTrace::getAlphaParameterTrace()
{
	return alphaParameterTrace;
}


std::vector<std::vector<std::vector<double>>> RFPTrace::getLambdaPrimeParameterTrace()
{
	return lambdaPrimeParameterTrace;
}


void RFPTrace::updateCodonSpecificParameterTrace(unsigned sample, std::string codon, std::vector<std::vector<double>> &curAlpParam, std::vector<std::vector<double>> &curLmPriParam)
{
    unsigned i = SequenceSummary::codonToIndex(codon);
	for(unsigned category = 0; category < alphaParameterTrace.size(); category++)
	{
        alphaParameterTrace[category][i][sample] = curAlpParam[category][i];
	}
	for(unsigned category = 0; category < lambdaPrimeParameterTrace.size(); category++)
	{
        lambdaPrimeParameterTrace[category][i][sample] = curLmPriParam[category][i];
	}
}

//----------------------------------------------------
//----------------------R WRAPPERS--------------------
//----------------------------------------------------
std::vector<double> RFPTrace::getAlphaParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon) 
{
	std::vector<double> RV;
	bool checkMixtureElement = checkIndex(mixtureElement, 1, getNumberOfMixtures());
	if (checkMixtureElement)
	{
		RV = getAlphaParameterTraceByMixtureElementForCodon(mixtureElement - 1, codon);  
	}
	return RV;
}


std::vector<double> RFPTrace::getLambdaPrimeParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon) 
{
	std::vector<double> RV;
	bool checkMixtureElement = checkIndex(mixtureElement, 1, getNumberOfMixtures());
	if (checkMixtureElement)
	{
		RV = getLambdaPrimeParameterTraceByMixtureElementForCodon(mixtureElement - 1, codon);  
	}
	return RV;
}


void RFPTrace::setAlphaParameterTrace(std::vector<std::vector<std::vector<double>>> _alphaParameterTrace)
{
    alphaParameterTrace = _alphaParameterTrace;
}

void RFPTrace::setLambdaPrimeParameterTrace(std::vector<std::vector<std::vector<double>>> _lambdaPrimeParameterTrace)
{
    lambdaPrimeParameterTrace = _lambdaPrimeParameterTrace;
}
