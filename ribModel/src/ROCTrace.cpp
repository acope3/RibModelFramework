#include "include/ROC/ROCTrace.h"
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

ROCTrace::ROCTrace() : Trace()
{
	//CTOR
}

void ROCTrace::initAllTraces(unsigned samples, unsigned num_genes, unsigned numMutationCategories, unsigned numSelectionCategories,
		unsigned numParam, unsigned numMixtures, std::vector<mixtureDefinition> &_categories, unsigned maxGrouping, unsigned numPhiGroupings)
{
	initBaseTraces(samples, num_genes, numSelectionCategories, numMixtures, _categories, maxGrouping);
	initROCTraces(samples, numMutationCategories, numSelectionCategories, numParam, numPhiGroupings);
}


void ROCTrace::initROCTraces(unsigned samples, unsigned numMutationCategories, unsigned numSelectionCategories, unsigned numParam, unsigned numPhiGroupings)
{
	initMutationParameterTrace(samples, numMutationCategories, numParam);
	initSelectionParameterTrace(samples, numSelectionCategories, numParam);
	initAphiTrace(samples, numPhiGroupings);
	initSepsilonTrace(samples, numPhiGroupings);
}


void ROCTrace::initMutationParameterTrace(unsigned samples, unsigned numMutationCategories, unsigned numParam)
{
	mutationParameterTrace.resize(numMutationCategories);
	for (unsigned category = 0; category < numMutationCategories; category++)
	{
		mutationParameterTrace[category].resize(numParam);
		for (unsigned i = 0; i < numParam; i++)
		{
			mutationParameterTrace[category][i].resize(samples);
			std::vector <double> temp(samples, 0.0);
			mutationParameterTrace[category][i] = temp;
		}
	}
}


void ROCTrace::initSelectionParameterTrace(unsigned samples, unsigned numSelectionCategories, unsigned numParam)
{
	selectionParameterTrace.resize(numSelectionCategories);
	for (unsigned category = 0; category < numSelectionCategories; category++)
	{
		selectionParameterTrace[category].resize(numParam);
		for (unsigned i = 0; i < numParam; i++)
		{
			selectionParameterTrace[category][i].resize(samples);
			std::vector <double> temp(samples, 0.0);
			selectionParameterTrace[category][i] = temp;
		}
	}
}

void ROCTrace::initAphiTrace(unsigned samples, unsigned numPhiGroupings)
{
	AphiTrace.resize(numPhiGroupings);
	for (unsigned i = 0; i < numPhiGroupings; i++) {
		AphiTrace[i].resize(samples);
	}
	
	AphiAcceptanceRatioTrace.resize(numPhiGroupings);
}

void ROCTrace::initSepsilonTrace(unsigned samples, unsigned numPhiGroupings)
{
	SepsilonTrace.resize(numPhiGroupings);
	for (unsigned i = 0; i < numPhiGroupings; i++) {
		SepsilonTrace[i].resize(samples);
	}
}


std::vector<double> ROCTrace::getMutationParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon)
{
	unsigned codonIndex = SequenceSummary::codonToIndex(codon, true);
	unsigned category = getMutationCategory(mixtureElement);
	return mutationParameterTrace[category][codonIndex];
}


std::vector<double> ROCTrace::getSelectionParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon)
{
	unsigned codonIndex = SequenceSummary::codonToIndex(codon, true);
	unsigned category = getSelectionCategory(mixtureElement);
	return selectionParameterTrace[category][codonIndex];
}

void ROCTrace::updateCodonSpecificParameterTrace(unsigned sample, std::string aa, std::vector<std::vector<double>> &curMutParam, std::vector<std::vector<double>> &curSelectParam)
{
	for(unsigned category = 0; category < mutationParameterTrace.size(); category++)
	{
		std::array <unsigned, 2> aaRange = SequenceSummary::AAToCodonRange(aa, true);
		for (unsigned i = aaRange[0]; i < aaRange[1]; i++)
		{
			mutationParameterTrace[category][i][sample] = curMutParam[category][i];
		}
	}
	for(unsigned category = 0; category < selectionParameterTrace.size(); category++)
	{
		std::array <unsigned, 2> aaRange = SequenceSummary::AAToCodonRange(aa, true);
		for (unsigned i = aaRange[0]; i < aaRange[1]; i++)
		{
			selectionParameterTrace[category][i][sample] = curSelectParam[category][i];
		}
	}
}



std::vector<std::vector<std::vector<double>>> ROCTrace::getMutationParameterTrace()
{
	return mutationParameterTrace;
}


std::vector<std::vector<std::vector<double>>> ROCTrace::getSelectionParameterTrace()
{
	return selectionParameterTrace;
}


std::vector<std::vector<double>> ROCTrace::getAphiAcceptanceRatioTrace()
{
	return AphiAcceptanceRatioTrace;
}


//----------------------------------------------------
//----------------------R WRAPPERS--------------------
//----------------------------------------------------

std::vector<double> ROCTrace::getMutationParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon) 
{
	std::vector<double> RV;
	bool checkMixtureElement = checkIndex(mixtureElement, 1, getNumberOfMixtures());
	if (checkMixtureElement)
	{
		RV = getMutationParameterTraceByMixtureElementForCodon(mixtureElement - 1, codon);  
	}
	return RV;
}


std::vector<double> ROCTrace::getSelectionParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon) 
{
	std::vector<double> RV;
	bool checkMixtureElement = checkIndex(mixtureElement, 1, getNumberOfMixtures());
	if (checkMixtureElement)
	{
		RV = getSelectionParameterTraceByMixtureElementForCodon(mixtureElement - 1, codon);  

	}
	return RV;
}

std::vector<std::vector<double>> ROCTrace::getAphiTraceR()
{
	return AphiTrace;
}

std::vector<std::vector<double>> ROCTrace::getSepsilonTraceR()
{
	return SepsilonTrace;
}







void ROCTrace::setAphiTrace(std::vector<std::vector <double> > _AphiTrace)
{
	AphiTrace = _AphiTrace;
}


void ROCTrace::setAphiAcceptanceRatioTrace(std::vector<std::vector <double> > _AphiAcceptanceRatioTrace)
{
	AphiAcceptanceRatioTrace = _AphiAcceptanceRatioTrace;
}


void ROCTrace::setSepsilonTrace(std::vector<std::vector <double> > _SepsilonTrace)
{
	SepsilonTrace = _SepsilonTrace;
}


void ROCTrace::setMutationParameterTrace(std::vector<std::vector<std::vector<double>>> _mutationParameterTrace)
{
	mutationParameterTrace = _mutationParameterTrace;
}


void ROCTrace::setSelectionParameterTrace(std::vector<std::vector<std::vector<double>>> _selectionParameterTrace)
{
	selectionParameterTrace = _selectionParameterTrace;
}
