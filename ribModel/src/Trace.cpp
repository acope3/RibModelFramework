#include "include/base/Trace.h"
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

Trace::Trace()
{
	//CTOR
	categories = NULL;
}

Trace::~Trace()
{
	//DTOR
}

void Trace::initAllTraces(unsigned samples, unsigned num_genes, unsigned numSelectionCategories, unsigned numMixtures, 
		std::vector<mixtureDefinition> &_categories, unsigned maxGrouping)
{
	initBaseTraces(samples, num_genes, numSelectionCategories, numMixtures, _categories, maxGrouping);
}

void Trace::initBaseTraces(unsigned samples, unsigned num_genes, unsigned numSelectionCategories, unsigned numMixtures, 
		std::vector<mixtureDefinition> &_categories, unsigned maxGrouping)
{
#ifndef STANDALONE
	Rprintf("maxGrouping: %d\n", maxGrouping);
#else
	std::cout << "maxGrouping: " << maxGrouping << "\n";
#endif
	//numSelectionCategories always == numSynthesisRateCategories, so only one is passed in for convience
	initSphiTrace(numSelectionCategories, samples);
	initSynthesisRateAcceptanceRatioTrace(num_genes, numSelectionCategories);
	cspAcceptanceRatioTrace.resize(maxGrouping);
	initSynthesisRateTrace(samples, num_genes, numSelectionCategories);
	initMixtureAssignmentTrace(samples, num_genes);
	initMixtureProbabilitesTrace(samples, numMixtures);

	categories = &_categories;
}


bool Trace::checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound)
{
  bool check = false;
  if (lowerbound <= index && index <= upperbound)
  {
    check = true;
  }
  else
  {
#ifndef STANDALONE
		Rf_error("Index: %d is out of bounds. Index must be between %d & %d\n", index, lowerbound, upperbound);
#else
		std::cerr << "Error with the index\nGIVEN: " << index << "\n";
		std::cerr << "MUST BE BETWEEN:	" << lowerbound << " & " << upperbound << "\n";
#endif
  }

  return check;
}


void Trace::initSphiTrace(unsigned numSelectionCategories, unsigned samples)
{
	sPhiTrace.resize(numSelectionCategories);
	for(unsigned i = 0u; i < numSelectionCategories; i++)
	{
		std::vector<double> temp(samples, 0.0);
		sPhiTrace[i] = temp;
	}
}


void Trace::initSynthesisRateAcceptanceRatioTrace(unsigned num_genes, unsigned numSynthesisRateCategories)
{
	synthesisRateAcceptanceRatioTrace.resize(numSynthesisRateCategories);
	for(unsigned category = 0; category < numSynthesisRateCategories; category++)
	{
		synthesisRateAcceptanceRatioTrace[category].resize(num_genes);
	}
	//NOTE: this is not sized to samples because push_back takes care of the initialization
}


void Trace::initSynthesisRateTrace(unsigned samples, unsigned num_genes, unsigned numSynthesisRateCategories)
{
	synthesisRateTrace.resize(numSynthesisRateCategories);
	for(unsigned category = 0; category < numSynthesisRateCategories; category++)
	{
		synthesisRateTrace[category].resize(num_genes);
		for(unsigned i = 0; i < num_genes; i++)
		{
			synthesisRateTrace[category][i].resize(samples);
			std::vector<double> tempExpr(samples, 0.0);
			synthesisRateTrace[category][i] = tempExpr;
		}
	}
}


void Trace::initMixtureAssignmentTrace(unsigned samples, unsigned num_genes)
{
	mixtureAssignmentTrace.resize(num_genes);
	for (unsigned i = 0u; i < num_genes; i ++)
	{
		mixtureAssignmentTrace[i].resize(samples);
	}
}


void Trace::initMixtureProbabilitesTrace(unsigned samples, unsigned numMixtures)
{
	mixtureProbabilitiesTrace.resize(numMixtures);
	for (unsigned i = 0u; i < numMixtures; i++)
	{
		mixtureProbabilitiesTrace[i].resize(samples, 0.0);
	}
}


std::vector<double> Trace::getExpectedPhiTrace()
{
	unsigned numGenes = synthesisRateTrace[0].size(); //number of genes
	unsigned samples = synthesisRateTrace[0][0].size(); //number of samples
	std::vector<double> RV(samples, 0.0);
	for (unsigned sample = 0; sample < samples; sample++)
	{
		for (unsigned geneIndex = 0; geneIndex < numGenes; geneIndex++)
		{
			unsigned mixtureElement = mixtureAssignmentTrace[geneIndex][sample];
      unsigned category = getSynthesisRateCategory(mixtureElement);
			RV[sample] += synthesisRateTrace[category][geneIndex][sample];
		}
		RV[sample] /= numGenes;
	}
	return RV;
}


std::vector<double> Trace::getSynthesisRateAcceptanceRatioTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex)
{
	unsigned category = getSynthesisRateCategory(mixtureElement);
	return synthesisRateAcceptanceRatioTrace[category][geneIndex];
}


std::vector<std::vector<std::vector <double>>> Trace::getSynthesisRateAcceptanceRatioTrace()
{
	return synthesisRateAcceptanceRatioTrace;
}


std::vector<double> Trace::getCspAcceptanceRatioTraceForAA(std::string aa)
{
	aa[0] = (char) std::toupper(aa[0]);
	unsigned aaIndex = SequenceSummary::aaToIndex.find(aa) -> second;
	return cspAcceptanceRatioTrace[aaIndex];
}


std::vector<std::vector<std::vector<double>>> Trace::getSynthesisRateTrace()
{
	return synthesisRateTrace;
}

std::vector<double> Trace::getSynthesisRateTraceForGene(unsigned geneIndex)
{
	unsigned traceLength = synthesisRateTrace[0][0].size();

	std::vector<double> returnVector(traceLength, 0.0);
	for(unsigned i = 0u; i < traceLength; i++)
	{
		unsigned mixtureElement = mixtureAssignmentTrace[geneIndex][i];
		unsigned category = getSynthesisRateCategory(mixtureElement);
		returnVector[i] =  synthesisRateTrace[category][geneIndex][i];
	}
	return returnVector;
}


std::vector<double> Trace::getSynthesisRateTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex)
{
	unsigned category = getSynthesisRateCategory(mixtureElement);
	return synthesisRateTrace[category][geneIndex];
}


void Trace::updateSynthesisRateAcceptanceRatioTrace(unsigned category, unsigned geneIndex, double acceptanceLevel)
{
	synthesisRateAcceptanceRatioTrace[category][geneIndex].push_back(acceptanceLevel);
}


void Trace::updateSynthesisRateTrace(unsigned sample, unsigned geneIndex, std::vector<std::vector <double>> &currentSynthesisRateLevel)
{
	for(unsigned category = 0; category < synthesisRateTrace.size(); category++)
	{
		synthesisRateTrace[category][geneIndex][sample] = currentSynthesisRateLevel[category][geneIndex];
	}
}


void Trace::updateMixtureProbabilitiesTrace(unsigned samples, std::vector<double> &categoryProbabilities)
{
	for(unsigned category = 0; category < mixtureProbabilitiesTrace.size(); category++)
	{
		mixtureProbabilitiesTrace[category][samples] = categoryProbabilities[category];
	}
}


std::vector<std::vector<unsigned>> Trace::getMixtureAssignmentTrace()
{
	return mixtureAssignmentTrace;
}
std::vector<std::vector<double>> Trace::getMixtureProbabilitiesTrace()
{
	return mixtureProbabilitiesTrace;
}
std::vector<std::vector<double>> Trace::getCspAcceptanceRatioTrace()
{
	return cspAcceptanceRatioTrace;
}
//----------------------------------------------------
//----------------------R WRAPPERS--------------------
//----------------------------------------------------


std::vector<double> Trace::getSynthesisRateAcceptanceRatioTraceByMixtureElementForGeneR(unsigned mixtureElement, unsigned geneIndex)
{
	std::vector<double> RV;
	bool checkGene = checkIndex(geneIndex, 1, synthesisRateAcceptanceRatioTrace.size());
	bool checkMixtureElement = checkIndex(mixtureElement, 1, mixtureProbabilitiesTrace.size());
	if (checkGene && checkMixtureElement)
	{
		RV = getSynthesisRateAcceptanceRatioTraceByMixtureElementForGene(mixtureElement - 1, geneIndex - 1);
	}
	return RV;
}


std::vector<double> Trace::getSynthesisRateTraceForGeneR(unsigned geneIndex)
{
	std::vector<double> RV;
	bool checkGene = checkIndex(geneIndex, 1, synthesisRateTrace[0].size());
	if (checkGene)
	{
		RV = getSynthesisRateTraceForGene(geneIndex - 1);
	}
	return RV;
}


std::vector<double> Trace::getSynthesisRateTraceByMixtureElementForGeneR(unsigned mixtureElement, unsigned geneIndex)
{
	std::vector<double> RV;
	bool checkMixtureElement = checkIndex(mixtureElement, 1, mixtureProbabilitiesTrace.size());
	bool checkGene = checkIndex(geneIndex, 1, synthesisRateTrace[0].size());
	if (checkMixtureElement && checkGene)
	{
		RV = getSynthesisRateTraceByMixtureElementForGene(mixtureElement - 1, geneIndex - 1);
	}
	return RV;
}


std::vector<unsigned> Trace::getMixtureAssignmentTraceForGeneR(unsigned geneIndex)
{
	std::vector <unsigned> RV;
	bool checkGene = checkIndex(geneIndex, 1, mixtureAssignmentTrace[0].size());
	if (checkGene)
	{
		RV = getMixtureAssignmentTraceForGene(geneIndex - 1);
	}
	return RV;
}


std::vector<double> Trace::getMixtureProbabilitiesTraceForMixtureR(unsigned mixtureIndex)
{
	std::vector<double> RV;
	bool check = checkIndex(mixtureIndex, 1, mixtureProbabilitiesTrace.size());
	if (check)
	{
		RV = getMixtureProbabilitiesTraceForMixture(mixtureIndex - 1);
	}
	return RV;
}


void Trace::setSphiTraces(std::vector<std::vector<double>> _sPhiTrace)
{
    sPhiTrace = _sPhiTrace;
}


void Trace::setSphiAcceptanceRatioTrace(std::vector<double> _sphiAcceptanceRatioTrace)
{
    sphiAcceptanceRatioTrace = _sphiAcceptanceRatioTrace;
}


void Trace::setSynthesisRateTrace(std::vector<std::vector<std::vector<double>>> _synthesisRateTrace)
{
    synthesisRateTrace = _synthesisRateTrace;
}


void Trace::setSynthesisRateAcceptanceRatioTrace(std::vector<std::vector<std::vector<double>>>_synthesisRateAcceptanceRatioTrace)
{
    synthesisRateAcceptanceRatioTrace = _synthesisRateAcceptanceRatioTrace;
}


void Trace::setMixtureAssignmentTrace(std::vector<std::vector<unsigned>> _mixtureAssignmentTrace)
{
    mixtureAssignmentTrace = _mixtureAssignmentTrace;
}


void Trace::setMixtureProbabilitiesTrace(std::vector<std::vector<double>> _mixtureProbabilitiesTrace)
{
    mixtureProbabilitiesTrace = _mixtureProbabilitiesTrace;
}


void Trace::setCspAcceptanceRatioTrace(std::vector<std::vector<double>> _cspAcceptanceRatioTrace)
{
    cspAcceptanceRatioTrace = _cspAcceptanceRatioTrace;
}

void Trace::setCategories(std::vector<mixtureDefinition> &_categories)
{
    categories = &_categories;
}
