#include "include/ROCTrace.h"
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

ROCTrace::ROCTrace()
{
	//CTOR
}

void ROCTrace::initAllTraces(unsigned samples, unsigned num_genes, unsigned adaptiveSamples, unsigned numMutationCategories, unsigned numSelectionCategories,
		unsigned numParam, unsigned numMixtures, std::vector<thetaK> &_categories)
{
	//numSelectionCategories always == numExpressionCategories, so only one is passed in for convience
	initSphiTrace(samples);
	initExpressionAcceptanceRatioTrace(samples, num_genes, numSelectionCategories);
	cspAcceptanceRatioTrace.resize(22);
	initExpressionTrace(samples, num_genes, numSelectionCategories);
	initMutationParameterTrace(samples, numMutationCategories, numParam);
	initSelectionParameterTrace(samples, numSelectionCategories, numParam);
	initMixtureAssignmentTrace(samples, num_genes);
	initMixtureProbabilitesTrace(samples, numMixtures);

	categories = &_categories;
}


bool ROCTrace::checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound)
{
  bool check = false;
  if (lowerbound <= index && index <= upperbound)
  {
    check = true;
  }
  else
  {
    std::cerr <<"Error with Index\nINDEX: " << index <<"\n";
    std::cerr <<"MUST BE BETWEEN " << lowerbound << " & " << upperbound <<"\n";
  }

  return check;
}


void ROCTrace::initSphiTrace(unsigned samples)
{
	sPhiTrace.resize(samples);
}


void ROCTrace::initExpressionAcceptanceRatioTrace(unsigned samples, unsigned num_genes, unsigned numExpressionCategories)
{
	expressionAcceptanceRatioTrace.resize(numExpressionCategories);
	for(unsigned category = 0; category < numExpressionCategories; category++)
	{
		expressionAcceptanceRatioTrace[category].resize(num_genes);
	}
}


void ROCTrace::initExpressionTrace(unsigned samples, unsigned num_genes, unsigned numExpressionCategories)
{
	expressionTrace.resize(numExpressionCategories);
	for(unsigned category = 0; category < numExpressionCategories; category++)
	{
		expressionTrace[category].resize(num_genes);
		for(unsigned i = 0; i < num_genes; i++)
		{
			expressionTrace[category][i].resize(samples);
			std::vector<double> tempExpr(samples, 0.0);
			expressionTrace[category][i] = tempExpr;
		}
	}
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


void ROCTrace::initMixtureAssignmentTrace(unsigned samples, unsigned num_genes)
{
	mixtureAssignmentTrace.resize(num_genes);
	for (unsigned i = 0u; i < num_genes; i ++)
	{
		mixtureAssignmentTrace[i].resize(samples);
	}
}


void ROCTrace::initMixtureProbabilitesTrace(unsigned samples, unsigned numMixtures)
{
	mixtureProbabilitiesTrace.resize(numMixtures);
	for (unsigned i = 0u; i < numMixtures; i++)
	{
		mixtureProbabilitiesTrace[i].resize(samples, 0.0);
	}
}


std::vector<double> ROCTrace::getExpectedPhiTrace()
{
	unsigned numGenes = expressionTrace[0].size(); //number of genes
	unsigned samples = expressionTrace[0][0].size(); //number of samples
	std::vector<double> RV(samples, 0.0);
	for (unsigned sample = 0; sample < samples; sample++)
	{
		for (unsigned geneIndex = 0; geneIndex < numGenes; geneIndex++)
		{
			unsigned category = mixtureAssignmentTrace[geneIndex][sample];
			RV[sample] += expressionTrace[category][geneIndex][sample];
		}
		RV[sample] /= numGenes;
	}
	return RV;
}


std::vector<double> ROCTrace::getExpressionAcceptanceRatioTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex)
{
	unsigned category = getExpressionCategory(mixtureElement);
	return expressionAcceptanceRatioTrace[category][geneIndex];
}


std::vector<double> ROCTrace::getCspAcceptanceRatioTraceForAA(char aa)
{
	aa = std::toupper(aa);
	unsigned aaIndex = SequenceSummary:: aaToIndex.find(aa) -> second;
	return cspAcceptanceRatioTrace[aaIndex];
}


std::vector<double> ROCTrace::getExpressionTraceForGene(unsigned geneIndex)
{
	unsigned traceLength = expressionTrace[0][0].size();

	std::vector<double> returnVector(traceLength, 0.0);
	for(unsigned i = 0u; i < traceLength; i++)
	{
		unsigned category = mixtureAssignmentTrace[geneIndex][i];
		returnVector[i] =  expressionTrace[category][geneIndex][i];
	}
	return returnVector;
}


std::vector<double> ROCTrace::getExpressionTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex)
{
	unsigned category = getExpressionCategory(mixtureElement);
	return expressionTrace[category][geneIndex];
}


std::vector<double> ROCTrace::getMutationParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon)
{
	unsigned codonIndex = SequenceSummary::CodonToIndex(codon, true);
	unsigned category = getMutationCategory(mixtureElement);
	return mutationParameterTrace[category][codonIndex];
}


std::vector<double> ROCTrace::getSelectionParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon)
{
	unsigned codonIndex = SequenceSummary::CodonToIndex(codon, true);
	unsigned category = getSelectionCategory(mixtureElement);
	return selectionParameterTrace[category][codonIndex];
}


void ROCTrace::updateExpressionAcceptanceRatioTrace(unsigned category, unsigned geneIndex, double acceptanceLevel)
{
	expressionAcceptanceRatioTrace[category][geneIndex].push_back(acceptanceLevel);
}

void ROCTrace::updateCodonSpecificParameterTrace(unsigned sample, char aa, std::vector<std::vector<double>> &curMutParam, std::vector<std::vector<double>> &curSelectParam)
{
	for(unsigned category = 0; category < mutationParameterTrace.size(); category++)
	{
		unsigned aaRange[2];
		SequenceSummary::AAToCodonRange(aa, true, aaRange);
		for (unsigned i = aaRange[0]; i < aaRange[1]; i++)
		{
			mutationParameterTrace[category][i][sample] = curMutParam[category][i];
		}
	}
	for(unsigned category = 0; category < selectionParameterTrace.size(); category++)
	{
		unsigned aaRange[2];
		SequenceSummary::AAToCodonRange(aa, true, aaRange);
		for (unsigned i = aaRange[0]; i < aaRange[1]; i++)
		{
			selectionParameterTrace[category][i][sample] = curSelectParam[category][i];
		}
	}
}


void ROCTrace::updateExpressionTrace(unsigned sample, unsigned geneIndex, std::vector<std::vector <double>> &currentExpressionLevel)
{
	for(unsigned category = 0; category < expressionTrace.size(); category++)
	{
		expressionTrace[category][geneIndex][sample] = currentExpressionLevel[category][geneIndex];
	}
}


void ROCTrace::updateMixtureProbabilitiesTrace(unsigned samples, std::vector<double> &categoryProbabilities)
{
	for(unsigned category = 0; category < mixtureProbabilitiesTrace.size(); category++)
	{
		mixtureProbabilitiesTrace[category][samples] = categoryProbabilities[category];
	}
}

//----------------------------------------------------
//----------------------R WRAPPERS--------------------
//----------------------------------------------------


std::vector<double> ROCTrace::getExpressionAcceptanceRatioTraceByMixtureElementForGeneR(unsigned mixtureElement, unsigned geneIndex)
{
	std::vector<double> RV;
	bool checkGene = checkIndex(geneIndex, 1, expressionAcceptanceRatioTrace.size());
	bool checkMixtureElement = checkIndex(mixtureElement, 1, mixtureProbabilitiesTrace.size());
	if (checkGene && checkMixtureElement)
	{
		RV = getExpressionAcceptanceRatioTraceByMixtureElementForGene(mixtureElement - 1, geneIndex - 1);
	}
	return RV;
}


std::vector<double> ROCTrace::getExpressionTraceForGeneR(unsigned geneIndex)
{
	std::vector<double> RV;
	bool checkGene = checkIndex(geneIndex, 1, expressionTrace[0].size());
	if (checkGene)
	{
		RV = getExpressionTraceForGene(geneIndex - 1);
	}
	return RV;
}


std::vector<double> ROCTrace::getExpressionTraceByMixtureElementForGeneR(unsigned mixtureElement, unsigned geneIndex)
{
	std::vector<double> RV;
	bool checkMixtureElement = checkIndex(mixtureElement, 1, mixtureProbabilitiesTrace.size());
	bool checkGene = checkIndex(geneIndex, 1, expressionTrace[0].size());
	if (checkMixtureElement && checkGene)
	{
		RV = getExpressionTraceByMixtureElementForGene(mixtureElement - 1, geneIndex - 1);
	}
	return RV;
}


std::vector<double> ROCTrace::getMutationParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon) 
{
	std::vector<double> RV;
	bool checkMixtureElement = checkIndex(mixtureElement, 1, mixtureProbabilitiesTrace.size());
	if (checkMixtureElement)
	{
		RV = getMutationParameterTraceByMixtureElementForCodon(mixtureElement - 1, codon);  
	}
	return RV;
}


std::vector<double> ROCTrace::getSelectionParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon) 
{
	std::vector<double> RV;
	bool checkMixtureElement = checkIndex(mixtureElement, 1, mixtureProbabilitiesTrace.size());
	if (checkMixtureElement)
	{
		RV = getSelectionParameterTraceByMixtureElementForCodon(mixtureElement - 1, codon);  
	}
	return RV;
}

std::vector<unsigned> ROCTrace::getMixtureAssignmentTraceForGeneR(unsigned geneIndex)
{
	std::vector <unsigned> RV;
	bool checkGene = checkIndex(geneIndex, 1, mixtureAssignmentTrace[0].size());
	if (checkGene)
	{
		RV = getMixtureAssignmentTraceForGene(geneIndex - 1);
	}
	return RV;
}

std::vector<double> ROCTrace::getMixtureProbabilitiesTraceForMixtureR(unsigned mixtureIndex)
{
	std::vector<double> RV;
	bool check = checkIndex(mixtureIndex, 1, mixtureProbabilitiesTrace.size());
	if (check)
	{
		RV = getMixtureProbabilitiesTraceForMixture(mixtureIndex - 1);
	}
	return RV;
}


#ifndef STANDALONE

RCPP_MODULE(ROCTrace_mod)
{
	class_<ROCTrace>( "ROCTrace" )
		//These methods have only a C++ implementation
		.method("getSPhiTrace", &ROCTrace::getSPhiTrace)
		.method("getAPhiTrace", &ROCTrace::getAPhiTrace)
		.method("getSphiAcceptanceRatioTrace", &ROCTrace::getSphiAcceptanceRatioTrace)
		.method("getCspAcceptanceRatioTraceForAA", &ROCTrace::getCspAcceptanceRatioTraceForAA)


		//These methods have specific R wrappers
		.method("getExpressionAcceptanceRatioTraceByMixtureElementForGene", &ROCTrace::getExpressionAcceptanceRatioTraceByMixtureElementForGeneR)
		.method("getExpressionTraceForGene", &ROCTrace::getExpressionTraceForGeneR)
		.method("getExpressionTraceByMixtureElementForGene", &ROCTrace::getExpressionTraceByMixtureElementForGeneR)
		.method("getMutationParameterTraceByMixtureElementForCodon", &ROCTrace::getMutationParameterTraceByMixtureElementForCodonR)
		.method("getSelectionParameterTraceByMixtureElementForCodon", &ROCTrace::getSelectionParameterTraceByMixtureElementForCodonR)
		.method("getMixtureAssignmentTraceForGene", &ROCTrace::getMixtureAssignmentTraceForGeneR)
		.method("getMixtureProbabilitiesTraceForMixture", &ROCTrace::getMixtureProbabilitiesTraceForMixtureR)
		;
}
#endif
