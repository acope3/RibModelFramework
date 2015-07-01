#ifndef RFPTRACE_H
#define RFPTRACE_H

#include <iostream>
#include <vector>
#include <cctype>
#include "../base/Trace.h"
class RFPTrace : public Trace
{
private:
	std::vector<std::vector<std::vector<double>>> alphaParameterTrace; //order: mutationcategoy, codons, samples
	std::vector<std::vector<std::vector<double>>> lambdaPrimeParameterTrace; //order: selectioncategoy, codons, samples

public:
	//Constructor & Destructors
	RFPTrace();

	//Init functions
	void initAllTraces(unsigned samples, unsigned num_genes, unsigned numMutationCategories,
		unsigned numSelectionCategories, unsigned numParam, unsigned numMixtures, std::vector<mixtureDefinition> &_categories);
	void initROCTraces(unsigned samples, unsigned numMutationCategories, unsigned numSelectionCategories, unsigned numParam);
	void initAlphaParameterTrace(unsigned samples, unsigned numMutationCategories, unsigned numParam);
	void initLambdaPrimeParameterTrace(unsigned samples, unsigned numSelectionCategories, unsigned numParam);


	//Getter functions
	std::vector<double> getAlphaParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon);
	std::vector<double> getLambdaPrimeParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon);

	unsigned getMutationCategory(unsigned mixtureElement) { return categories->at(mixtureElement).delM; }
	unsigned getSelectionCategory(unsigned mixtureElement) { return categories->at(mixtureElement).delEta; }
	//Update functions	
	void updateCodonSpecificParameterTrace(unsigned sample, char aa, std::vector<std::vector<double>> &curAlphaParam, std::vector<std::vector<double>> &curLambdaParam);

	//R WRAPPER FUNCTIONS

	//Getter functions
	std::vector<double> getAlphaParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon); //R WRAPPER 
	std::vector<double> getLambdaParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon); //R WRAPPER
};

#endif //RFPTRACE_H
