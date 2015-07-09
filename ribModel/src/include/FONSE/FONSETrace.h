#ifndef FONSETRACE_H
#define FONSETRACE_H

#include <iostream>
#include <vector>
#include <cctype>
#include "../base/Trace.h"

class FONSETrace : public Trace
{
private:
	std::vector< std::vector< std::vector<double> > > mutationParameterTrace; //order: mutationcategoy, numparam, samples
	std::vector< std::vector< std::vector<double> > > selectionParameterTrace; //order: selectioncategoy, numparam, samples

public:
	//Constructor & Destructors
	FONSETrace();

	//Init functions
	void initAllTraces(unsigned samples, unsigned num_genes, unsigned numMutationCategories,
		unsigned numSelectionCategories, unsigned numParam, unsigned numMixtures, std::vector<mixtureDefinition> &_categories, unsigned maxGrouping);
	void initFONSETraces(unsigned samples, unsigned numMutationCategories, unsigned numSelectionCategories, unsigned numParam);
	void initMutationParameterTrace(unsigned samples, unsigned numMutationCategories, unsigned numParam);
	void initSelectionParameterTrace(unsigned samples, unsigned numSelectionCategories, unsigned numParam);


	//Getter functions
	std::vector <double> getMutationParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon);
	std::vector <double> getSelectionParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon);

	unsigned getMutationCategory(unsigned mixtureElement) { return categories->at(mixtureElement).delM; }
	unsigned getSelectionCategory(unsigned mixtureElement) { return categories->at(mixtureElement).delEta; }
	//Update functions	
	void updateCodonSpecificParameterTrace(unsigned sample, std::string aa, std::vector<std::vector<double>> &curMutParam, std::vector<std::vector<double>> &curSelectParam);

	//R WRAPPER FUNCTIONS

	//Getter functions
	std::vector <double> getMutationParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon); //R WRAPPER 
	std::vector <double> getSelectionParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon); //R WRAPPER
};

#endif //ROCTRACE_H
