#ifndef ROCTRACE_H
#define ROCTRACE_H

#include <iostream>
#include <vector>
#include <cctype>
#include "../base/Trace.h"
class ROCTrace : public Trace
{
	private:
		std::vector<std::vector<std::vector<double>>> mutationParameterTrace; //order: mutationcategoy, numparam, samples
		std::vector<std::vector<std::vector<double>>> selectionParameterTrace; //order: selectioncategoy, numparam, samples
		std::vector<std::vector <double> > AphiTrace;
		std::vector<std::vector <double> > AphiAcceptanceRatioTrace;
		std::vector<std::vector <double> > SepsilonTrace;

	public:
		//Constructor & Destructors
		ROCTrace();

		//Init functions
		void initAllTraces(unsigned samples, unsigned num_genes, unsigned numMutationCategories,
				unsigned numSelectionCategories, unsigned numParam, unsigned numMixtures, std::vector<mixtureDefinition> &_categories, unsigned maxGrouping, unsigned numPhiGroupings);
		void initROCTraces(unsigned samples, unsigned numMutationCategories, unsigned numSelectionCategories, unsigned numParam, unsigned numPhiGroupings);
		void initMutationParameterTrace(unsigned samples, unsigned numMutationCategories, unsigned numParam); 
		void initSelectionParameterTrace(unsigned samples, unsigned numSelectionCategories, unsigned numParam); 
		void initAphiTrace(unsigned samples, unsigned numPhiGroupings);
		void initSepsilonTrace(unsigned samples, unsigned numPhiGroupings);

		//Getter functions
		std::vector<double> getMutationParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon);
		std::vector<double> getSelectionParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon);
		std::vector<double> getAphiTrace(unsigned index) { return AphiTrace[index]; }
		std::vector<double> getAphiAcceptanceRatioTrace(unsigned index) { return AphiAcceptanceRatioTrace[index]; }
		std::vector<double> getSepsilonTrace(unsigned index) { return SepsilonTrace[index]; }
		

		unsigned getMutationCategory(unsigned mixtureElement) {return categories->at(mixtureElement).delM;}
		unsigned getSelectionCategory(unsigned mixtureElement) {return categories->at(mixtureElement).delEta;}
		//Update functions	
		void updateCodonSpecificParameterTrace(unsigned sample, std::string aa, std::vector<std::vector<double>> &curMutParam, std::vector<std::vector<double>> &curSelectParam);
		void updateAphiTrace(unsigned index, unsigned sample, double value) { AphiTrace[index][sample] = value; }
		void updateAphiAcceptanceRatioTrace(unsigned index, double value) { AphiAcceptanceRatioTrace[index].push_back(value); }
		void updateSepsilonTrace(unsigned index, unsigned sample, double value) { SepsilonTrace[index][sample] = value; }
		//R WRAPPER FUNCTIONS

		//Getter functions
		std::vector<double> getMutationParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon); //R WRAPPER
		std::vector<double> getSelectionParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon); //R WRAPPER
		std::vector<double> getAphiTraceR(unsigned index);
		std::vector<double> getSepsilonTraceR(unsigned index);
};

#endif //ROCTRACE_H
