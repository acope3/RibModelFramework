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
		std::vector<double> AphiTrace;
		std::vector<double> AphiAcceptanceRatioTrace;
		std::vector<double> SepsilonTrace;

	public:
		//Constructor & Destructors
		ROCTrace();

		//Init functions
		void initAllTraces(unsigned samples, unsigned num_genes, unsigned numMutationCategories,
				unsigned numSelectionCategories, unsigned numParam, unsigned numMixtures, std::vector<mixtureDefinition> &_categories, unsigned maxGrouping);
		void initROCTraces(unsigned samples, unsigned numMutationCategories, unsigned numSelectionCategories, unsigned numParam);
		void initMutationParameterTrace(unsigned samples, unsigned numMutationCategories, unsigned numParam); 
		void initSelectionParameterTrace(unsigned samples, unsigned numSelectionCategories, unsigned numParam); 
		void initAphiTrace(unsigned samples);
		void initSepsilonTrace(unsigned samples);

		//Getter functions
		std::vector<double> getMutationParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon);
		std::vector<double> getSelectionParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon);
		std::vector<double> getAphiTrace() { return AphiTrace; }
		std::vector<double> getAphiAcceptanceRatioTrace() { return AphiAcceptanceRatioTrace; }
		std::vector<double> getSepsilonTrace() { return SepsilonTrace; }
		

		unsigned getMutationCategory(unsigned mixtureElement) {return categories->at(mixtureElement).delM;}
		unsigned getSelectionCategory(unsigned mixtureElement) {return categories->at(mixtureElement).delEta;}
		//Update functions	
		void updateCodonSpecificParameterTrace(unsigned sample, std::string aa, std::vector<std::vector<double>> &curMutParam, std::vector<std::vector<double>> &curSelectParam);
		void updateAphiTrace(unsigned sample, double value) { AphiTrace[sample] = value; }
		void updateAphiAcceptanceRatioTrace(double value) { AphiAcceptanceRatioTrace.push_back(value); }
		void updateSepsilonTrace(unsigned sample, double value) { SepsilonTrace[sample] = value; }
		//R WRAPPER FUNCTIONS

		//Getter functions
		std::vector<double> getMutationParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon); //R WRAPPER 
		std::vector<double> getSelectionParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon); //R WRAPPER
};

#endif //ROCTRACE_H
