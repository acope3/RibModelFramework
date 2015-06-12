#ifndef ROCTRACE_H
#define ROCTRACE_H

#include <iostream>
#include <vector>
#include <cctype>
#include "SequenceSummary.h"
#include "thetaK.h"
class ROCTrace
{
	private:
		std::vector<double> sPhiTrace; //samples
		std::vector<double> aPhiTrace; //Not used yet
		std::vector<double> sphiAcceptanceRatioTrace; //samples
		std::vector<std::vector<std::vector<double>>> expressionAcceptanceRatioTrace; //order: expressionCategory, gene, sample
		std::vector<std::vector<double>> cspAcceptanceRatioTrace; //order: codon, sample
		std::vector<std::vector<std::vector<double>>> expressionTrace; //order: expressioncategoy, gene, samples
		std::vector<std::vector<std::vector<double>>> mutationParameterTrace; //order: mutationcategoy, numparam, samples
		std::vector<std::vector<std::vector<double>>> selectionParameterTrace; //order: selectioncategoy, numparam, samples
		std::vector<std::vector<unsigned>> mixtureAssignmentTrace; //order: numGenes, samples
		std::vector<std::vector<double>> mixtureProbabilitiesTrace; //order: numMixtures, samples
		

		std::vector<thetaK> *categories;

	public:
		//Constructor & Destructors
		ROCTrace();

		//Init functions
		void initAllTraces(unsigned samples, unsigned num_genes, unsigned adaptiveSamples, unsigned numMutationCategories,
				unsigned numSelectionCategories, unsigned numParam, unsigned numMixtures, std::vector<thetaK> &_categories);
		bool checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound);
		
		void initSphiTrace(unsigned samples);
		//aPhi will go here once implemented.
		void initExpressionAcceptanceRatioTrace(unsigned samples, unsigned num_genes, unsigned numExpressionCategories);
		void initExpressionTrace(unsigned samples, unsigned num_genes, unsigned numExpressionCategories); 
		void initMutationParameterTrace(unsigned samples, unsigned numMutationCategories, unsigned numParam); 
		void initSelectionParameterTrace(unsigned samples, unsigned numSelectionCategories, unsigned numParam); 
		void initMixtureAssignmentTrace(unsigned samples, unsigned num_genes); 
		void initMixtureProbabilitesTrace(unsigned samples, unsigned numMixtures); 


		//Getter functions
		std::vector<double> getSPhiTrace() {return sPhiTrace;}
		std::vector<double> getExpectedPhiTrace();
		std::vector<double> getAPhiTrace() {return aPhiTrace;}
		std::vector<double> getSphiAcceptanceRatioTrace() {return sphiAcceptanceRatioTrace;}
		std::vector<double> getExpressionAcceptanceRatioTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex); 
		std::vector<double> getCspAcceptanceRatioTraceForAA(char aa);
		std::vector<double> getExpressionTraceForGene(unsigned geneIndex);//for a given geneIndex, will build the trace appropriately based on what cat you are in
		std::vector<double> getExpressionTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex);
		std::vector<double> getMutationParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon);
		std::vector<double> getSelectionParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon);
		std::vector<unsigned> getMixtureAssignmentTraceForGene(unsigned geneIndex) {return mixtureAssignmentTrace[geneIndex];}
		std::vector<double> getMixtureProbabilitiesTraceForMixture(unsigned mixtureIndex) {return mixtureProbabilitiesTrace[mixtureIndex];}

    unsigned getMutationCategory(unsigned mixtureElement) {return categories->at(mixtureElement).delM;}
    unsigned getSelectionCategory(unsigned mixtureElement) {return categories->at(mixtureElement).delEta;}
    unsigned getExpressionCategory(unsigned mixtureElement) {return categories->at(mixtureElement).delEta;}
		//Update functions	
		void updateSphiTrace(unsigned sample, double Sphi) {sPhiTrace[sample] = Sphi;}
		void updateSphiAcceptanceRatioTrace(double acceptanceLevel) {sphiAcceptanceRatioTrace.push_back(acceptanceLevel);}
		void updateExpressionAcceptanceRatioTrace(unsigned category, unsigned geneIndex, double acceptanceLevel);
		void updateCspAcceptanceRatioTrace(unsigned codonIndex, double acceptanceLevel) {cspAcceptanceRatioTrace[codonIndex].push_back(acceptanceLevel);}
		void updateExpressionTrace(unsigned sample, unsigned geneIndex, std::vector<std::vector <double>> &currentExpressionLevel);
		void updateCodonSpecificParameterTrace(unsigned sample, char aa, std::vector<std::vector<double>> &curMutParam, std::vector<std::vector<double>> &curSelectParam);

		void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex, unsigned value) {mixtureAssignmentTrace[geneIndex][sample] = value;}
		void updateMixtureProbabilitiesTrace(unsigned samples, std::vector<double> &categoryProbabilities);


		//R WRAPPER FUNCTIONS

		//Getter functions
		std::vector<double> getExpressionAcceptanceRatioTraceByMixtureElementForGeneR(unsigned mixtureElement, unsigned geneIndex); //R WRAPPER
		std::vector<double> getExpressionTraceForGeneR(unsigned geneIndex); //R WRAPPER
		std::vector<double> getExpressionTraceByMixtureElementForGeneR(unsigned mixtureElement, unsigned geneIndex); //R WRAPPER
		std::vector<double> getMutationParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon); //R WRAPPER 
		std::vector<double> getSelectionParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon); //R WRAPPER
		std::vector<unsigned> getMixtureAssignmentTraceForGeneR(unsigned geneIndex); //R WRAPPER
		std::vector<double> getMixtureProbabilitiesTraceForMixtureR(unsigned mixtureIndex); //R WRAPPER
};

#endif //ROCTRACE_H
