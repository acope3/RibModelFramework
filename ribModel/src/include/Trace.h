#ifndef TRACE_H
#define TRACE_H

#include <iostream>
#include <vector>
#include <cctype>
#include "SequenceSummary.h"
#include "thetaK.h"
#include "Trace.h"

class Trace
{
	private:
		std::vector<double> sPhiTrace; //samples
		std::vector<double> aPhiTrace; //Not used yet
		std::vector<double> sphiAcceptanceRatioTrace; //samples
		std::vector<std::vector<std::vector<double>>> synthesisRateAcceptanceRatioTrace; //order: expressionCategory, gene, sample
		std::vector<std::vector<double>> cspAcceptanceRatioTrace; //order: codon, sample
		std::vector<std::vector<std::vector<double>>> synthesisRateTrace; //order: expressioncategoy, gene, samples
		std::vector<std::vector<unsigned>> mixtureAssignmentTrace; //order: numGenes, samples
		std::vector<std::vector<double>> mixtureProbabilitiesTrace; //order: numMixtures, samples

	public:
		//Constructor & Destructors
		Trace();
		~Trace();


		//Init functions
		void initAllTraces(unsigned samples, unsigned num_genes, unsigned numSelectionCategories, unsigned numMixtures, std::vector<thetaK> &_categories);
		void initBaseTraces(unsigned samples, unsigned num_genes, unsigned numSelectionCategories, unsigned numMixtures, std::vector<thetaK> &_categories);
		bool checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound);

		void initSphiTrace(unsigned samples);
		//aPhi will go here once implemented.
		void initSynthesisRateAcceptanceRatioTrace(unsigned samples, unsigned num_genes, unsigned numExpressionCategories);
		void initSynthesisRateTrace(unsigned samples, unsigned num_genes, unsigned numExpressionCategories); 
		void initMixtureAssignmentTrace(unsigned samples, unsigned num_genes); 
		void initMixtureProbabilitesTrace(unsigned samples, unsigned numMixtures);


		//Getter functions
		std::vector<double> getSPhiTrace() {return sPhiTrace;}
		std::vector<double> getExpectedPhiTrace();
		std::vector<double> getAPhiTrace() {return aPhiTrace;}
		std::vector<double> getSphiAcceptanceRatioTrace() {return sphiAcceptanceRatioTrace;}
		std::vector<double> getSynthesisRateAcceptanceRatioTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex); 
		std::vector<double> getCspAcceptanceRatioTraceForAA(char aa);
		std::vector<double> getSynthesisRateTraceForGene(unsigned geneIndex);//will build the trace appropriately based on what cat you are in
		std::vector<double> getSynthesisRateTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex);
		std::vector<unsigned> getMixtureAssignmentTraceForGene(unsigned geneIndex) {return mixtureAssignmentTrace[geneIndex];}
		std::vector<double> getMixtureProbabilitiesTraceForMixture(unsigned mixtureIndex) {return mixtureProbabilitiesTrace[mixtureIndex];}		

		unsigned getSynthesisRateCategory(unsigned mixtureElement) {return categories->at(mixtureElement).delEta;}
		//Update functions	
		void updateSphiTrace(unsigned sample, double Sphi) {sPhiTrace[sample] = Sphi;}
		void updateSphiAcceptanceRatioTrace(double acceptanceLevel) {sphiAcceptanceRatioTrace.push_back(acceptanceLevel);}
		void updateSynthesisRateAcceptanceRatioTrace(unsigned category, unsigned geneIndex, double acceptanceLevel);
		void updateCspAcceptanceRatioTrace(unsigned codonIndex, double acceptanceLevel) {cspAcceptanceRatioTrace[codonIndex].push_back(acceptanceLevel);}
		void updateSynthesisRateTrace(unsigned sample, unsigned geneIndex, std::vector<std::vector <double>> &currentExpressionLevel);
		void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex, unsigned value) {mixtureAssignmentTrace[geneIndex][sample] = value;}
		void updateMixtureProbabilitiesTrace(unsigned samples, std::vector<double> &categoryProbabilities);


		//R WRAPPER FUNCTIONS

		//Getter functions
		std::vector<double> getSynthesisRateAcceptanceRatioTraceByMixtureElementForGeneR(unsigned mixtureElement, unsigned geneIndex); //R WRAPPER
		std::vector<double> getSynthesisRateTraceForGeneR(unsigned geneIndex); //R WRAPPER
		std::vector<double> getSynthesisRateTraceByMixtureElementForGeneR(unsigned mixtureElement, unsigned geneIndex); //R WRAPPER
		std::vector<unsigned> getMixtureAssignmentTraceForGeneR(unsigned geneIndex); //R WRAPPER
		std::vector<double> getMixtureProbabilitiesTraceForMixtureR(unsigned mixtureIndex); //R WRAPPER
		unsigned getNumberOfMixtures() {return mixtureProbabilitiesTrace.size();}

	protected:
		std::vector<thetaK> *categories;
};

#endif //TRACE_H
