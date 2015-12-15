#ifndef TRACE_H
#define TRACE_H

#include <iostream>
#include <vector>
#include <cctype>
#include "../SequenceSummary.h"
#include "../mixtureDefinition.h"

class Trace {
	private:
		std::vector<std::vector<double>> sPhiTrace; //samples
		std::vector<double> sphiAcceptanceRatioTrace; //samples
		std::vector<std::vector<std::vector<double>>>synthesisRateAcceptanceRatioTrace; //order: expressionCategory, gene, sample
		std::vector<std::vector<double>> cspAcceptanceRatioTrace;//order: codon, sample
		std::vector<std::vector<std::vector<double>>> synthesisRateTrace;//order: expressioncategoy, gene, samples
		std::vector<std::vector<unsigned>> mixtureAssignmentTrace;//order: numGenes, samples
		std::vector<std::vector<double>> mixtureProbabilitiesTrace;//order: numMixtures, samples

		public:
		//Constructor & Destructors
		Trace();
		~Trace();

		//Init functions
		void initAllTraces(unsigned samples, unsigned num_genes, unsigned numSelectionCategories, unsigned numMixtures, 
				std::vector<mixtureDefinition> &_categories, unsigned maxGrouping);
		void initBaseTraces(unsigned samples, unsigned num_genes, unsigned numSelectionCategories, unsigned numMixtures, 
				std::vector<mixtureDefinition> &_categories, unsigned maxGrouping);
		bool checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound);

		void initSphiTrace(unsigned numSelectionCategories, unsigned samples);
		//aPhi will go here once implemented.
		void initSynthesisRateAcceptanceRatioTrace(unsigned num_genes, unsigned numExpressionCategories);
		void initSynthesisRateTrace(unsigned samples, unsigned num_genes, unsigned numExpressionCategories);
		void initMixtureAssignmentTrace(unsigned samples, unsigned num_genes);
		void initMixtureProbabilitesTrace(unsigned samples, unsigned numMixtures);

		//Getter functions
		std::vector<double> getSphiTrace(unsigned selectionCategory) { return sPhiTrace[selectionCategory]; }
		std::vector<double> getExpectedPhiTrace();
		std::vector<double> getSphiAcceptanceRatioTrace()
		{	return sphiAcceptanceRatioTrace;}
		std::vector<std::vector<std::vector<double>>> getSynthesisRateTrace();
		std::vector<double> getSynthesisRateAcceptanceRatioTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex);
		std::vector<std::vector<std::vector<double>>> getSynthesisRateAcceptanceRatioTrace();
		std::vector<double> getCspAcceptanceRatioTraceForAA(std::string aa);
		std::vector<double> getSynthesisRateTraceForGene(unsigned geneIndex); //will build the trace appropriately based on what cat you are in
		std::vector<double> getSynthesisRateTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex);
		std::vector<unsigned> getMixtureAssignmentTraceForGene(unsigned geneIndex)
		{	return mixtureAssignmentTrace[geneIndex];}
		std::vector<double> getMixtureProbabilitiesTraceForMixture(unsigned mixtureIndex)
		{	return mixtureProbabilitiesTrace[mixtureIndex];}

		std::vector<std::vector<unsigned>> getMixtureAssignmentTrace();
		std::vector<std::vector<double>> getMixtureProbabilitiesTrace();
		std::vector<std::vector<double>> getCspAcceptanceRatioTrace();

		unsigned getSynthesisRateCategory(unsigned mixtureElement)
		{	return categories->at(mixtureElement).delEta;}
		//Update functions	
		void updateSphiTrace(unsigned sample, double Sphi, unsigned selectionCategory)
		{	sPhiTrace[selectionCategory][sample] = Sphi;}
		void updateSphiAcceptanceRatioTrace(double acceptanceLevel)
		{	sphiAcceptanceRatioTrace.push_back(acceptanceLevel);}
		void updateSynthesisRateAcceptanceRatioTrace(unsigned category, unsigned geneIndex, double acceptanceLevel);
		void updateCspAcceptanceRatioTrace(unsigned codonIndex, double acceptanceLevel)
		{	cspAcceptanceRatioTrace[codonIndex].push_back(acceptanceLevel);}
		void updateSynthesisRateTrace(unsigned sample, unsigned geneIndex, std::vector<std::vector <double>> &currentExpressionLevel);
		void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex, unsigned value)
		{	mixtureAssignmentTrace[geneIndex][sample] = value;}
		void updateMixtureProbabilitiesTrace(unsigned samples, std::vector<double> &categoryProbabilities);

		//R WRAPPER FUNCTIONS

		//Getter functions
		std::vector<double> getSynthesisRateAcceptanceRatioTraceByMixtureElementForGeneR(unsigned mixtureElement, unsigned geneIndex);//R WRAPPER
		std::vector<double> getSynthesisRateTraceForGeneR(unsigned geneIndex);//R WRAPPER
		std::vector<double> getSynthesisRateTraceByMixtureElementForGeneR(unsigned mixtureElement, unsigned geneIndex);//R WRAPPER
		std::vector<unsigned> getMixtureAssignmentTraceForGeneR(unsigned geneIndex);//R WRAPPER
		std::vector<double> getMixtureProbabilitiesTraceForMixtureR(unsigned mixtureIndex);//R WRAPPER
		std::vector<std::vector<double>> getSphiTraces()
		{
			return sPhiTrace;
		}
		unsigned getNumberOfMixtures()
		{	return mixtureProbabilitiesTrace.size();}
    
    
        //Only use R:
        void setSphiTraces(std::vector<std::vector<double>> _sPhiTrace);
        void setSphiAcceptanceRatioTrace(std::vector<double> _sphiAcceptanceRatioTrace);
        void setSynthesisRateTrace(std::vector<std::vector<std::vector<double>>> _synthesisRateTrace);
        void setSynthesisRateAcceptanceRatioTrace(std::vector<std::vector<std::vector<double>>>_synthesisRateAcceptanceRatioTrace);
        void setMixtureAssignmentTrace(std::vector<std::vector<unsigned>> _mixtureAssignmentTrace);
        void setMixtureProbabilitiesTrace(std::vector<std::vector<double>> _mixtureProbabilitiesTrace);
        void setCspAcceptanceRatioTrace(std::vector<std::vector<double>> _cspAcceptanceRatioTrace);
        void setCategories(std::vector<mixtureDefinition> &_categories);
    

		protected:
		std::vector<mixtureDefinition> *categories;
	};

#endif //TRACE_H
