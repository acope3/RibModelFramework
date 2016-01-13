

#ifndef TRACE_H
#define TRACE_H

#include <iostream>
#include <vector>
#include <cctype>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

#include "../mixtureDefinition.h"

class Trace {
	private:
    
		unsigned numCodonSpecificParamTypes;

		std::vector<std::vector<double>> stdDevSynthesisRateTrace; //mixture, samples
		std::vector<double> stdDevSynthesisRateAcceptanceRatioTrace; //samples TODO: Correctly sized for the time being,
        //however, it will need to be changed at some point when there are some adjustments to hyper parameter acceptance/rejection
		std::vector<std::vector<std::vector<double>>>synthesisRateAcceptanceRatioTrace; //order: expressionCategory, gene, sample
		std::vector<std::vector<double>> codonSpecificAcceptanceRatioTrace;//order: codon, sample
		std::vector<std::vector<std::vector<double>>> synthesisRateTrace;//order: expressioncategoy, gene, samples
		std::vector<std::vector<unsigned>> mixtureAssignmentTrace;//order: numGenes, samples
		std::vector<std::vector<double>> mixtureProbabilitiesTrace;//order: numMixtures, samples
		std::vector<std::vector<std::vector<std::vector<double>>>> codonSpecificParameterTrace; //order: paramType, category, numparam, samples
		//std::vector<std::vector<std::vector<double>>> codonSpecificParameterTraceTwo; //order: category, numparam, samples
		std::vector<mixtureDefinition> *categories;



		//ROC Trace:
		std::vector<std::vector <double>> synthesisOffsetTrace;
		std::vector<std::vector <double>> synthesisOffsetAcceptanceRatioTrace;
		std::vector<std::vector <double>> observedSynthesisNoiseTrace;

		//FONSE Trace:



		//RFP Trace:


		//Initialization Functions:
		void initializeSharedTraces(unsigned samples, unsigned num_genes, unsigned numSelectionCategories, unsigned numMixtures,
			std::vector<mixtureDefinition> &_categories, unsigned maxGrouping);

		void initStdDevSynthesisRateTrace(unsigned numSelectionCategories, unsigned samples);
		void initSynthesisRateAcceptanceRatioTrace(unsigned num_genes, unsigned numExpressionCategories);
		void initSynthesisRateTrace(unsigned samples, unsigned num_genes, unsigned numExpressionCategories);
		void initMixtureAssignmentTrace(unsigned samples, unsigned num_genes);
		void initMixtureProbabilitesTrace(unsigned samples, unsigned numMixtures);
		void initCodonSpecificParameterTrace(unsigned samples, unsigned numMutationCategories, unsigned numParam, unsigned paramType);


		//ROC Specific:
		void initSynthesisOffsetTrace(unsigned samples, unsigned numPhiGroupings);
		void initObservedSynthesisNoiseTrace(unsigned samples, unsigned numPhiGroupings);


		//FONSE Specific:


		//RFP Specific:

public:
	//Constructors & Destructors:
	Trace();
	virtual ~Trace();
	Trace(unsigned _numCodonSpecificParamTypes);


	//Initialization Functions:
	void initializeRFPTrace(unsigned samples, unsigned num_genes, unsigned numAlphaCategories,
		unsigned numLambdaPrimeCategories, unsigned numParam, unsigned numMixtures,
		std::vector<mixtureDefinition> &_categories, unsigned maxGrouping);
        
        
	void initializeROCTrace(unsigned samples, unsigned num_genes, unsigned numMutationCategories,
		unsigned numSelectionCategories, unsigned numParam, unsigned numMixtures, std::vector<mixtureDefinition> &_categories,
		unsigned maxGrouping, unsigned numObservedPhiSets);
        
        
	void initializeFONSETrace(unsigned samples, unsigned num_genes, unsigned numMutationCategories,
		unsigned numSelectionCategories, unsigned numParam, unsigned numMixtures,
		std::vector<mixtureDefinition> &_categories, unsigned maxGrouping);
        
        
	void initializePANSETrace(unsigned samples, unsigned num_genes, unsigned numAlphaCategories,
		unsigned numLambdaPrimeCategories, unsigned numParam, unsigned numMixtures,
		std::vector<mixtureDefinition> &_categories, unsigned maxGrouping);


        //Getter Functions:
        std::vector<double> getStdDevSynthesisRateTrace(unsigned selectionCategory);
        std::vector<double> getExpectedSynthesisRateTrace();
        std::vector<double> getStdDevSynthesisRateAcceptanceRatioTrace();
        std::vector<std::vector<std::vector<double>>> getSynthesisRateTrace();
        std::vector<double> getSynthesisRateAcceptanceRatioTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex);
        std::vector<std::vector<std::vector<double>>> getSynthesisRateAcceptanceRatioTrace();
        std::vector<double> getCodonSpecficAcceptanceRatioTraceForAA(std::string aa);
        std::vector<double> getSynthesisRateTraceForGene(unsigned geneIndex); //will build the trace appropriately based on what cat you are in
        std::vector<double> getSynthesisRateTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex);
        std::vector<unsigned> getMixtureAssignmentTraceForGene(unsigned geneIndex);
        std::vector<double> getMixtureProbabilitiesTraceForMixture(unsigned mixtureIndex);
        std::vector<std::vector<unsigned>> getMixtureAssignmentTrace();
        std::vector<std::vector<double>> getMixtureProbabilitiesTrace();
        std::vector<std::vector<double>> getCodonSpecificAcceptanceRatioTrace();
        unsigned getSynthesisRateCategory(unsigned mixtureElement);
        unsigned getCodonSpecificCategory(unsigned mixtureElement, unsigned paramType);


        //ROC Specific:
        std::vector<double> getCodonSpecificParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon, unsigned paramType);
        std::vector<double> getSynthesisOffsetTrace(unsigned index);
        std::vector<double> getSynthesisOffsetAcceptanceRatioTraceForIndex(unsigned index);
        std::vector<double> getObservedSynthesisNoiseTrace(unsigned index);
        std::vector<std::vector<std::vector<double>>> getCodonSpecificParameterTrace(unsigned paramType);
        std::vector<std::vector<double>> getSynthesisOffsetAcceptanceRatioTrace();

        //FONSE Specific:

        //RFP Specific:



        //Update Functions:
        void updateStdDevSynthesisRateTrace(unsigned sample, double stdDevSynthesisRate, unsigned synthesisRateCategory);
        void updateStdDevSynthesisRateAcceptanceRatioTrace(double acceptanceLevel);
        void updateSynthesisRateAcceptanceRatioTrace(unsigned category, unsigned geneIndex, double acceptanceLevel);
        void updateCodonSpecificAcceptanceRatioTrace(unsigned codonIndex, double acceptanceLevel);
        void updateSynthesisRateTrace(unsigned sample, unsigned geneIndex, std::vector<std::vector <double>> &currentExpressionLevel);
        void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex, unsigned value);
        void updateMixtureProbabilitiesTrace(unsigned samples, std::vector<double> &categoryProbabilities);


        //ROC Specific:
        void updateCodonSpecificParameterTraceForAA(unsigned sample, std::string aa, std::vector<std::vector<double>> &curParam, unsigned paramType);
        void updateSynthesisOffsetTrace(unsigned index, unsigned sample, double value);
        void updateSynthesisOffsetAcceptanceRatioTrace(unsigned index, double value);
        void updateObservedSynthesisNoiseTrace(unsigned index, unsigned sample, double value);

        //FONSE Specific:

        //RFP Specific:
        void updateCodonSpecificParameterTraceForCodon(unsigned sample, std::string codon, std::vector<std::vector<double>> &curParam, unsigned paramType);



        //R Section:
#ifndef STANDALONE
        //Getter Functions:
        std::vector<double> getSynthesisRateAcceptanceRatioTraceByMixtureElementForGeneR(unsigned mixtureElement, unsigned geneIndex);//R WRAPPER
        std::vector<double> getSynthesisRateTraceForGeneR(unsigned geneIndex);//R WRAPPER
        std::vector<double> getSynthesisRateTraceByMixtureElementForGeneR(unsigned mixtureElement, unsigned geneIndex);//R WRAPPER
        std::vector<unsigned> getMixtureAssignmentTraceForGeneR(unsigned geneIndex);//R WRAPPER
        std::vector<double> getMixtureProbabilitiesTraceForMixtureR(unsigned mixtureIndex);//R WRAPPER
        std::vector<std::vector<double>> getStdDevSynthesisRateTraces();
        unsigned getNumberOfMixtures();



        //Setter Functions:
        void setStdDevSynthesisRateTraces(std::vector<std::vector<double>> _stdDevSynthesisRateTrace);
        void setStdDevSynthesisRateAcceptanceRatioTrace(std::vector<double> _stdDevSynthesisRateAcceptanceRatioTrace);
        void setSynthesisRateTrace(std::vector<std::vector<std::vector<double>>> _synthesisRateTrace);
        void setSynthesisRateAcceptanceRatioTrace(std::vector<std::vector<std::vector<double>>>_synthesisRateAcceptanceRatioTrace);
        void setMixtureAssignmentTrace(std::vector<std::vector<unsigned>> _mixtureAssignmentTrace);
        void setMixtureProbabilitiesTrace(std::vector<std::vector<double>> _mixtureProbabilitiesTrace);
        void setCodonSpecificAcceptanceRatioTrace(std::vector<std::vector<double>> _cspAcceptanceRatioTrace);
        void setCategories(std::vector<mixtureDefinition> &_categories);


        //ROC Specific:
		std::vector<double> getCodonSpecificParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon, unsigned paramType);
        std::vector<std::vector<double>> getSynthesisOffsetTraceR();
        std::vector<std::vector<double>> getObservedSynthesisNoiseTraceR();


        void setSynthesisOffsetTrace(std::vector<std::vector <double> > _NoiseOffsetTrace);
        void setSynthesisOffsetAcceptanceRatioTrace(std::vector<std::vector <double> > _NoiseOffsetAcceptanceRatioTrace);
        void setObservedSynthesisNoiseTrace(std::vector<std::vector <double> > _ObservedSynthesisNoiseTrace);
        void setCodonSpecificParameterTrace(std::vector<std::vector<std::vector<double>>> _parameterTrace, unsigned paramType);


        bool checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound);

#endif
};
#endif // TRACE_H
