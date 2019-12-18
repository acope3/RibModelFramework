#ifndef PANSEPARAMETER_H
#define PANSEPARAMETER_H


#include "../base/Trace.h"
#include "../base/Parameter.h"
#include "../SequenceSummary.h"


#include <vector>
#include <random>
#include <string>
#include <iostream>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

class PANSEParameter: public Parameter {
	private:
        std::vector <double> partitionFunction_proposed;
        std::vector <double> partitionFunction;

        bool fix_alpha=false;
        bool fix_lp=false;
        bool fix_nse=false;
        
        double std_partitionFunction;
        unsigned numAcceptForPartitionFunction;

		double bias_csp;

	public:
        //Testing Functions
        std::vector<double> oneMixLambda();
        std::vector<double> oneMixAlpha();
        std::vector<double> oneMixNSE();

		//Constructors & Destructors:
		explicit PANSEParameter();
		PANSEParameter(std::string filename);
		PANSEParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
				std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer = true,
				std::string _mutationSelectionState = "allUnique");
		PANSEParameter& operator=(const PANSEParameter& rhs);
		virtual ~PANSEParameter();


		//Initialization, Restart, Index Checking:
		void initPANSEParameterSet();
		void initPANSEValuesFromFile(std::string filename);
		void writeEntireRestartFile(std::string filename);
		void writePANSERestartFile(std::string filename);
		void initFromRestartFile(std::string filename);

		void initAllTraces(unsigned samples, unsigned num_genes, bool estimateSynthesisRate=true);
		void initAlpha(double alphaValue, unsigned mixtureElement, std::string codon); //R?
		void initLambdaPrime(double lambdaPrimeValue, unsigned mixtureElement, std::string codon); //R?
		void initNonsenseErrorRate(double nonsenseErrorRateValue, unsigned mixtureElement, std::string codon);
		void initMutationSelectionCategories(std::vector<std::string> files, unsigned numCategories,
				unsigned paramType); //TODO: function needs to be changed
		void fixAlpha();
		void fixLambdaPrime();
		void fixNSERate();

        //CSP Read Functions:
        void readAlphaValues(std::string filename);
        void readLambdaValues(std::string filename);
        void readNSEValues(std::string filename);

		//Trace Functions:
		void updateCodonSpecificParameterTrace(unsigned sample, std::string codon);
        void updatePartitionFunctionTrace(unsigned sample);


		//CSP Functions:
		double getCurrentCodonSpecificProposalWidth(unsigned index);
		void proposeCodonSpecificParameter();
		void updateCodonSpecificParameter(std::string grouping); //Adds to an update queue
        void completeUpdateCodonSpecificParameter();

        //partitionFunction Functions: Mostly tested, see comments.
        double getPartitionFunction(unsigned mixtureCategory, bool proposed); //TODO: test
        virtual void proposePartitionFunction(); //TODO: test
        void setPartitionFunction(double newPartitionFunction, unsigned mixtureCategory); //TODO: test
        double getCurrentPartitionFunctionProposalWidth(); //TODO: test
        unsigned getNumAcceptForPartitionFunction(); //Only for unit testing.
        void updatePartitionFunction(); //TODO: test

		//Adaptive Width Functions:
		void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth, unsigned lastIteration, bool adapt); //may make virtual
		void adaptPartitionFunctionProposalWidth(unsigned adaptationWidth, bool adapt);

		//Other functions:
		double getParameterForCategory(unsigned category, unsigned paramType, std::string codon, bool proposal);







		//R Section:

#ifndef STANDALONE

		//Constructors & Destructors:
		PANSEParameter(std::vector<double> stdDevSynthesisRate, std::vector<unsigned> geneAssignment,
			std::vector<unsigned> _matrix, bool splitSer = true);
		PANSEParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
			bool splitSer = true, std::string _mutationSelectionState = "allUnique");



		//Initialization, Restart, Index Checking:
		void initAlphaR(double alphaValue, unsigned mixtureElement, std::string codon);
		void initLambdaPrimeR(double lambdaPrimeValue, unsigned mixtureElement, std::string codon);
        void initNSERateR(double NSETRateValue, unsigned mixtureElement, std::string codon);
		void initMutationSelectionCategoriesR(std::vector<std::string> files, unsigned numCategories, std::string paramType);

		//CSP Functions:
		std::vector<std::vector<double>> getProposedAlphaParameter();
		std::vector<std::vector<double>> getProposedLambdaPrimeParameter();
		std::vector<std::vector<double>> getCurrentAlphaParameter();
		std::vector<std::vector<double>> getCurrentLambdaPrimeParameter();
        std::vector<std::vector<double>> getProposedNSERateParameter();
        std::vector<std::vector<double>> getCurrentNSERateParameter();

		void setProposedAlphaParameter(std::vector<std::vector<double>> alpha);
        void setProposedLambdaPrimeParameter(std::vector<std::vector<double>> lambdaPrime);
        void setProposedNSERateParameter(std::vector<std::vector<double>> nseRate);
		void setCurrentAlphaParameter(std::vector<std::vector<double>> alpha);
		void setCurrentLambdaPrimeParameter(std::vector<std::vector<double>> lambdaPrime);
        void setCurrentNSERateParameter(std::vector<std::vector<double>> nseRate);

		//Posterior, Variance, and Estimates Functions:
		double getAlphaPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon);
		double getLambdaPrimePosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon);

		double getAlphaVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased);
		double getLambdaPrimeVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased);
        double getNSERateVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased);



		//Other Functions:
		double getParameterForCategoryR(unsigned mixtureElement, unsigned paramType, std::string codon, bool proposal);

#endif //STANDALONE

	protected:
};

#endif // PANSEPARAMETER_H
