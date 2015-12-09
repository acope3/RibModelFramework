#ifndef RFPPARAMETER_H
#define RFPPARAMETER_H


#include "RFPTrace.h"
#include "../base/Parameter.h"


#include <vector>
#include <random>
#include <string>
#include <iostream>


#ifndef STANDALONE
#include <Rcpp.h>
#endif


class RFPParameter: public Parameter {
	private:

		RFPTrace traces;


		//Codon Specific Parameters:
		std::vector<std::vector<double>> currentAlphaParameter;
		std::vector<std::vector<double>> proposedAlphaParameter;
		std::vector<std::vector<double>> currentLambdaPrimeParameter;
		std::vector<std::vector<double>> proposedLambdaPrimeParameter;
		std::vector<std::vector<double>> lambdaValues;
		std::vector<unsigned> numAcceptForAlphaAndLambdaPrime;
		double bias_csp;
		std::vector<double> std_csp; //TODO: proposal bias and std for codon specific parameter, probably can move up to Parameter

		//Functions:
		std::vector<double> propose(std::vector<double> currentParam, double (*proposal)(double a, double b), double A,
				std::vector<double> B);
		std::vector <double> tmp;

	public:
		//Constructors & destructors:
		explicit RFPParameter();
		RFPParameter(std::string filename);
		RFPParameter(std::vector<double> sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
				std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer = true,
				std::string _mutationSelectionState = "allUnique");
		virtual ~RFPParameter();
		RFPParameter& operator=(const RFPParameter& rhs);
#ifndef STANDALONE
		SEXP calculateSelectionCoefficientsR(unsigned sample, unsigned mixture);
		RFPParameter(std::vector<double> sphi, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix,
			bool splitSer = true);
		RFPParameter(std::vector<double> sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, bool splitSer = true,
			std::string _mutationSelectionState = "allUnique");
#endif

		std::vector <double> getTmp(); 
		void updateTmp() { tmp.push_back(currentAlphaParameter[0][0]); }
		//Initialization functions:
		void initRFPParameterSet();
		void initAlpha(double alphaValue, unsigned mixtureElement, std::string codon);
		void initLambdaPrime(double lambdaPrimeValue, unsigned mixtureElement, std::string codon);
		void initMutationSelectionCategories(std::vector<std::string> files, unsigned numCategories,
				unsigned paramType); //TODO: function needs to be changed
		virtual void setNumPhiGroupings(unsigned _phiGroupings);

		//Restart file functions:
		void writeEntireRestartFile(std::string filename);
		void writeRFPRestartFile(std::string filename);
		void initFromRestartFile(std::string filename);
		void initRFPValuesFromFile(std::string filename);


		//Trace functions:
		RFPTrace& getTraceObject();
		void initAllTraces(unsigned samples, unsigned num_genes);
		virtual void updateSphiTrace(unsigned sample);
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned geneIndex);
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex);
		virtual void updateMixtureProbabilitiesTrace(unsigned samples);
		void updateCodonSpecificParameterTrace(unsigned sample, std::string codon);


		//Codon specific parameter functions:
		void updateCodonSpecificParameter(std::string grouping);
		double getCurrentCodonSpecificProposalWidth(unsigned index);
		std::vector<std::vector<double>> getCurrentAlphaParameter();
		std::vector<std::vector<double>> getCurrentLambdaPrimeParameter();
		void proposeCodonSpecificParameter();

		//Proposal widths:
		void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth); //may make virtual
		virtual void adaptSphiProposalWidth(unsigned adaptationWidth);
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptationWidth);


		//Posterior mean functions:
		virtual double getSphiPosteriorMean(unsigned samples, unsigned mixture);
		virtual double getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement);
		double getAlphaPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);
		double getLambdaPrimePosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);


		//Variance functions:
		virtual double getSphiVariance(unsigned samples, unsigned mixture, bool unbiased = true);
		virtual double getSynthesisRateVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement,
				bool unbiased = true);
		double getAlphaVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);
		double getLambdaPrimeVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased =
				true);


		//Other functions:
		double getParameterForCategory(unsigned category, unsigned paramType, std::string codon, bool proposal);
		virtual std::vector<double> getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex);
		void calculateRFPMean(Genome& genome);
		std::vector<std::vector<double>> getProposedAlphaParameter();
		std::vector<std::vector<double>> getProposedLambdaPrimeParameter();

		//Statics:
		static const unsigned alp;
		static const unsigned lmPri;


		std::vector<std::vector<double>> calculateSelectionCoefficients(unsigned sample, unsigned mixture); //TODO: need?



		//R wrappers:
		void initAlphaR(double alphaValue, unsigned mixtureElement, std::string codon);
		void initLambdaPrimeR(double lambdaPrimeValue, unsigned mixtureElement, std::string codon);
		double getParameterForCategoryR(unsigned mixtureElement, unsigned paramType, std::string codon, bool proposal);
		void initMutationSelectionCategoriesR(std::vector<std::string> files, unsigned numCategories, std::string paramType);
		double getAlphaPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon);
		double getLambdaPrimePosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon);
		double getAlphaVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased);
		double getLambdaPrimeVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased);
    
        //Functions that should only be used in R:
    #ifndef STANDALONE
        void setCurrentAlphaParameter(Rcpp::List _alpha);
        void setProposedAlphaParameter(SEXP _alpha);
        void setCurrentLambdaPrimeParameter(SEXP _lambdaPrime);
        void setProposedLambdaPrimeParameter(SEXP _lambdaPrime);
    #endif




	protected:

};

#endif // RFPPARAMETER_H
