#ifndef RFPPARAMETER_H
#define RFPPARAMETER_H

#include <vector>
#include <random>
#include <string>
#include <iostream>


#ifndef STANDALONE
#include <Rcpp.h>
#endif

#include "RFPTrace.h"
#include "../base/Parameter.h"

class RFPParameter: public Parameter {
	private:

		RFPTrace traces;

		std::vector<std::vector<double>> proposedAlphaParameter;
		std::vector<std::vector<double>> currentAlphaParameter;
		std::vector<std::vector<double>> proposedLambdaPrimeParameter;
		std::vector<std::vector<double>> currentLambdaPrimeParameter;
		std::vector<std::vector<double>> lambdaValues; //Currently not used.
		std::vector<unsigned> numAcceptForAlphaAndLambdaPrime;
		double bias_csp;
		std::vector<double> std_csp; //TODO: proposal bias and std for codon specific parameter, probably can move up to Parameter



	public:
		static const unsigned alp;
		static const unsigned lmPri;



		//Constructors & Destructors:
		explicit RFPParameter();
		RFPParameter(std::string filename);
		RFPParameter(std::vector<double> sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
				std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer = true,
				std::string _mutationSelectionState = "allUnique");
		RFPParameter& operator=(const RFPParameter& rhs);
		virtual ~RFPParameter();



		//Initialization, Restart, Index Checking:
		void initRFPParameterSet();
		void initRFPValuesFromFile(std::string filename);
		void writeEntireRestartFile(std::string filename);
		void writeRFPRestartFile(std::string filename);
		void initFromRestartFile(std::string filename);

		void initAllTraces(unsigned samples, unsigned num_genes);
		void initAlpha(double alphaValue, unsigned mixtureElement, std::string codon); //R?
		void initLambdaPrime(double lambdaPrimeValue, unsigned mixtureElement, std::string codon); //R?
		void initMutationSelectionCategories(std::vector<std::string> files, unsigned numCategories,
				unsigned paramType); //TODO: function needs to be changed



		//Trace Functions:
		RFPTrace& getTraceObject();
		virtual void updateSphiTrace(unsigned sample);
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned geneIndex);
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex);
		virtual void updateMixtureProbabilitiesTrace(unsigned samples);
		void updateCodonSpecificParameterTrace(unsigned sample, std::string codon);



		//CSP Functions:
		double getCurrentCodonSpecificProposalWidth(unsigned index);
		void proposeCodonSpecificParameter();
		void updateCodonSpecificParameter(std::string grouping);



		//Posterior, Variance, and Estimates Functions:
		virtual double getSphiPosteriorMean(unsigned samples, unsigned mixture);
		virtual double getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement);
		double getAlphaPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);
		double getLambdaPrimePosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);

		virtual double getSphiVariance(unsigned samples, unsigned mixture, bool unbiased = true);
		virtual double getSynthesisRateVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement,
				bool unbiased = true);
		double getAlphaVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);
		double getLambdaPrimeVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased =
				true);

		virtual std::vector<double> getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex);



		//Adaptive Width Functions:
		virtual void adaptSphiProposalWidth(unsigned adaptationWidth);
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptationWidth);
		void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth); //may make virtual



		//Other functions:
		virtual void setNumObservedPhiSets(unsigned _phiGroupings);
		double getParameterForCategory(unsigned category, unsigned paramType, std::string codon, bool proposal);
		void calculateRFPMean(Genome& genome);





		//R Section:

#ifndef STANDALONE


		//Constructors & Destructors:
		RFPParameter(std::vector<double> sphi, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix,
			bool splitSer = true);
		RFPParameter(std::vector<double> sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, bool splitSer = true,
			std::string _mutationSelectionState = "allUnique");



		//Initialization, Restart, Index Checking:
		void initAlphaR(double alphaValue, unsigned mixtureElement, std::string codon);
		void initLambdaPrimeR(double lambdaPrimeValue, unsigned mixtureElement, std::string codon);
		void initMutationSelectionCategoriesR(std::vector<std::string> files, unsigned numCategories, std::string paramType);



		//Trace Functions:
		void setTraceObject(RFPTrace _trace);
		void setCategoriesForTrace();



		//CSP Functions:
		std::vector<std::vector<double>> getProposedAlphaParameter();
		std::vector<std::vector<double>> getProposedLambdaPrimeParameter();
		std::vector<std::vector<double>> getCurrentAlphaParameter();
		std::vector<std::vector<double>> getCurrentLambdaPrimeParameter();
		void setProposedAlphaParameter(std::vector<std::vector<double>> alpha);
		void setProposedLambdaPrimeParameter(std::vector<std::vector<double>> lambdaPrime);
		void setCurrentAlphaParameter(std::vector<std::vector<double>> alpha);
		void setCurrentLambdaPrimeParameter(std::vector<std::vector<double>> lambdaPrime);



		//Posterior, Variance, and Estimates Functions:
		double getAlphaPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon);
		double getLambdaPrimePosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon);

		double getAlphaVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased);
		double getLambdaPrimeVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased);



		//Other Functions:
		double getParameterForCategoryR(unsigned mixtureElement, unsigned paramType, std::string codon, bool proposal);

#endif //STANDALONE

	protected:

};

#endif // RFPPARAMETER_H
