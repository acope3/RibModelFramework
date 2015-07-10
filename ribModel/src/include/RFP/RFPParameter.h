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
		//members

		//Codon Specific Parameters
		std::vector<std::vector<double>> currentAlphaParameter;
		std::vector<std::vector<double>> proposedAlphaParameter;
		std::vector<std::vector<double>> currentLambdaPrimeParameter;
		std::vector<std::vector<double>> proposedLambdaPrimeParameter;

		std::vector<std::vector<double>> lambdaValues;
		std::vector<unsigned> numAcceptForAlphaAndLambdaPrime;

		// proposal bias and std for codon specific parameter -- probably can move up to Parameter
		double bias_csp;
		std::vector<double> std_csp;

		// functions
		std::vector<double> propose(std::vector<double> currentParam, double (*proposal)(double a, double b), double A,
				std::vector<double> B);
		//might change at the moment

	public:
		//Constructors & Destructors:
		RFPParameter();
		explicit RFPParameter(std::string filename);
		RFPParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
				std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer = true,
				std::string _mutationSelectionState = "allUnique");
		virtual ~RFPParameter();
		RFPParameter& operator=(const RFPParameter& rhs);
#ifndef STANDALONE
		RFPParameter(double sphi, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix, bool splitSer = true);
		RFPParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, bool splitSer = true, std::string _mutationSelectionState = "allUnique");
#endif

		//Inititalization Functions:
		void initRFPParameterSet();
		void initAlpha(double alphaValue, unsigned mixtureElement, std::string codon);
		void initLambdaPrime(double lambdaPrimeValue, unsigned mixtureElement, std::string codon);
		void initMutationSelectionCategories(std::vector<std::string> files, unsigned numCategories,
				unsigned paramType); //still need to alter, not sure how

		//Restart File functions:	
		void writeEntireRestartFile(std::string filename); //maybe move up?
		void writeRFPRestartFile(std::string filename);
		void initFromRestartFile(std::string filename); //maybe move up?
		void initRFPValuesFromFile(std::string filename);

		//Trace functions:
		RFPTrace& getTraceObject()
		{
			return traces;
		}
		void initAllTraces(unsigned samples, unsigned num_genes)
		{
			traces.initAllTraces(samples, num_genes, numMutationCategories, numSelectionCategories, numParam,
					numMixtures, categories, groupList.size());
		}
		//update trace functions
		virtual void updateSphiTrace(unsigned sample)
		{
			traces.updateSphiTrace(sample, Sphi);
		}
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned geneIndex)
		{
			traces.updateSynthesisRateTrace(sample, geneIndex, currentSynthesisRateLevel);
		}
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex)
		{
			traces.updateMixtureAssignmentTrace(sample, geneIndex, mixtureAssignment[geneIndex]);
		}
		virtual void updateMixtureProbabilitiesTrace(unsigned samples)
		{
			traces.updateMixtureProbabilitiesTrace(samples, categoryProbabilities);
		}
		void updateCodonSpecificParameterTrace(unsigned sample, std::string codon)
		{
			traces.updateCodonSpecificParameterTrace(sample, codon, currentLambdaPrimeParameter, currentAlphaParameter);
		}

		// functions to manage codon specific parameter
		void updateCodonSpecificParameter(std::string grouping);
		double getCurrentCodonSpecificProposalWidth(unsigned index);
		std::vector<std::vector<double>> getCurrentAlphaParameter()
		{
			return currentAlphaParameter;
		}
		std::vector<std::vector<double>> getCurrentLambdaPrimeParameter()
		{
			return currentLambdaPrimeParameter;
		}
		void proposeCodonSpecificParameter();

		//Proposal Widths:
		void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth); //may make virtual
		virtual void adaptSphiProposalWidth(unsigned adaptationWidth);
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptationWidth);

		//Posterior Mean functions:
		virtual double getSphiPosteriorMean(unsigned samples);
		virtual double getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement);
		double getAlphaPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);
		double getLambdaPrimePosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);

		//Variance Functions:
		virtual double getSphiVariance(unsigned samples, bool unbiased = true);
		virtual double getSynthesisRateVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement,
				bool unbiased = true);
		double getAlphaVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);
		double getLambdaPrimeVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased =
				true);

		// Other functions
		double getParameterForCategory(unsigned category, unsigned paramType, std::string codon, bool proposal);
		virtual std::vector<double> getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex);
		virtual std::string getGrouping(unsigned index)
		{
			return groupList[index];
		}

		//Statics:
		static const unsigned alp;
		static const unsigned lmPri;

		//Still need????
		std::vector<std::vector<double>> calculateSelectionCoefficients(unsigned sample, unsigned mixture);
#ifndef STANDALONE
		SEXP calculateSelectionCoefficientsR(unsigned sample, unsigned mixture);
#endif

		//More R Wrappers:
		void initAlphaR(double alphaValue, unsigned mixtureElement, std::string codon);
		void initLambdaPrimeR(double lambdaPrimeValue, unsigned mixtureElement, std::string codon);
		double getParameterForCategoryR(unsigned mixtureElement, unsigned paramType, std::string codon, bool proposal);

		//Need to alter these functions, not going to deal with these at this time-----------------------------------
		/*//R wrapper functions
		 void initMutationSelectionCategoriesR(std::vector<std::string> files, unsigned numCategories, std::string paramType);
		 double getMutationPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon)
		 {
		 double rv = -1.0;
		 bool check = checkIndex(mixtureElement, 1, numMixtures);
		 if (check)
		 {
		 rv = getMutationPosteriorMean(mixtureElement - 1, samples, codon);
		 }
		 return rv;
		 }
		 double getSelectionPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon)
		 {
		 double rv = -1.0;
		 bool check = checkIndex(mixtureElement, 1, numMixtures);
		 if (check)
		 {
		 rv = getSelectionPosteriorMean(mixtureElement - 1, samples, codon);
		 }
		 return rv;
		 }
		 double getMutationVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased)
		 {
		 double rv = -1.0;
		 bool check = checkIndex(mixtureElement, 1, numMixtures);
		 if (check)
		 {
		 rv = getMutationVariance(mixtureElement - 1, samples, codon, unbiased);
		 }
		 return rv;
		 }
		 double getSelectionVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased)
		 {
		 double rv = -1.0;
		 bool check = checkIndex(mixtureElement, 1, numMixtures);
		 if (check)
		 {
		 rv = getSelectionVariance(mixtureElement - 1, samples, codon, unbiased);
		 }
		 return rv;
		 }
		 */

	protected:

};

#endif // RFPPARAMETER_H
