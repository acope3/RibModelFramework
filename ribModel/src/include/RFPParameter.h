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
#include "Parameter.h"

class RFPParameter : public Parameter
{
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
    
    const std::string groupList[61] = {"GCA", "GCC", "GCG", "GCT", "TGC", "TGT", "GAC", "GAT", "GAA", "GAG",
        "TTC", "TTT", "GGA", "GGC", "GGG", "GGT", "CAC", "CAT", "ATA", "ATC",
        "ATT", "AAA", "AAG", "CTA", "CTC", "CTG", "CTT", "TTA", "TTG", "ATG",
        "AAC", "AAT", "CCA", "CCC", "CCG", "CCT", "CAA", "CAG", "AGA", "AGG",
        "CGA", "CGC", "CGG", "CGT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC",
        "ACG", "ACT", "GTA", "GTC", "GTG", "GTT", "TGG", "TAC", "TAT", "AGC",
        "AGT"};

		std::vector<double> proposediidSum; //might go away
		std::vector<double> currentiidSum; //might go away

		// proposal bias and std for codon specific parameter -- probably can move up to Parameter
		double bias_csp;
		std::vector<double> std_csp;
		std::vector<double> prev_std_csp;		
		// functions
		std::vector<double> propose(std::vector<double> currentParam, double (*proposal)(double a, double b), double A, std::vector<double> B);
		//might change at the moment

	public:
		//Constructors & Destructors:
		explicit RFPParameter(std::string filename);
		RFPParameter(double sphi, unsigned _numMixtures,
				std::vector<unsigned> geneAssignment, std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer = true, std::string _mutationSelectionState = "allUnique");
		virtual ~RFPParameter();
		RFPParameter(const RFPParameter& other);
		RFPParameter& operator=(const RFPParameter& rhs);
#ifndef STANDALONE
		RFPParameter(double sphi, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix, bool splitSer = true);
		RFPParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, bool splitSer = true, std::string _mutationSelectionState = "allUnique");
#endif
	
		//Inititalization Functions:
		void initRFPParameterSet();
		void initAlpha(std::vector<double> alphaValues, unsigned mixtureElement, char aa);
		void initLambdaPrime(std::vector<double> lambdaPrimeValues, unsigned mixtureElement, char aa);
		void initMutationSelectionCategories(std::vector<std::string> files, unsigned numCategories, unsigned paramType); //still need to alter, not sure how
		

		//Restart File functions:	
		void writeEntireRestartFile(std::string filename); //maybe move up?
		void writeRFPRestartFile(std::string filename);
		void initFromRestartFile(std::string filename); //maybe move up?
		void initRFPValuesFromFile(std::string filename);
		
		//Trace functions:
		RFPTrace& getTraceObject() {return traces;}
		void initAllTraces(unsigned samples, unsigned num_genes, unsigned adaptiveSamples) {traces.initAllTraces(samples, num_genes, adaptiveSamples, 
				numMutationCategories, numSelectionCategories, numParam, numMixtures, categories);}
		//update trace functions
		virtual void updateSphiTrace(unsigned sample) {traces.updateSphiTrace(sample, Sphi);}
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned geneIndex){traces.updateSynthesisRateTrace(sample, geneIndex, currentSynthesisRateLevel);}
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex) {traces.updateMixtureAssignmentTrace(sample, geneIndex, mixtureAssignment[geneIndex]);}
		virtual void updateMixtureProbabilitiesTrace(unsigned samples) {traces.updateMixtureProbabilitiesTrace(samples, categoryProbabilities);}
		void updateCodonSpecificParameterTrace(unsigned sample, char aa) {traces.updateCodonSpecificParameterTrace(sample, aa, currentLambdaPrimeParameter, 
				currentAlphaParameter);}
		

		// functions to manage codon specific parameter
		void updateCodonSpecificParameter(char aa);
		double getPreviousCodonSpecificProposalWidth(unsigned aa);
		double getCurrentCodonSpecificProposalWidth(unsigned aa);
		std::vector<std::vector<double>> getCurrentAlphaParameter() {return currentAlphaParameter;}
		std::vector<std::vector<double>> getCurrentLambdaPrimeParameter() {return currentLambdaPrimeParameter;}
		void proposeCodonSpecificParameter();
		

		//Proposal Widths:
		void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth);
		virtual void adaptSphiProposalWidth(unsigned adaptationWidth);
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptationWidth);


		//Posterior Mean functions:
		virtual double getSphiPosteriorMean(unsigned samples);
		virtual double getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement);
		double getAlphaPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);
		double getLambdaPrimePosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);

		//Variance Functions:
		virtual double getSphiVariance(unsigned samples, bool unbiased = true );
		virtual double getSynthesisRateVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement, bool unbiased = true);
		double getAlphaVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);
		double getLambdaPrimeVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);
	

		// Other functions
		void getParameterForCategory(unsigned category, unsigned parameter, char aa, bool proposal, double* returnValue);
		virtual std::vector <double> getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex);
        virtual std::string getGrouping(unsigned index) {return groupList[index];}

		//Still need????
		std::vector<std::vector<double>> calculateSelectionCoefficients(unsigned sample, unsigned mixture);
#ifndef STANDALONE
		SEXP calculateSelectionCoefficientsR(unsigned sample, unsigned mixture);
#endif
		double getCurrentIidSum(unsigned aaindex) {return currentiidSum[aaindex];}
		double getProposedIidSum(unsigned aaindex) {return proposediidSum[aaindex];}





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
