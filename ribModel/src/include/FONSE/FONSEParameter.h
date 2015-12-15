#ifndef FONSEPARAMETER_H
#define FONSEPARAMETER_H

#include <vector>
#include <random>
#include <string>
#include <iostream>


#ifndef STANDALONE
#include <Rcpp.h>
#endif

#include "../base/Parameter.h"
#include "FONSETrace.h"

class FONSEParameter : public Parameter
{
	private:
		FONSETrace traces;

		std::vector<CovarianceMatrix> covarianceMatrix;

		double phiEpsilon_proposed;
		double phiEpsilon;

		std::vector <std::vector <double> > proposedMutationParameter;
		std::vector <std::vector <double> > currentMutationParameter;
		std::vector <std::vector <double> > proposedSelectionParameter;
		std::vector <std::vector <double> > currentSelectionParameter;
		std::vector <unsigned> numAcceptForMutationAndSelection;
		double bias_csp;
		std::vector <double> std_csp;


		std::vector <double> propose(std::vector <double> currentParam, double(*proposal)(double a, double b), double A, std::vector <double> B);




	public:
		static const unsigned dM;
		static const unsigned dOmega;



		//Constructors & Destructors:
		FONSEParameter();
		explicit FONSEParameter(std::string filename);
		FONSEParameter(std::vector<double> sphi, unsigned _numMixtures, std::vector <unsigned> geneAssignment,
				   std::vector <std::vector <unsigned> > thetaKmatrix, bool splitSer = true,
				   std::string _mutationSelectionState = "allUnique");
		FONSEParameter& operator=(const FONSEParameter& rhs);
		FONSEParameter(const FONSEParameter &other); //TODO: No longer needed?
		virtual ~FONSEParameter();



		//Initialization, Restart, Index Checking:
		void initFONSEParameterSet();
		void initFONSEValuesFromFile(std::string filename);
		void writeEntireRestartFile(std::string filename);
		void writeFONSERestartFile(std::string filename);
		void initFromRestartFile(std::string filename);

		void initAllTraces(unsigned samples, unsigned num_genes);
		void initMutationSelectionCategories(std::vector <std::string> files, unsigned numCategories, unsigned paramType);



		//Trace Functions:
		FONSETrace& getTraceObject();
		virtual void updateSphiTrace(unsigned sample);
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned geneIndex);
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex);
		virtual void updateMixtureProbabilitiesTrace(unsigned samples);
		void updateCodonSpecificParameterTrace(unsigned sample, std::string grouping);



		//Covariance Functions:
		CovarianceMatrix& getCovarianceMatrixForAA(std::string aa);



		//Phi Epsilon Functions:
		double getPhiEpsilon();



		//CSP Functions:
		double getPreviousCodonSpecificProposalWidth(unsigned aa); //TODO: Implement?
		double getCurrentCodonSpecificProposalWidth(unsigned aa);
		void proposeCodonSpecificParameter();
		void updateCodonSpecificParameter(std::string grouping);



		//Posterior, Variance, and Estimates Functions:
		virtual double getSphiPosteriorMean(unsigned samples, unsigned mixture);
		virtual double getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement);
		double getMutationPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);
		double getSelectionPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);

		virtual double getSphiVariance(unsigned samples, unsigned mixture, bool unbiased = true);
		virtual double getSynthesisRateVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement, bool unbiased = true);
		double getMutationVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);
		double getSelectionVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);

		virtual std::vector <double> getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex);



		//Adaptive Width Functions:
		virtual void adaptSphiProposalWidth(unsigned adaptationWidth);
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptationWidth);
		void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth);



		//Other functions:
		virtual void setNumPhiGroupings(unsigned _phiGroupings);
		void getParameterForCategory(unsigned category, unsigned paramType, std::string aa, bool proposal, double *returnSet);
		std::vector <std::vector <double> > calculateSelectionCoefficients(unsigned sample, unsigned mixture);
		void proposeHyperParameters();












		
		void initMutationSelectionCategoriesR(std::vector<std::string> files, unsigned numCategories, std::string paramType);
#ifndef STANDALONE
		FONSEParameter(std::vector<double> sphi, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix, bool splitSer = true);
		FONSEParameter(std::vector<double> sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, bool splitSer = true, std::string _mutationSelectionState = "allUnique");
		void initCovarianceMatrix(SEXP matrix, std::string aa);
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
		SEXP calculateSelectionCoefficientsR(unsigned sample, unsigned mixture);
#endif


	void initSelection(std::vector <double> selectionValues, unsigned mixtureElement, std::string aa);
	void initMutation(std::vector <double> mutationValues, unsigned mixtureElement, std::string aa);
        void setTraceObject(FONSETrace _trace);
	std::vector< std::vector <double> > getCurrentMutationParameter() { return currentMutationParameter; }
	std::vector< std::vector <double> > getCurrentSelectionParameter() { return currentSelectionParameter; }
};
#endif // FONSEPARAMETER_H
