// FONSEParameter.h

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
		double phiEpsilon;
		double phiEpsilon_proposed;

		std::vector <std::vector <double> > currentMutationParameter;
		std::vector <std::vector <double> > proposedMutationParameter;

		std::vector <std::vector <double> > currentSelectionParameter;
		std::vector <std::vector <double> > proposedSelectionParameter;
		std::vector <unsigned> numAcceptForMutationAndSelection;

		double bias_csp;
		std::vector <double> std_csp;

		std::vector <double> propose(std::vector <double> currentParam, double(*proposal)(double a, double b), double A, std::vector <double> B);

	public:

		static const unsigned dM;
		static const unsigned dOmega;

		FONSEParameter();
		FONSEParameter(const FONSEParameter &other);
		explicit FONSEParameter(std::string filename);
		FONSEParameter(std::vector<double> sphi, unsigned _numMixtures, std::vector <unsigned> geneAssignment,
			std::vector <std::vector <unsigned> > thetaKmatrix, bool splitSer = true,
			std::string _mutationSelectionState = "allUnique");

		virtual ~FONSEParameter();
#ifndef STANDALONE
		FONSEParameter(std::vector<double> sphi, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix, bool splitSer = true);
		FONSEParameter(std::vector<double> sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, bool splitSer = true, std::string _mutationSelectionState = "allUnique");
		void initCovarianceMatrix(SEXP matrix, std::string aa);
#endif

		FONSEParameter& operator=(const FONSEParameter& rhs);
		FONSETrace& getTraceObject() { return traces; }
		CovarianceMatrix& getCovarianceMatrixForAA(std::string aa);

		void writeEntireRestartFile(std::string filename);
		void writeFONSERestartFile(std::string filename);
		void initFromRestartFile(std::string filename);
		void initFONSEValuesFromFile(std::string filename);
		void initFONSEParameterSet();
		void initSelection(std::vector <double> selectionValues, unsigned mixtureElement, std::string aa);
		void initMutation(std::vector <double> mutationValues, unsigned mixtureElement, std::string aa);
		std::vector <std::vector <double> > calculateSelectionCoefficients(unsigned sample, unsigned mixture);

#ifndef STANDALONE
		SEXP calculateSelectionCoefficientsR(unsigned sample, unsigned mixture);
#endif

		void initAllTraces(unsigned samples, unsigned num_genes) {
			traces.initAllTraces(samples, num_genes, numMutationCategories, numSelectionCategories, numParam, numMixtures, categories, maxGrouping);
		}
		void initMutationSelectionCategories(std::vector <std::string> files, unsigned numCategories, unsigned paramType);

		double getCurrentCodonSpecificProposalWidth(unsigned aa);
		std::vector< std::vector <double> > getCurrentMutationParameter() { return currentMutationParameter; }
		std::vector< std::vector <double> > getCurrentSelectionParameter() { return currentSelectionParameter; }

		double getPreviousCodonSpecificProposalWidth(unsigned aa);

		double getPhiEpsilon() { return phiEpsilon; }

		void updateCodonSpecificParameter(std::string grouping);

		virtual void updateSphiTrace(unsigned sample)
		{
			for(unsigned i = 0u; i < numSelectionCategories; i++)
			{
				traces.updateSphiTrace(sample, Sphi[i], i);
			}
		}
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned geneIndex){ traces.updateSynthesisRateTrace(sample, geneIndex, currentSynthesisRateLevel); }
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex) { traces.updateMixtureAssignmentTrace(sample, geneIndex, mixtureAssignment[geneIndex]); }
		virtual void updateMixtureProbabilitiesTrace(unsigned samples) { traces.updateMixtureProbabilitiesTrace(samples, categoryProbabilities); }
		void updateCodonSpecificParameterTrace(unsigned sample, std::string grouping) {
			traces.updateCodonSpecificParameterTrace(sample, grouping,
				currentMutationParameter, currentSelectionParameter);
		}

		void proposeCodonSpecificParameter();
		void proposeHyperParameters();
		void getParameterForCategory(unsigned category, unsigned paramType, std::string aa, bool proposal, double *returnSet);

		void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth);


		virtual void setNumPhiGroupings(unsigned _phiGroupings);
		virtual void adaptSphiProposalWidth(unsigned adaptationWidth);
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptationWidth);
		virtual double getSphiPosteriorMean(unsigned samples, unsigned mixture);
		virtual std::vector <double> getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex);
		virtual double getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement);
		virtual double getSphiVariance(unsigned samples, unsigned mixture, bool unbiased = true);
		virtual double getSynthesisRateVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement, bool unbiased = true);
		double getMutationPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);
		double getSelectionPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);
		double getMutationVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);
		double getSelectionVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);
		
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

        void setTraceObject(FONSETrace _trace);
};
#endif // FONSEPARAMETER_H
