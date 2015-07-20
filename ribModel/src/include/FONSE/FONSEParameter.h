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
		explicit FONSEParameter(std::string filename);
		FONSEParameter(double sphi, unsigned _numMixtures, std::vector <unsigned> geneAssignment,
			std::vector <std::vector <unsigned> > thetaKmatrix, bool splitSer = true,
			std::string _mutationSelectionState = "allUnique");

		virtual ~FONSEParameter();
#ifndef STANDALONE
		FONSEParameter(double sphi, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix, bool splitSer = true);
		FONSEParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, bool splitSer = true, std::string _mutationSelectionState = "allUnique");

#endif

		FONSEParameter(const FONSEParameter& other);
		FONSEParameter& operator=(const FONSEParameter& rhs);
		FONSETrace& getTraceObject() { return traces; }

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

		virtual void updateSphiTrace(unsigned sample) { traces.updateSphiTrace(sample, Sphi); }
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned geneIndex){ traces.updateSynthesisRateTrace(sample, geneIndex, currentSynthesisRateLevel); }
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex) { traces.updateMixtureAssignmentTrace(sample, geneIndex, mixtureAssignment[geneIndex]); }
		virtual void updateMixtureProbabilitiesTrace(unsigned samples) { traces.updateMixtureProbabilitiesTrace(samples, categoryProbabilities); }
		void updateCodonSpecificParameterTrace(unsigned sample, std::string grouping) {
			traces.updateCodonSpecificParameterTrace(sample, grouping,
				currentMutationParameter, currentSelectionParameter);
		}

		void proposeCodonSpecificParameter();
		std::vector <double> *getParameterForCategory(unsigned category, unsigned parameter, bool proposal);

		void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth);


		virtual void adaptSphiProposalWidth(unsigned adaptationWidth);
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptationWidth);
		virtual double getSphiPosteriorMean(unsigned samples);
		virtual std::vector <double> getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex);
		virtual double getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement);
		virtual double getSphiVariance(unsigned samples, bool unbiased = true);
		virtual double getSynthesisRateVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement, bool unbiased = true);
		double getMutationPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);
		double getSelectionPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);
		double getMutationVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);
		double getSelectionVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);
		
		void initMutationSelectionCategoriesR(std::vector<std::string> files, unsigned numCategories, std::string paramType);

		double getMutationPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon);
		double getSelectionPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon);
		double getMutationVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased);
		double getSelectionVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased);

};
#endif // FONSEPARAMETER_H