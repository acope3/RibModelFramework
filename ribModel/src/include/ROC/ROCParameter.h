#ifndef ROCPARAMETER_H
#define ROCPARAMETER_H

#include <vector>
#include <random>
#include <string>
#include <iostream>
#include <array>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

#include "ROCTrace.h"
#include "../base/Parameter.h"

class ROCParameter : public Parameter
{
	private:

		ROCTrace traces;

		//members
		std::vector<CovarianceMatrix> covarianceMatrix;

		double phiEpsilon;
		double phiEpsilon_proposed;

		std::vector <double> Aphi;
		std::vector <double> Aphi_proposed;
		std::vector <double> Sepsilon;
		std::vector <double> std_Aphi;
		std::vector <double> numAcceptForAphi;

		std::vector<std::vector<double>> currentMutationParameter;
		std::vector<std::vector<double>> proposedMutationParameter;

		std::vector<std::vector<double>> currentSelectionParameter;
		std::vector<std::vector<double>> proposedSelectionParameter;
		std::vector<unsigned> numAcceptForMutationAndSelection;


		// proposal bias and std for codon specific parameter
		double bias_csp;
		std::vector<double> std_csp;


		// functions TODO: never used?
		std::vector<double> propose(std::vector<double> currentParam, double (*proposal)(double a, double b), double A, std::vector<double> B);

	public:
		//static const members
		static const unsigned dM;
		static const unsigned dEta;

		//Constructors & Destructors:
		ROCParameter();
		explicit ROCParameter(std::string filename);
		ROCParameter(double sphi, unsigned _numMixtures,
				std::vector<unsigned> geneAssignment, std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer = true, std::string _mutationSelectionState = "allUnique");
		virtual ~ROCParameter();
#ifndef STANDALONE
		ROCParameter(double sphi, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix, bool splitSer = true);
		ROCParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, bool splitSer = true, std::string _mutationSelectionState = "allUnique");
		void initCovarianceMatrix(SEXP matrix, std::string aa);
#endif
		ROCParameter& operator=(const ROCParameter& rhs);
		ROCTrace& getTraceObject() {return traces;}
		CovarianceMatrix& getCovarianceMatrixForAA(std::string aa);
		
		void writeEntireRestartFile(std::string filename);
		void writeROCRestartFile(std::string filename);
		void initFromRestartFile(std::string filename);
		void initROCValuesFromFile(std::string filename);
		void initROCParameterSet();
		void initSelection(std::vector<double> selectionValues, unsigned mixtureElement, std::string aa);
		void initMutation(std::vector<double> mutationValues, unsigned mixtureElement, std::string aa);
		std::vector<std::vector<double>> calculateSelectionCoefficients(unsigned sample, unsigned mixture);
#ifndef STANDALONE
		SEXP calculateSelectionCoefficientsR(unsigned sample, unsigned mixture);
#endif
		void initAllTraces(unsigned samples, unsigned num_genes) {traces.initAllTraces(samples, num_genes,
				numMutationCategories, numSelectionCategories, numParam, numMixtures, categories, maxGrouping, phiGroupings);}

		void initMutationCategories(std::vector<std::string> files, unsigned numCategories);
		void initSelectionCategories(std::vector<std::string> files, unsigned numCategories);
		virtual void setNumPhiGroupings(unsigned _phiGroupings);
		double getCurrentCodonSpecificProposalWidth(unsigned aa);
		std::vector<std::vector<double>> getCurrentMutationParameter() {return currentMutationParameter;}
		std::vector<std::vector<double>> getCurrentSelectionParameter() {return currentSelectionParameter;}


		double getCurrentAphiProposalWidth(unsigned index) { return std_Aphi[index]; }
		// Phi epsilon functions
		double getPhiEpsilon() { return phiEpsilon; }


		// functions to manage codon specific parameter
		void updateCodonSpecificParameter(std::string grouping);

		// functions to manage Aphi
		double getAphi(unsigned index, bool proposed = false) { return (proposed ? Aphi_proposed[index] : Aphi[index]); }
		void setAphi(unsigned index, double aPhi) { Aphi[index] = aPhi; }
		void updateAphi(unsigned index)
		{
			Aphi[index] = Aphi_proposed[index];
			numAcceptForAphi[index]++;
		}

		// functions to manage Sepsilon
		double getSepsilon(unsigned index) { return Sepsilon[index]; }
		void setSepsilon(unsigned index, double se) { Sepsilon[index] = se; }
		//update trace functions
		virtual void updateSphiTrace(unsigned sample) {traces.updateSphiTrace(sample, Sphi);}
		void updateAphiTraces(unsigned sample) { for (unsigned i = 0; i < Aphi.size(); i++) traces.updateAphiTrace(i, sample, Aphi[i]); }
		void updateSepsilonTraces(unsigned sample) { for (unsigned i = 0; i < Sepsilon.size(); i++) traces.updateSepsilonTrace(i, sample, Sepsilon[i]); }
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned geneIndex){traces.updateSynthesisRateTrace(sample, geneIndex, currentSynthesisRateLevel);}
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex) {traces.updateMixtureAssignmentTrace(sample, geneIndex, mixtureAssignment[geneIndex]);}
		virtual void updateMixtureProbabilitiesTrace(unsigned samples) {traces.updateMixtureProbabilitiesTrace(samples, categoryProbabilities);}
		void updateCodonSpecificParameterTrace(unsigned sample, std::string grouping) {traces.updateCodonSpecificParameterTrace(sample, grouping, 
			currentMutationParameter, currentSelectionParameter, codonTable);}

		// proposal functions
		void proposeCodonSpecificParameter();
		void proposeAphi();
		void getParameterForCategory(unsigned category, unsigned parameter, std::string aa, bool proposal, double *returnValue);

		void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth);


		virtual void adaptSphiProposalWidth(unsigned adaptationWidth);
		void adaptAphiProposalWidth(unsigned adaptationWidth);
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptationWidth);
		virtual double getSphiPosteriorMean(unsigned samples);
		double getAphiPosteriorMean(unsigned index, unsigned samples);
		virtual std::vector <double> getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex);
		virtual double getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement);
		virtual double getSphiVariance(unsigned samples, bool unbiased = true);
		double getAphiVariance(unsigned index, unsigned samples, bool unbiased = true);
		virtual double getSynthesisRateVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement, bool unbiased = true);
		double getMutationPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);
		double getSelectionPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);
		double getMutationVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);
		double getSelectionVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);

		//R wrapper functions
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

	protected:
};

#endif // ROCPARAMETER_H
