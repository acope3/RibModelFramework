#ifndef ROCPARAMETER_H
#define ROCPARAMETER_H

#include <vector>
#include <random>
#include <string>
#include <iostream>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

#include "Genome.h"
#include "CovarianceMatrix.h"
#include "ROCTrace.h"
#include "Trace.h"
#include "Parameter.h"

class ROCParameter : public Parameter
{
	private:

		ROCTrace traces;
		//members

		double phiEpsilon;
		double phiEpsilon_proposed;

		std::vector<std::vector<double>> currentMutationParameter;
		std::vector<std::vector<double>> proposedMutationParameter;

		std::vector<std::vector<double>> currentSelectionParameter;
		std::vector<std::vector<double>> proposedSelectionParameter;
		std::vector<unsigned> numAcceptForMutationAndSelection;
		std::vector<double> proposediidSum;
		std::vector<double> currentiidSum;

		// proposal bias and std for codon specific parameter
		double bias_csp;
		std::vector<double> std_csp;
		std::vector<double> prev_std_csp;		
		// functions
		std::vector<double> propose(std::vector<double> currentParam, double (*proposal)(double a, double b), double A, std::vector<double> B);

	public:
		//static const members
		explicit ROCParameter();
		ROCParameter(double sphi, unsigned _numMixtures,
				std::vector<unsigned> geneAssignment, std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer = true, std::string _mutationSelectionState = "allUnique");
		virtual ~ROCParameter();
#ifndef STANDALONE
		ROCParameter(double sphi, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix, bool splitSer = true);
		ROCParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, bool splitSer = true, std::string _mutationSelectionState = "allUnique");
		 
#endif
		ROCParameter(const ROCParameter& other);
		ROCParameter& operator=(const ROCParameter& rhs);
		ROCTrace& getTraceObject() {return traces;}
		void initROCParameterSet();
		void initSelection(std::vector<double> selectionValues, unsigned mixtureElement, char aa);
		void initMutation(std::vector<double> mutationValues, unsigned mixtureElement, char aa);
		std::vector<std::vector<double>> calculateSelectionCoefficients(unsigned sample, unsigned mixture);
#ifndef STANDALONE
		SEXP calculateSelectionCoefficientsR(unsigned sample, unsigned mixture);
#endif
		void initAllTraces(unsigned samples, unsigned num_genes, unsigned adaptiveSamples) {traces.initAllTraces(samples, num_genes, adaptiveSamples, 
				numMutationCategories, numSelectionCategories, numParam, numMixtures, categories);}

		void initMutationSelectionCategories(std::vector<std::string> files, unsigned numCategories, unsigned paramType);

		double getCodonSpecificProposalWidth(unsigned aa);
		std::vector<std::vector<double>> getCurrentMutationParameter() {return currentMutationParameter;}
		std::vector<std::vector<double>> getCurrentSelectionParameter() {return currentSelectionParameter;}

		// Phi epsilon functions
		double getPhiEpsilon() {return phiEpsilon;}


		// functions to manage codon specific parameter
		void updateCodonSpecificParameter(char aa);



		//update trace functions
		virtual void updateSphiTrace(unsigned sample) {traces.updateSphiTrace(sample, Sphi);}
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned geneIndex){traces.updateSynthesisRateTrace(sample, geneIndex, currentSynthesisRateLevel);}
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex) {traces.updateMixtureAssignmentTrace(sample, geneIndex, mixtureAssignment[geneIndex]);}
		virtual void updateMixtureProbabilitiesTrace(unsigned samples) {traces.updateMixtureProbabilitiesTrace(samples, categoryProbabilities);}
		void updateCodonSpecificParameterTrace(unsigned sample, char aa) {traces.updateCodonSpecificParameterTrace(sample, aa, currentMutationParameter, 
				currentSelectionParameter);}

		// poposal functions
		void proposeCodonSpecificParameter();
		void getParameterForCategory(unsigned category, unsigned parameter, char aa, bool proposal, double* returnValue);

		void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth);


		virtual void adaptSphiProposalWidth(unsigned adaptationWidth);
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptationWidth);
		virtual double getSphiPosteriorMean(unsigned samples);
		virtual std::vector <double> getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex);
		virtual double getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement);
		virtual double getSphiVariance(unsigned samples, bool unbiased = true );
		virtual double getSynthesisRateVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement, bool unbiased = true);
		double getMutationPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);
		double getSelectionPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);
		double getMutationVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);
		double getSelectionVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);


		//R wrapper functions
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

	protected:

};

#endif // ROCPARAMETER_H
