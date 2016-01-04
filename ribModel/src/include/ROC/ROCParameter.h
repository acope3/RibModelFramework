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

#include "../base/Trace.h"
#include "../base/Parameter.h"

class ROCParameter : public Parameter
{
	private:

		Trace traces;
		std::vector<CovarianceMatrix> covarianceMatrix;

		double phiEpsilon_proposed;
		double phiEpsilon;

		std::vector <double> Sepsilon;

		std::vector <double> Aphi_proposed;
		std::vector <double> Aphi;
		std::vector <double> std_Aphi;
		std::vector <double> numAcceptForAphi;

		std::vector<std::vector<double>> proposedMutationParameter;
		std::vector<std::vector<double>> currentMutationParameter;
		std::vector<std::vector<double>> proposedSelectionParameter;
		std::vector<std::vector<double>> currentSelectionParameter;
		std::vector<unsigned> numAcceptForMutationAndSelection;


		double bias_csp;
		std::vector<double> std_csp;
		double mutation_prior_sd;


		// functions TODO: never used?
		std::vector<double> propose(std::vector<double> currentParam, double (*proposal)(double a, double b), double A, std::vector<double> B);

	public:
		static const unsigned dM;
		static const unsigned dEta;




		//Constructors & Destructors:
		ROCParameter();
		explicit ROCParameter(std::string filename);
		ROCParameter(std::vector<double> sphi, unsigned _numMixtures,
					std::vector<unsigned> geneAssignment, std::vector<std::vector<unsigned>> thetaKMatrix,
					bool splitSer = true, std::string _mutationSelectionState = "allUnique");
		ROCParameter& operator=(const ROCParameter& rhs);
		virtual ~ROCParameter();



		//Initialization, Restart, Index Checking:
		void initROCParameterSet();
		void initROCValuesFromFile(std::string filename);
		void writeEntireRestartFile(std::string filename);
		void writeROCRestartFile(std::string filename);
		void initFromRestartFile(std::string filename);

		void initAllTraces(unsigned samples, unsigned num_genes);
		void initMutationCategories(std::vector<std::string> files, unsigned numCategories);
		void initSelectionCategories(std::vector<std::string> files, unsigned numCategories);



		//Trace Functions:
		Trace& getTraceObject();
		virtual void updateSphiTrace(unsigned sample);
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned geneIndex);
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex);
		virtual void updateMixtureProbabilitiesTrace(unsigned samples);
		void updateSepsilonTraces(unsigned sample);
		void updateAphiTraces(unsigned sample);
		void updateCodonSpecificParameterTrace(unsigned sample, std::string grouping);



		//Covariance Functions:
		CovarianceMatrix& getCovarianceMatrixForAA(std::string aa);



		// Phi Epsilon Functions:
		double getPhiEpsilon();



		//Sepsilon Functions:
		double getSepsilon(unsigned index);
		void setSepsilon(unsigned index, double se);



		//Aphi Functions:
		double getAphi(unsigned index, bool proposed = false);
		double getCurrentAphiProposalWidth(unsigned index);
		void proposeAphi();
		void setAphi(unsigned index, double aPhi);
		void updateAphi(unsigned index);



		//CSP Functions:
		double getCurrentCodonSpecificProposalWidth(unsigned aa);
		void proposeCodonSpecificParameter();
		void updateCodonSpecificParameter(std::string grouping);



		//Prior Functions:
		double getMutationPriorStandardDeviation();
		void setMutationPriorStandardDeviation(double _mutation_prior_sd);



		//Posterior, Variance, and Estimates Functions:
		virtual double getSphiPosteriorMean(unsigned samples, unsigned mixture);
		virtual double getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement);
		double getAphiPosteriorMean(unsigned index, unsigned samples);
		double getMutationPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);
		double getSelectionPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);

		virtual double getSphiVariance(unsigned samples, unsigned mixture, bool unbiased = true);
		virtual double getSynthesisRateVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement, bool unbiased = true);
		double getAphiVariance(unsigned index, unsigned samples, bool unbiased = true);
		double getMutationVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);
		double getSelectionVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);

		virtual std::vector <double> getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex);



		//Adaptive Width Functions:
		virtual void adaptSphiProposalWidth(unsigned adaptationWidth);
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptationWidth);
		void adaptAphiProposalWidth(unsigned adaptationWidth);
		void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth);



		//Other Functions:
		virtual void setNumObservedPhiSets(unsigned _phiGroupings);
		void getParameterForCategory(unsigned category, unsigned parameter, std::string aa, bool proposal, double *returnValue);
		std::vector<std::vector<double>> calculateSelectionCoefficients(unsigned sample, unsigned mixture);




		//R Section:

#ifndef STANDALONE

		//Constructors & Destructors:
		ROCParameter(std::vector<double> sphi, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix,
					bool splitSer = true);
		ROCParameter(std::vector<double> sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
					bool splitSer = true, std::string _mutationSelectionState = "allUnique");



		//Initialization, Restart, Index Checking:
		void initCovarianceMatrix(SEXP matrix, std::string aa);
		void initMutation(std::vector<double> mutationValues, unsigned mixtureElement, std::string aa);
		void initSelection(std::vector<double> selectionValues, unsigned mixtureElement, std::string aa);



 		//Trace Functions:
 		void setTraceObject(ROCTrace _trace);
 		void setCategoriesForTrace();



		//CSP Functions:
		std::vector<std::vector<double>> getProposedMutationParameter();
		std::vector<std::vector<double>> getCurrentMutationParameter();
		std::vector<std::vector<double>> getProposedSelectionParameter();
		std::vector<std::vector<double>> getCurrentSelectionParameter();


		void setProposedMutationParameter(std::vector<std::vector<double>> _proposedMutationParameter);
		void setCurrentMutationParameter(std::vector<std::vector<double>> _currentMutationParameter);
		void setProposedSelectionParameter(std::vector<std::vector<double>> _proposedSelectionParameter);
		void setCurrentSelectionParameter(std::vector<std::vector<double>> _currentSelectionParameter);



		//Posterior, Variance, and Estimates Functions:
		double getMutationPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon);
		double getSelectionPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon);
		double getMutationVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased);
		double getSelectionVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased);



		//Other Functions:
		SEXP calculateSelectionCoefficientsR(unsigned sample, unsigned mixture);
#endif //STADNALONE

	protected:
};

#endif // ROCPARAMETER_H
