#ifndef MCMCALGORITHM_H
#define MCMCALGORITHM_H

#include <vector>
#include "ROC/ROCModel.h"
#include "RFP/RFPModel.h"
#include "FONSE/FONSEModel.h"

class MCMCAlgorithm
{
	private:
		unsigned samples;
		unsigned thining;
		unsigned adaptiveWidth;

		unsigned lastConvergenceTest;

		bool estimateSynthesisRate;
		bool estimateCodonSpecificParameter;
		bool estimateHyperParameter;
		bool estimateMixtureAssignment;
		bool writeRestartFile;
		std::vector<double> likelihoodTrace;
		std::vector<double> tmp;
		std::string file;
		unsigned fileWriteInterval;
		bool multipleFiles;
 
		double acceptRejectSynthesisRateLevelForAllGenes(Genome& genome, Model& model, int iteration);
		void acceptRejectCodonSpecificParameter(Genome& genome, Model& model, int iteration);
		void acceptRejectHyperParameter(Genome &genome, Model& model, int iteration);

	public:
		explicit MCMCAlgorithm();
		MCMCAlgorithm(unsigned samples, unsigned thining, unsigned _adaptiveWidth = 100, bool _estimateSynthesisRate = true, bool _estimateCodonSpecificParameter = true,
		bool _estimateHyperParameter = true);
		virtual ~MCMCAlgorithm();
	

		void run(Genome& genome, Model& model, unsigned numCores = 1u, unsigned divergenceIterations = 0u);
		void varyInitialConditions(Genome& genome, Model& model, unsigned divergenceIterations);
		double calculateGewekeScore(unsigned current_iteration);

		bool isEstimateSynthesisRate() {return estimateSynthesisRate;}
		bool isEstimateCodonSpecificParameter() {return estimateCodonSpecificParameter;}
		bool isEstimateHyperParameter() {return estimateHyperParameter;}
		bool isEstimateMixtureAssignment() {return estimateMixtureAssignment;}

		void setEstimateSynthesisRate(bool in) {estimateSynthesisRate = in;}
		void setEstimateCodonSpecificParameter(bool in) {estimateCodonSpecificParameter = in;}
		void setEstimateHyperParameter(bool in) {estimateHyperParameter = in;}
		void setEstimateMixtureAssignment(bool in) {estimateMixtureAssignment = in;}

		void setRestartFileSettings(std::string filename, unsigned interval, bool multiple);

		std::vector<double> getLogLikelihoodTrace() {return likelihoodTrace;}
		double getLogLikelihoodPosteriorMean(unsigned samples);

		static std::vector<double> acf(std::vector<double>& x, int nrows, int ncols, int lagmax, bool correlation, bool demean);
		static std::vector<std::vector<double>> solveToeplitzMatrix(int lr, std::vector<double> r, std::vector<double> g);

	protected:
};

#endif // MCMCALGORITHM_H
