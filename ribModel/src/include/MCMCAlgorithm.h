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

		bool estimateSynthesisRate;
		bool estimateCodonSpecificParameter;
		bool estimateHyperParameter;
		bool writeRestartFile;
		std::vector<double> likelihoodTrace;

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

		void setEstimateSynthesisRate(bool in) {estimateSynthesisRate = in;}
		void setEstimateCodonSpecificParameter(bool in) {estimateCodonSpecificParameter = in;}
		void setEstimateHyperParameter(bool in) {estimateHyperParameter = in;}

		void setRestartFileSettings(std::string filename, unsigned interval, bool multiple);

		std::vector<double> getLogLikelihoodTrace() {return likelihoodTrace;}
		double getLogLikelihoodPosteriorMean(unsigned samples);

	protected:
};

#endif // MCMCALGORITHM_H
