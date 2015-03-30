#ifndef MCMCALGORITHM_H
#define MCMCALGORITHM_H

#include <vector>

#include "../include/Genome.h"
#include "../include/ROCModel.h"
#include "../include/ROCParameter.h"

class MCMCAlgorithm
{
    private:
        unsigned samples;
        unsigned thining;
        unsigned adaptiveWidth;

        bool estimateExpression;
        bool estimateCodonSpecificParameter;
        bool estimateHyperParameter;

        std::vector<double> likelihoodTrace;

        double acceptRejectExpressionLevelForAllGenes(Genome& genome, ROCParameter& parameter, ROCModel& model, int iteration);
        void acceptRejectCodonSpecificParameter(Genome& genome, ROCParameter& parameter, ROCModel& model, int iteration);
        void acceptRejectHyperParameter(int numGenes, ROCParameter& parameter, ROCModel& model, int iteration);

    public:
        explicit MCMCAlgorithm();
        MCMCAlgorithm(int samples, int thining, bool _estimateExpression = true, bool _estimateCodonSpecificParameter = true, bool _estimateHyperParameter = true);
        virtual ~MCMCAlgorithm();
        MCMCAlgorithm(const MCMCAlgorithm& other);

        void run(Genome& genome, ROCModel& model, ROCParameter& parameter);

        bool isEstimateExpression() {return estimateExpression;}
        bool isEstimateCodonSpecificParameter() {return estimateCodonSpecificParameter;}
        bool isEstimateHyperParameter() {return estimateHyperParameter;}

        void setEstimateExpression(bool in) {estimateExpression = in;}
        void setEstimateCodonSpecificParameter(bool in) {estimateCodonSpecificParameter = in;}
        void setEstimateHyperParameter(bool in) {estimateHyperParameter = in;}


    protected:
};

#endif // MCMCALGORITHM_H
