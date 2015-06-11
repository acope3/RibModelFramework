#ifndef ROCMODEL_H
#define ROCMODEL_H

//#include "../include/Gene.h"
//#include "../include/Genome.h"
#include "../include/ROCParameter.h"

class ROCModel
{
    private:


        void obtainCodonCount(SequenceSummary& seqsum, char curAA, int codonCount[]);
        double calculateLogLikelihoodPerAAPerGene(unsigned numCodons, int codonCount[], double mutation[], double selection[], double phiValue);
    public:
        explicit ROCModel();
        virtual ~ROCModel();
        ROCModel(const ROCModel& other);

        void calculateCodonProbabilityVector(unsigned numCodons, double* mutation, double* selection, double phi, double* codonProb);
        // Likelihood ratio functions
        void calculateLogLiklihoodRatioPerGene(Gene& gene, int geneIndex, ROCParameter& parameter, unsigned k, double* logProbabilityRatio);
        void calculateLogLikelihoodRatioPerAAPerCategory(char curAA, Genome& genome, ROCParameter& parameter, double& logAcceptanceRatioForAllMixtures);
        void calculateLogLikelihoodRatioPerAAPerCategory2(char curAA, Genome& genome, ROCParameter& parameter, double& logAcceptanceRatioForAllMixtures);

        // R wrapper
        std::vector<double> CalculateProbabilitiesForCodons(std::vector<double> mutation, std::vector<double> selection, double phi);

    protected:
};

#endif // ROCMODEL_H
