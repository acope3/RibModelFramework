#ifndef ROCMODEL_H
#define ROCMODEL_H

#include "../include/Gene.h"
#include "../include/Genome.h"
#include "../include/ROCParameter.h"

class ROCModel
{
    private:


        void calculateCodonProbabilityVector(unsigned numCodons, double* mutation, double* selection, double phi, double* codonProb);
        void obtainCodonCount(SequenceSummary& seqsum, char curAA, int codonCount[]);
        double calculateLogLikelihoodPerAAPerGene(unsigned numCodons, int codonCount[], SequenceSummary& seqsum, double mutation[], double selection[], double phiValue);

    public:
        explicit ROCModel();
        virtual ~ROCModel();
        ROCModel(const ROCModel& other);

        // Likelihood functions
        double calculateLogLiklihoodPerGene(Gene& gene, int geneIndex, ROCParameter& parameter, bool proposed);

        // Likelihood ratio functions
        double* calculateLogLiklihoodRatioPerGene(Gene& gene, int geneIndex, ROCParameter& parameter);
        void calculateLogLikelihoodRatioPerAAPerCategory(char curAA, Genome& genome, ROCParameter& parameter, double* logAcceptanceRatioPerCategory);
//        void calculateLogLikelihoodRatioPerCategory(Genome& genome, ROCParameter& parameter, double* logAcceptanceRatioPerCategory);

        //double calculateLikelihoodPerAA(char aa, Genome& genome, ROCParameter& parameter, bool pdM, bool pdEta);
        //double calculateLikelihoodPerGene(Gene gene, ROCParameter parameter, bool proposed);
        //double calculateLogLikelihood(Genome genome, ROCParameter parameter, unsigned category, unsigned paramType, bool proposed);
        //double calculateLogLikRatioPerGene(Gene gene, ROCParameter parameter);

    protected:
};

#endif // ROCMODEL_H
