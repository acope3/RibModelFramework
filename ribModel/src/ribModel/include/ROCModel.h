#ifndef ROCMODEL_H
#define ROCMODEL_H

#include "../include/Gene.h"
#include "../include/Genome.h"
#include "../include/ROCParameter.h"

class ROCModel
{
    private:


        void calculateCodonProbabilityVector(unsigned numCodons, double* mutation, double* selection, double phi, double* codonProb);
    public:
        explicit ROCModel();
        virtual ~ROCModel();
        ROCModel(const ROCModel& other);

        double calculateLogLiklihoodPerGene(Gene& gene, int geneIndex, ROCParameter& parameter, bool proposed);

        //double calculateLikelihoodPerAA(char aa, Genome& genome, ROCParameter& parameter, bool pdM, bool pdEta);
        //double calculateLikelihoodPerGene(Gene gene, ROCParameter parameter, bool proposed);
        //double calculateLogLikelihood(Genome genome, ROCParameter parameter, unsigned category, unsigned paramType, bool proposed);
        //double calculateLogLikRatioPerGene(Gene gene, ROCParameter parameter);

    protected:
};

#endif // ROCMODEL_H
