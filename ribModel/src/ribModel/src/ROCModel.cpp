#include "../include/ROCModel.h"

#include <vector>
#include <math.h>
#include <iostream>
ROCModel::ROCModel()
{
    //ctor
}

ROCModel::~ROCModel()
{
    //dtor
}

ROCModel::ROCModel(const ROCModel& other)
{
    //copy ctor
}

double ROCModel::calculateLogLiklihoodPerGene(Gene& gene, int geneIndex, ROCParameter& parameter, bool proposed)
{
    double logLikelihood = 0.0;
    double phiValue = parameter.getExpression(geneIndex, proposed);

    SequenceSummary seqsum = gene.getSequenceSummary();
    std::string id = gene.getId();
    for(int i = 0; i < 22; i++)
    {
        char curAA = seqsum.AminoAcidArray[i];
        // skip amino acids with only one codon or stop codons
        if(curAA == 'X' || curAA == 'M' || curAA == 'W') continue;
        // skip amino acids which do not occur in current gene. Avoid useless calculations and multiplying by 0
        if(seqsum.getAAcountForAA(i) == 0) continue;

        // get codon count (total count not parameter count)
        int numCodons = seqsum.GetNumCodonsForAA(curAA);
        // get mutation and selection parameter for gene
        double mutation[numCodons - 1];
        parameter.getParameterForCategory(gene.getMutationCategory(), ROCParameter::dM, curAA, false, mutation);
        double selection[numCodons - 1];
        parameter.getParameterForCategory(gene.getDeltaEtaCategory(), ROCParameter::dEta, curAA, false, selection);

        int codonCount[numCodons];
        // prepare array for codon counts for AA
        unsigned* codonRange = SequenceSummary::AAToCodonRange(curAA);
        // get codon counts for AA
        unsigned j = 0u;
        for(unsigned i = codonRange[0]; i < codonRange[1]; i++, j++)
        {
            codonCount[j] = seqsum.getCodonCountForCodon(i);
        }

        logLikelihood += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, seqsum, mutation, selection, phiValue);
    }
    double sPhi = parameter.getSphi(false);
    return logLikelihood + std::log(ROCParameter::densityLogNorm(phiValue, (-(sPhi * sPhi) / 2), sPhi)); // + std::log(ROCParameter::densityLogNorm(phiObserved, std::log(phiValue), parameter.getPhiEpsilon() ) )
}

// for proposinng phi values
double* ROCModel::calculateLogLiklihoodRatioPerGene(Gene& gene, int geneIndex, ROCParameter& parameter)
{
    double logLikelihood = 0.0;
    double logLikelihood_proposed = 0.0;
    double phiValue = parameter.getExpression(geneIndex, false);
    double phiValue_proposed = parameter.getExpression(geneIndex, true);

    SequenceSummary seqsum = gene.getSequenceSummary();
    for(int i = 0; i < 22; i++)
    {
        char curAA = seqsum.AminoAcidArray[i];
        // skip amino acids with only one codon or stop codons
        if(curAA == 'X' || curAA == 'M' || curAA == 'W') continue;
        // skip amino acids which do not occur in current gene. Avoid useless calculations and multiplying by 0
        if(seqsum.getAAcountForAA(i) == 0) continue;

        // get codon count (total count not parameter count)
        int numCodons = seqsum.GetNumCodonsForAA(curAA);
        // get mutation and selection parameter for gene
        double mutation[numCodons - 1];
        parameter.getParameterForCategory(gene.getMutationCategory(), ROCParameter::dM, curAA, false, mutation);
        double selection[numCodons - 1];
        parameter.getParameterForCategory(gene.getDeltaEtaCategory(), ROCParameter::dEta, curAA, false, selection);

        // prepare array for codon counts for AA
        int codonCount[numCodons];
        obtainCodonCount(seqsum, curAA, codonCount);

        logLikelihood += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, seqsum, mutation, selection, phiValue);
        logLikelihood_proposed += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, seqsum, mutation, selection, phiValue_proposed);
    }
    double sPhi = parameter.getSphi(false);
    double logPhiProbability = std::log(ROCParameter::densityLogNorm(phiValue, (-(sPhi * sPhi) / 2), sPhi));
    double logPhiProbability_proposed = std::log(ROCParameter::densityLogNorm(phiValue_proposed, (-(sPhi * sPhi) / 2), sPhi));

    double currentLogLikelihood = (logLikelihood + logPhiProbability);
    double proposedLogLikelihood = (logLikelihood_proposed + logPhiProbability_proposed);

    double returnArray[3];
    returnArray[0] = proposedLogLikelihood - currentLogLikelihood - (std::log(phiValue) - std::log(phiValue_proposed));
    returnArray[1] = currentLogLikelihood;
    returnArray[2] = proposedLogLikelihood;
    return returnArray;
    //return logLikelihood + a;
}

void ROCModel::calculateCodonProbabilityVector(unsigned numCodons, double mutation[], double selection[], double phi, double codonProb[])
{
    // calculate numerator and denominator for codon probabilities
    double denominator = 1.0;
    for(unsigned i = 0; i < (numCodons - 1); i++)
    {
        codonProb[i] = std::exp( -mutation[i] - (selection[i] * phi) );
        denominator += codonProb[i];
    }
    // alphabetically last codon is reference codon!
    codonProb[numCodons - 1] = 1.0;
    // normalize codon probabilities
    for(unsigned i = 0; i < numCodons; i++) { codonProb[i] = codonProb[i] / denominator;}
}

double ROCModel::calculateLogLikelihoodPerAAPerGene(unsigned numCodons, int codonCount[], SequenceSummary& seqsum, double mutation[], double selection[], double phiValue)
{
    double logLikelihood = 0.0;
    // calculate codon probabilities
    double codonProbabilities[numCodons];
    calculateCodonProbabilityVector(numCodons, mutation, selection, phiValue, codonProbabilities);

    // calculate likelihood for current AA for this combination of selection and mutation category
    for(unsigned i = 0; i < numCodons; i++)
    {
        logLikelihood += std::log(codonProbabilities[i]) * codonCount[i];
    }
    return logLikelihood;
}

void ROCModel::calculateLogLikelihoodRatioPerAAPerCategory(char curAA, Genome& genome, ROCParameter& parameter, double* logAcceptanceRatioPerCategory)
{
    //std::cout << "=======================================" << std::endl;
    //std::cout << "mutation: " << logAcceptanceRatioPerCategory[0] << "\t selection: " << logAcceptanceRatioPerCategory[1] << std::endl;
    unsigned numMutationCat = parameter.getNumMutationCategories();
    unsigned numSelectionCat = parameter.getNumSelectionCategories();

    double logLikelihoodMutationCurrent[numMutationCat];
    double logLikelihoodMutationProposed[numMutationCat];
    double logLikelihoodSelectionCurrent[numSelectionCat];
    double logLikelihoodSelectionProposed[numSelectionCat];
    int numGenes = genome.getGenomeSize();

    for(int i = 0; i < numGenes; i++)
    {
        Gene gene = genome.getGene(i);
        SequenceSummary seqsum = gene.getSequenceSummary();
        if(seqsum.getAAcountForAA(curAA) == 0) continue;

        int numCodons = seqsum.GetNumCodonsForAA(curAA);

        int codonCount[numCodons];
        obtainCodonCount(seqsum, curAA, codonCount);

        unsigned mutCat = gene.getMutationCategory();
        unsigned selCat = gene.getDeltaEtaCategory();

        // get mutation and selection parameter for gene
        double mutation[numCodons - 1];
        double selection[numCodons - 1];
        double phiValue = parameter.getExpression(i, false);
        // calculate current likelihood
        parameter.getParameterForCategory(mutCat, ROCParameter::dM, curAA, false, mutation);
        parameter.getParameterForCategory(selCat, ROCParameter::dEta, curAA, false, selection);
        double logLikelihood = calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, seqsum, mutation, selection, phiValue);
        logLikelihoodMutationCurrent[mutCat] += logLikelihood;
        logLikelihoodSelectionCurrent[selCat] += logLikelihood;

        // calculate proposed mutation likelihood
        parameter.getParameterForCategory(mutCat, ROCParameter::dM, curAA, true, mutation);
        parameter.getParameterForCategory(selCat, ROCParameter::dEta, curAA, false, selection);
        logLikelihood = calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, seqsum, mutation, selection, phiValue);
        logLikelihoodMutationProposed[mutCat] += logLikelihood;

        // calculate proposed selection likelihood
        parameter.getParameterForCategory(mutCat, ROCParameter::dM, curAA, false, mutation);
        parameter.getParameterForCategory(selCat, ROCParameter::dEta, curAA, true, selection);
        logLikelihood = calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, seqsum, mutation, selection, phiValue);
        logLikelihoodSelectionProposed[selCat] += logLikelihood;
    }


    for(unsigned i = 0; i < numMutationCat; i++)
    {
        logAcceptanceRatioPerCategory[i] = logLikelihoodMutationProposed[i] - logLikelihoodMutationCurrent[i];
    }
    unsigned j = 0;
    for(unsigned i = numMutationCat; i < (numMutationCat + numSelectionCat); i++, j++)
    {
        logAcceptanceRatioPerCategory[i] = logLikelihoodSelectionProposed[j] - logLikelihoodSelectionCurrent[j];
    }
    //std::cout << "mutation: " << logAcceptanceRatioPerCategory[0] << "\t selection: " << logAcceptanceRatioPerCategory[1] << std::endl;
}

void ROCModel::obtainCodonCount(SequenceSummary& seqsum, char curAA, int codonCount[])
{
    unsigned* codonRange = SequenceSummary::AAToCodonRange(curAA);
    // get codon counts for AA
    unsigned j = 0u;
    for(unsigned i = codonRange[0]; i < codonRange[1]; i++, j++)
    {
        codonCount[j] = seqsum.getCodonCountForCodon(i);
    }
}





