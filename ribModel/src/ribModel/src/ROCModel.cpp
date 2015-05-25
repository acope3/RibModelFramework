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

void ROCModel::calculateLogLiklihoodRatioPerGene(Gene& gene, int geneIndex, ROCParameter& parameter, unsigned k, double* logProbabilityRatio)
{
    double logLikelihood = 0.0;
    double logLikelihood_proposed = 0.0;

    SequenceSummary seqsum = gene.getSequenceSummary();

    // get correct index for everything
    unsigned mutationCategory = parameter.getMutationCategory(k);
    unsigned selectionCategory = parameter.getSelectionCategory(k);
    unsigned expressionCategory = parameter.getExpressionCategory(k);

    double phiValue = parameter.getExpression(geneIndex, expressionCategory, false);
    double phiValue_proposed = parameter.getExpression(geneIndex, expressionCategory, true);
	//#pragma omp parallel for
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
        parameter.getParameterForCategory(mutationCategory, ROCParameter::dM, curAA, false, mutation);
        double selection[numCodons - 1];
        parameter.getParameterForCategory(selectionCategory, ROCParameter::dEta, curAA, false, selection);

        // prepare array for codon counts for AA
        int codonCount[numCodons];
        obtainCodonCount(seqsum, curAA, codonCount);

        //#pragma omp parallel num_threads(2)
        {
            logLikelihood += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation, selection, phiValue);
            logLikelihood_proposed += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation, selection, phiValue_proposed);
        }
    }

    double sPhi = parameter.getSphi(false);
    double logPhiProbability = std::log(ROCParameter::densityLogNorm(phiValue, (-(sPhi * sPhi) / 2), sPhi));
    double logPhiProbability_proposed = std::log(ROCParameter::densityLogNorm(phiValue_proposed, (-(sPhi * sPhi) / 2), sPhi));

    double currentLogLikelihood = (logLikelihood + logPhiProbability);
    double proposedLogLikelihood = (logLikelihood_proposed + logPhiProbability_proposed);

    logProbabilityRatio[0] = (proposedLogLikelihood - currentLogLikelihood) - (std::log(phiValue) - std::log(phiValue_proposed));
    logProbabilityRatio[1] = currentLogLikelihood - std::log(phiValue_proposed);
    logProbabilityRatio[2] = proposedLogLikelihood - std::log(phiValue);

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

double ROCModel::calculateLogLikelihoodPerAAPerGene(unsigned numCodons, int codonCount[], double mutation[], double selection[], double phiValue)
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

void ROCModel::calculateLogLikelihoodRatioPerAAPerCategory(char curAA, Genome& genome, ROCParameter& parameter, double& logAcceptanceRatioForAllMixtures)
{

    int numGenes = genome.getGenomeSize();
    double likelihood = 0.0;
    double likelihood_proposed = 0.0;
    int numCodons = SequenceSummary::GetNumCodonsForAA(curAA);

    for(int i = 0; i < numGenes; i++)
    {
        Gene gene = genome.getGene(i);
        SequenceSummary seqsum = gene.getSequenceSummary();
        if(seqsum.getAAcountForAA(curAA) == 0) continue;

        // which mixture element does this gene belong to
        unsigned mixtureElement = parameter.getMixtureAssignment(i);
        //std::cout << "Gene " << i << ": mixtureElement = " << mixtureElement << "\n";
        // how is the mixture element defined. Which categories make it up
        unsigned mutationCategory = parameter.getMutationCategory(mixtureElement);
        unsigned selectionCategory = parameter.getSelectionCategory(mixtureElement);
        unsigned expressionCategory = parameter.getExpressionCategory(mixtureElement);
        //std::cout << "Gene " << i << ": mutCat = " << mutationCategory << ", selCat = " << selectionCategory << ", exprCat = " << expressionCategory << "\n";
        // get phi value, calculate likelihood conditional on phi
        double phiValue = parameter.getExpression(i, expressionCategory, false);

        // get current mutation and selection parameter
        double mutation[numCodons - 1];
        parameter.getParameterForCategory(mutationCategory, ROCParameter::dM, curAA, false, mutation);
        double selection[numCodons - 1];
        parameter.getParameterForCategory(selectionCategory, ROCParameter::dEta, curAA, false, selection);

        // get proposed mutation and selection parameter
        double mutation_proposed[numCodons - 1];
        parameter.getParameterForCategory(mutationCategory, ROCParameter::dM, curAA, true, mutation_proposed);
        double selection_proposed[numCodons - 1];
        parameter.getParameterForCategory(selectionCategory, ROCParameter::dEta, curAA, true, selection_proposed);

/*        if(i == 0)
        {
//            std::cout << curAA << " selection array: ";
//            for(int k = 0; k < numCodons - 1; k++){std::cout << selection[k] << " ";}
//            std::cout << "\n";
//            std::cout << curAA << " selection proposed array: ";
//            for(int k = 0; k < numCodons - 1; k++){std::cout << selection_proposed[k] << " ";}
//            std::cout << "\n";
//
//            std::cout << curAA << " mutation array: ";
//            for(int k = 0; k < numCodons - 1; k++){std::cout << mutation[k] << " ";}
//            std::cout << "\n";
//            std::cout << curAA << " mutation proposed array: ";
//            for(int k = 0; k < numCodons - 1; k++){std::cout << mutation_proposed[k] << " ";}
//            std::cout << "\n";
        }
*/
        int codonCount[numCodons];
        obtainCodonCount(seqsum, curAA, codonCount);

        // get probability of current mixture assignment, calculate likelihood conditional on current mixture assignment
        double mixtureElementProbability = parameter.getCategoryProbability(mixtureElement);
        //std::cout << "Gene " << i << ": mixtureElement = " << mixtureElement << "\t mixtureElementProbability = " << mixtureElementProbability << "\n";
        double a = calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation, selection, phiValue);
        double b = calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation_proposed, selection_proposed, phiValue);
        //std::cout << "curLogLike: " << a << "\t propLogLike: " << b << "\n";
        likelihood += a; //std::exp( calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, seqsum, mutation, selection, phiValue) );
        likelihood_proposed += b; //std::exp( calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, seqsum, mutation_proposed, selection_proposed, phiValue) );
        //std::cout << "curMixedLike: " << likelihood << "\t propMixedLike: " << likelihood_proposed << "\n";
    }
    //std::cout << "( likelihood_proposed / likelihood ) = " << likelihood_proposed / likelihood  << "\n";
    logAcceptanceRatioForAllMixtures = likelihood_proposed - likelihood;
}

void ROCModel::obtainCodonCount(SequenceSummary& seqsum, char curAA, int codonCount[])
{
    unsigned codonRange[2];
    SequenceSummary::AAToCodonRange(curAA, false, codonRange);
    // get codon counts for AA
    unsigned j = 0u;
    //std::cout << curAA << " " << codonRange[0] << " to " << codonRange[1] << "\n";
    for(unsigned i = codonRange[0]; i < codonRange[1]; i++, j++)
    {
        codonCount[j] = seqsum.getCodonCountForCodon(i);
    }
}





