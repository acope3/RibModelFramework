#include "../include/ROCModel.h"

#include <vector>
#include <math.h>
#include <cfloat>
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
		int numParam = parameter.getNumParam();
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
            logLikelihood += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation, selection, phiValue, numParam);
						//std::cout <<"calling now!\n";
            logLikelihood_proposed += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation, selection, phiValue_proposed, numParam);
						//if (std::isnan(logLikelihood_proposed)) std::cout <<"WRONG HERE\n";
						//if (!std::isfinite(logLikelihood_proposed)) std::cout <<"Finished AA: " << i <<"\n";
						//else std::cout <<"At AA: " << i <<" logLikelihood_proposed is " << logLikelihood_proposed <<"\n";
				}
    }

    double sPhi = parameter.getSphi(false);
    double logPhiProbability = std::log(ROCParameter::densityLogNorm(phiValue, (-(sPhi * sPhi) / 2), sPhi));
    double logPhiProbability_proposed = std::log(ROCParameter::densityLogNorm(phiValue_proposed, (-(sPhi * sPhi) / 2), sPhi));
    double currentLogLikelihood = (logLikelihood + logPhiProbability);
    double proposedLogLikelihood = (logLikelihood_proposed + logPhiProbability_proposed);
	/*	if (std::isnan(proposedLogLikelihood))
		{
			std::cout <<"proposedLogLikelihood: " << proposedLogLikelihood <<"\n";
			std::cout <<"logLikelihood_proposed: " << logLikelihood_proposed <<"\n";
			std::cout <<"logPhiProbability_proposed: " << logPhiProbability_proposed <<"\n";
		}
		*/
    logProbabilityRatio[0] = (proposedLogLikelihood - currentLogLikelihood) - (std::log(phiValue) - std::log(phiValue_proposed));
    logProbabilityRatio[1] = currentLogLikelihood - std::log(phiValue_proposed);
    logProbabilityRatio[2] = proposedLogLikelihood - std::log(phiValue);
		/*if (std::isnan(logProbabilityRatio[2]))
		{
			std::cout <<"logProbabilityRatio[2]: " << logProbabilityRatio[2] <<"\n";
			std::cout <<"proposedLogLikelihood: " << proposedLogLikelihood <<"\n";
			std::cout <<"std::log(phiValue): " << std::log(phiValue) <<"\n";
			std::cout <<"phiValue: " << phiValue << "\n\n\n";
		}*/
    //return logLikelihood + a;
}
/*
void ROCModel::calculateCodonProbabilityVector(unsigned numCodons, double mutation[], double selection[], double phi, double codonProb[])
{
    // calculate numerator and denominator for codon probabilities
    double denominator = 1.0;
		double max = (DBL_MAX - 5000 ) / 6.0;
    for(unsigned i = 0; i < (numCodons - 1); i++)
    {
        codonProb[i] = std::exp( -mutation[i] - (selection[i] * phi) );
        if (codonProb[i] > max)
				{
					codonProb[i] = max;
				}
				denominator += codonProb[i];
    }
		if (std::isinf(denominator)) std::cout <<"At infinity\n";
    // alphabetically last codon is reference codon!
    codonProb[numCodons - 1] = 1.0;
    // normalize codon probabilities
    for(unsigned i = 0; i < numCodons; i++)
		{
			if (std::isinf(denominator)) std::cout <<"BEFORE: " << codonProb[i] <<"\n";
			codonProb[i] = codonProb[i] / denominator;
		  if (std::isinf(denominator))
			{
				std::cout <<"codonProb = " << codonProb[i] <<"\n";
			}
		}
}
*/
void ROCModel::calculateCodonProbabilityVector(unsigned numCodons, double mutation[], double selection[], double phi, double codonProb[], int numParam)
{
    // calculate numerator and denominator for codon probabilities
    double denominator = 0.0;
		int minIndexVal = 0;
		//double oldProb[numCodons];

		for (unsigned i = 1u; i < numCodons - 1; i++)
		{
			if (selection[minIndexVal] > selection[i])
			{
				minIndexVal = i;
			}
		}

    for(unsigned i = 0; i < (numCodons - 1); i++)
    {
        codonProb[i] = std::exp( -(mutation[i] - mutation[minIndexVal]) - ((selection[i] - selection[minIndexVal]) * phi) );
     //   oldProb[i] = std::exp( -mutation[i] - (selection[i] * phi) );
				denominator += codonProb[i];
				if (std::isinf(codonProb[i]))
				{
					std::cout <<" mutation[i]: " <<mutation[i] <<"\n";
					std::cout <<" mutation[min]: " <<mutation[minIndexVal] << "\n";
					std::cout <<"(mutation[i] - mutation[minIndexVal]): " << (mutation[i] - mutation[minIndexVal]) <<"\n";
    			std::cout <<"selection[i] - selection[minIndexVal]: " << selection[i] - selection[minIndexVal] <<"\n";
					std::cout <<"phi : " << phi <<"\n";
				}
		}

    // alphabetically last codon is reference codon!
    codonProb[numCodons - 1] = std::exp(mutation[minIndexVal] + selection[minIndexVal] * phi);
    denominator += codonProb[numCodons - 1];
		if (std::isinf(denominator)) std::cout <<"infinity\n";
		//oldProb[numCodons - 1] = 1.0;

		// normalize codon probabilities
    for(unsigned i = 0; i < numCodons; i++)
		{
			codonProb[i] = codonProb[i] / denominator;
			//oldProb[i] = oldProb[i] / oldDenominator;
		//	std::cout << "codonProb vs oldProb: " << codonProb[i] << " vs " << oldProb[i] <<"\n";
		}
}

double ROCModel::calculateLogLikelihoodPerAAPerGene(unsigned numCodons, int codonCount[], double mutation[], double selection[], double phiValue, int numParam)
{
    double logLikelihood = 0.0;
    // calculate codon probabilities
    double codonProbabilities[numCodons];
    calculateCodonProbabilityVector(numCodons, mutation, selection, phiValue, codonProbabilities, numParam);

    // calculate likelihood for current AA for this combination of selection and mutation category
    for(unsigned i = 0; i < numCodons; i++)
    {
				if (codonCount[i] == 0) continue;
        logLikelihood += std::log(codonProbabilities[i]) * codonCount[i];
			//if (std::isnan(logLikelihood))
				{
					//std::cout <<"logLikelihood: " << logLikelihood <<"\n";
					//std::cout <<"codonCount[i]: " << codonCount[i] <<"\n";
					//std::cout <<"i:           : " << i <<"\n";
					//std::cout <<"codonProbabilities[i]: " << codonProbabilities[i] <<"\n\n";
				}
    }
    return logLikelihood;
}

void ROCModel::calculateLogLikelihoodRatioPerAAPerCategory(char curAA, Genome& genome, ROCParameter& parameter, double& logAcceptanceRatioForAllMixtures)
{

    int numGenes = genome.getGenomeSize();
    int numCodons = SequenceSummary::GetNumCodonsForAA(curAA);
    double likelihood = 0.0;
    double likelihood_proposed = 0.0;
		int numParam = parameter.getNumParam();
    for(int i = 0; i < numGenes; i++)
    {
        int codonCount[numCodons];
        double mutation[numCodons - 1];
        double selection[numCodons - 1];
        double mutation_proposed[numCodons - 1];
        double selection_proposed[numCodons - 1];
        Gene gene = genome.getGene(i);
        SequenceSummary seqsum = gene.getSequenceSummary();
        if(seqsum.getAAcountForAA(curAA) == 0) continue;

        // which mixture element does this gene belong to
        unsigned mixtureElement = parameter.getMixtureAssignment(i);
        // how is the mixture element defined. Which categories make it up
        unsigned mutationCategory = parameter.getMutationCategory(mixtureElement);
        unsigned selectionCategory = parameter.getSelectionCategory(mixtureElement);
        unsigned expressionCategory = parameter.getExpressionCategory(mixtureElement);
        // get phi value, calculate likelihood conditional on phi
        double phiValue = parameter.getExpression(i, expressionCategory, false);

        // get current mutation and selection parameter
        parameter.getParameterForCategory(mutationCategory, ROCParameter::dM, curAA, false, mutation);
        parameter.getParameterForCategory(selectionCategory, ROCParameter::dEta, curAA, false, selection);

        // get proposed mutation and selection parameter
        parameter.getParameterForCategory(mutationCategory, ROCParameter::dM, curAA, true, mutation_proposed);
        parameter.getParameterForCategory(selectionCategory, ROCParameter::dEta, curAA, true, selection_proposed);

//        if(i == 100 || i == 600)
//        {
//            std::cout << "Gene " << i << ": mixtureElement = " << mixtureElement << "\n";
//            std::cout << "Gene " << i << ": mutCat = " << mutationCategory << ", selCat = " << selectionCategory << ", exprCat = " << expressionCategory << "\n";
//            std::cout << curAA << " selection array: ";
//            for(int k = 0; k < numCodons - 1; k++){std::cout << selection[k] << " ";}
//            std::cout << "\n";
//            std::cout << curAA << " selection proposed array: ";
//            for(int k = 0; k < numCodons - 1; k++){std::cout << selection_proposed[k] << " ";}
//            std::cout << "\n\n";
//
//            std::cout << curAA << " mutation array: ";
//            for(int k = 0; k < numCodons - 1; k++){std::cout << mutation[k] << " ";}
//            std::cout << "\n";
//            std::cout << curAA << " mutation proposed array: ";
//            for(int k = 0; k < numCodons - 1; k++){std::cout << mutation_proposed[k] << " ";}
//            std::cout << "\n\n\n\n";
//        }

        obtainCodonCount(seqsum, curAA, codonCount);
        //std::cout << "Gene " << i << ": mixtureElement = " << mixtureElement << "\t mixtureElementProbability = " << mixtureElementProbability << "\n";
        double a = calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation, selection, phiValue, numParam);
        double b = calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation_proposed, selection_proposed, phiValue, numParam);
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





