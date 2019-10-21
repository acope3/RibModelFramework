#include "include/PANSE/PANSEModel.h"

//R runs only
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

//--------------------------------------------------//
//----------- Constructors & Destructors ---------- //
//--------------------------------------------------//

PANSEModel::PANSEModel(unsigned _RFPCountColumn) : Model()
{
    parameter = NULL;
    RFPCountColumn = _RFPCountColumn - 1;
    //ctor
}


PANSEModel::~PANSEModel()
{
    //dtor
    //TODO: call Parent's deconstructor
    //delete parameter;
}


double PANSEModel::calculateLogLikelihoodPerCodonPerGene(double currAlpha, double currLambdaPrime,
        unsigned currRFPObserved, double phiValue, double prevSigma)
{

    // unsigned y_g_fact = 0;
    // for (unsigned j = 1; j <= currRFPObserved;j++)
    // {
    //     y_g_fact += std::log(j);
    // }
	//my_print("Sigma %\n",prevSigma);
    double term1 = std::lgamma(currAlpha + currRFPObserved) - (std::lgamma(currAlpha));// + y_g_fact);
    double term2 = std::log(phiValue) + std::log(prevSigma) - std::log(currLambdaPrime + (phiValue * prevSigma));
    double term3 = std::log(currLambdaPrime) - std::log(currLambdaPrime + (phiValue * prevSigma));

    //my_print("term 1: % term 2: % term 3: %\n",term1,term2,term3);
    term2 *= currRFPObserved;
    term3 *= currAlpha;

    //my_print("term 1: % term 2: % term 3: %\n\n",term1,term2,term3);
    double rv = term1 + term2 + term3;
   // my_print("LogProbability %\n",rv);
    return rv;
}





//------------------------------------------------//
//---------- Likelihood Ratio Functions ----------//
//------------------------------------------------//


void PANSEModel::calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k, double* logProbabilityRatio)
{

    double logLikelihood = 0.0;
    double logLikelihood_proposed = 0.0;


    unsigned alphaCategory = parameter->getMutationCategory(k);
    unsigned lambdaPrimeCategory = parameter->getSelectionCategory(k);
    unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(k);
    double currU = getPartitionFunction(k, false);
    double currSigma = 0.0;
    double totalSigma = 0.0;

    std::vector<double> codonSigmas;
    codonSigmas.resize(getGroupListSize(),INFINITY);

    double phiValue = parameter->getSynthesisRate(geneIndex, synthesisRateCategory, false);
    double phiValue_proposed = parameter->getSynthesisRate(geneIndex, synthesisRateCategory, true);

    std::string geneID = gene.getId();
/*
#ifdef _OPENMP
    //#ifndef __APPLE__
#pragma omp parallel for reduction(+:logLikelihood,logLikelihood_proposed)
#endif
    for (unsigned index = 0; index < getGroupListSize(); index++) //number of codons, without the stop codons
    {
        std::string codon = getGrouping(index);

        double currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
        double currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, false);
        unsigned currRFPObserved = gene.geneData.getCodonSpecificSumRFPCount(index, RFPCountColumn);

        unsigned currNumCodonsInMRNA = gene.geneData.getCodonCountForCodon(index);
        if (currNumCodonsInMRNA == 0) continue;

        logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, currRFPObserved, currNumCodonsInMRNA, phiValue);
        logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, currRFPObserved, currNumCodonsInMRNA, phiValue_proposed);
    }*/

    std::vector <unsigned> positions = gene.geneData.getPositionCodonID();
    std::vector <int> rfpCounts = gene.geneData.getRFPCount(0);


    for (unsigned index = 0; index < positions.size();index++)
    {
        std::string codon = gene.geneData.indexToCodon(positions[index]);
        int codonIndex = positions[index];
        double currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
        double currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, false);
        double currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);
        if(std::isinf(codonSigmas[codonIndex]))
        {
            // This was previously 1 - currAlpha, but looks like this should be currAlpha
            codonSigmas[codonIndex] = elongationProbabilityLog(currAlpha, currLambdaPrime/currU, 1/currNSERate);
        }
        totalSigma += codonSigmas[codonIndex];
    }
// #ifdef _OPENMP
//     //#ifndef __APPLE__
// #pragma omp parallel for private() reduction(+:logLikelihood,logLikelihood_proposed)
// #endif
    for (unsigned index = 0; index < positions.size(); index++)
    {
        std::string codon = gene.geneData.indexToCodon(positions[index]);
        int codonIndex = positions[index];
        std::string codon = gene.geneData.indexToCodon(positions[index]);
        int codonIndex = positions[index];
        double currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
        double currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, false);
        double currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);
        if(std::isinf(codonSigmas[codonIndex]))
        {
            // This was previously 1 - currAlpha, but looks like this should be currAlpha
            codonSigmas[codonIndex] = elongationProbabilityLog(currAlpha, currLambdaPrime/currU, 1/currNSERate);
            
            //my_print("Elongation probability of codon % is %",codonSigmas[codonIndex]);
        }

        //Have to redo the math because rfp observed has changed
        totalSigma += codonSigmas[codonIndex];
    }
// #ifdef _OPENMP
//     //#ifndef __APPLE__
// #pragma omp parallel for private() reduction(+:logLikelihood,logLikelihood_proposed)
// #endif
    for (unsigned index = 0; index < positions.size(); index++)
    {
        std::string codon = gene.geneData.indexToCodon(positions[index]);
        int codonIndex = positions[index];

        double currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
        double currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, false);
        double currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);
        //Should be rfp value at position not all of codon
        unsigned currRFPObserved = rfpCounts[index];
        //my_print("Index % RFP %\n",index,currRFPObserved);
        unsigned currNumCodonsInMRNA = gene.geneData.getCodonCountForCodon(codon);

        //sigma *= elongationProbability(currAlpha, currLambdaPrime/U, 1/currNSERate);
        //This line will never execute
        if (currNumCodonsInMRNA == 0) continue;
        
        logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, currRFPObserved, phiValue/std::exp(totalSigma),std::exp(currSigma));
        logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, currRFPObserved, phiValue_proposed/std::exp(totalSigma), std::exp(currSigma));

        // floating point representation makes it nearly impossible to use == correctly to compare floats.
        // If comparing floats like this, should use (x - y) < epsilon, where epsilon is some small number determined to fit your needs
        // Instead, I'm initializing this vector with infinity and checking if the value is infinite, update accordingly
        if(std::isinf(codonSigmas[codonIndex]))
        {
            // This was previously 1 - currAlpha, but looks like this should be currAlpha
            codonSigmas[codonIndex] = elongationProbabilityLog(currAlpha, currLambdaPrime/currU, 1/currNSERate);

        }
        // Want currSigma to reflect up to position i-1 when calculating loglikelihood for position i, update after loglikelihood calculation
 		currSigma += codonSigmas[codonIndex];
    }

    //Double check math here
    double stdDevSynthesisRate = parameter->getStdDevSynthesisRate(lambdaPrimeCategory, false);
    double logPhiProbability = Parameter::densityLogNorm(phiValue, (-(stdDevSynthesisRate * stdDevSynthesisRate) * 0.5), stdDevSynthesisRate, true);
    double logPhiProbability_proposed = Parameter::densityLogNorm(phiValue_proposed, (-(stdDevSynthesisRate * stdDevSynthesisRate) * 0.5), stdDevSynthesisRate, true);
	double currentLogPosterior = (logLikelihood + logPhiProbability);
	double proposedLogPosterior = (logLikelihood_proposed + logPhiProbability_proposed);

	logProbabilityRatio[0] = (proposedLogPosterior - currentLogPosterior) - (std::log(phiValue) - std::log(phiValue_proposed));//Is recalulcated in MCMC
	logProbabilityRatio[1] = currentLogPosterior - std::log(phiValue_proposed);
	logProbabilityRatio[2] = proposedLogPosterior - std::log(phiValue);
	logProbabilityRatio[3] = currentLogPosterior;
	logProbabilityRatio[4] = proposedLogPosterior;
	logProbabilityRatio[5] = logLikelihood;
	logProbabilityRatio[6] = logLikelihood_proposed;
    //my_print("Gene: Proposed % Current %\n",logLikelihood_proposed,logLikelihood);
}


void PANSEModel::calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome, std::vector<double> &logAcceptanceRatioForAllMixtures)
{
    double logLikelihood = 0.0;
    double logLikelihood_proposed = 0.0;
    double logLikelihood_adjusted = 0.0;
    double logLikelihood_proposed_adjusted = 0.0;
    double currAdjustmentTerm = 0;
    double propAdjustmentTerm = 0;
    Gene *gene;
    unsigned index = SequenceSummary::codonToIndex(grouping);

    double propCodonSigma = INFINITY;
    double propSigma, currSigma,currTotalSigma,propTotalSigma;
    std::vector<double> codonSigmas;

    codonSigmas.resize(getGroupListSize(),INFINITY);

#ifdef _OPENMP
    //#ifndef __APPLE__
#pragma omp parallel for private(gene,currSigma,propSigma,currTotalSigma,propTotalSigma) reduction(+:logLikelihood,logLikelihood_proposed)
#endif
    for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
    {
        gene = &genome.getGene(i);
        currSigma = 0.0;
        propSigma = 0.0;
        currTotalSigma = 0.0;
        propTotalSigma = 0.0;

        // which mixture element does this gene belong to
        unsigned mixtureElement = parameter->getMixtureAssignment(i);

        // how is the mixture element defined. Which categories make it up
        unsigned alphaCategory = parameter->getMutationCategory(mixtureElement);
        unsigned lambdaPrimeCategory = parameter->getSelectionCategory(mixtureElement);
        unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(mixtureElement);
        // get non codon specific values, calculate likelihood conditional on these
        double phiValue = parameter->getSynthesisRate(i, synthesisRateCategory, false);
        unsigned currNumCodonsInMRNA = gene->geneData.getCodonCountForCodon(index);
        double currU = getPartitionFunction(mixtureElement, false);
        if (currNumCodonsInMRNA == 0) continue;

        std::vector <unsigned> positions = gene->geneData.getPositionCodonID();
        std::vector <int> rfpCounts = gene->geneData.getRFPCount(/*RFPCountColumn*/ 0);
        for (unsigned positionIndex = 0; positionIndex < positions.size();positionIndex++)
        {
            
            int codonIndex = positions[positionIndex];
            std::string codon = gene->geneData.indexToCodon(codonIndex);
            double currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
            double currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, false);
            double currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);
            double propAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, true);
            double propLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, true);
            double propNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, true);
            if(std::isinf(codonSigmas[codonIndex]))
            {
                // This was previously 1 - currAlpha, but looks like this should be currAlpha
                codonSigmas[codonIndex] = elongationProbabilityLog(currAlpha, currLambdaPrime/currU, 1/currNSERate);
                //my_print("Elongation probability of codon % is %",codonSigmas[codonIndex]);
            }

            //Have to redo the math because rfp observed has changed
            currTotalSigma += codonSigmas[codonIndex];
            if(codon == grouping)
            {
                if (std::isinf(propCodonSigma)) 
                {   
                    propCodonSigma = elongationProbabilityLog(propAlpha, propLambdaPrime/currU, 1/propNSERate);
                }
                propTotalSigma += propCodonSigma;
                
            }
            else
            {
                propTotalSigma += codonSigmas[codonIndex];
            }

        }
        for (unsigned positionIndex = 0; positionIndex < positions.size(); positionIndex++)
        {

            int positionalRFPCount = rfpCounts[positionIndex];
            int codonIndex = positions[positionIndex];

            //my_print("%\n",positionalRFPCount);
            std::string codon = gene->geneData.indexToCodon(codonIndex);

            double currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
            double currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, false);
            double currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);
            
            if (codon == grouping)
            {
              double propAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, true); 
              double propLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, true);
              double propNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, true);
              logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(propAlpha, propLambdaPrime, positionalRFPCount, phiValue/std::exp(propTotalSigma), std::exp(propSigma));
              if (std::isinf(propCodonSigma)) 
              {   
                propCodonSigma = elongationProbabilityLog(propAlpha, propLambdaPrime/currU, 1/propNSERate);
              }
              propSigma += propCodonSigma;
            }
            else
            {
              logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, positionalRFPCount, phiValue/std::exp(propTotalSigma), std::exp(propSigma));
              propSigma += codonSigmas[codonIndex];
            }
            
            // double propAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, true); 
            // double propLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, true);
            // double propNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, true);
              
        
            //Based on equation 25 in the write-up, at position i, the sigma used should be up to position i + 1
			
            logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, positionalRFPCount, phiValue/std::exp(currTotalSigma), std::exp(currSigma));
            //logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(propAlpha, propLambdaPrime, positionalRFPCount, phiValue/std::exp(propTotalSigma), std::exp(propSigma));
              
            if(std::isinf(codonSigmas[codonIndex]))
            {
                // This was previously 1 - currAlpha, but looks like this should be currAlpha
                codonSigmas[codonIndex] = elongationProbabilityLog(currAlpha, currLambdaPrime/currU, 1/currNSERate);
            }
            currSigma += codonSigmas[codonIndex];
            // if (codon == grouping)
            // {
            //   if (std::isinf(propCodonSigma)) 
            //   {   
            //     propCodonSigma = elongationProbabilityLog(propAlpha, propLambdaPrime/currU, 1/propNSERate);
            //   }
            //   propSigma += propCodonSigma;
            // }
            // else
            // {
            //   propSigma += codonSigmas[codonIndex];
            // }
        }

    }

    unsigned n = getNumMixtureElements();
    for (unsigned j = 0; j < n; j++)
    {
        unsigned alphaCategory = parameter->getMutationCategory(j);
        unsigned lambdaPrimeCategory = parameter->getSelectionCategory(j);
        double currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, grouping, false);
        double currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, grouping, false);
        double currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, grouping, false);
        double propAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, grouping, true);
        double propLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, grouping, true);
        double propNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, grouping, true);
        currAdjustmentTerm += std::log(currAlpha) + std::log(currLambdaPrime) + std::log(currNSERate);
        propAdjustmentTerm += std::log(propAlpha) + std::log(propLambdaPrime) + std::log(propNSERate);
    }
    logAcceptanceRatioForAllMixtures[0] = logLikelihood_proposed - logLikelihood - (currAdjustmentTerm - propAdjustmentTerm);
	logAcceptanceRatioForAllMixtures[1] = logLikelihood - propAdjustmentTerm;
	logAcceptanceRatioForAllMixtures[2] = logLikelihood_proposed - currAdjustmentTerm;
	logAcceptanceRatioForAllMixtures[3] = logLikelihood; //Currently, no prior on CSP parameters for PANSE, so just return loglikelihood
	logAcceptanceRatioForAllMixtures[4] = logLikelihood_proposed;
   // my_print("CSP: Proposed % Current %\n",logLikelihood_proposed,logLikelihood);
}


void PANSEModel::calculateLogLikelihoodRatioForHyperParameters(Genome &genome, unsigned iteration, std::vector <double> & logProbabilityRatio)
{

    double lpr = 0.0; // this variable is only needed because OpenMP doesn't allow variables in reduction clause to be reference

    unsigned selectionCategory = getNumSynthesisRateCategories();
    std::vector<double> currentStdDevSynthesisRate(selectionCategory, 0.0);
    std::vector<double> currentMphi(selectionCategory, 0.0);
    std::vector<double> proposedStdDevSynthesisRate(selectionCategory, 0.0);
    std::vector<double> proposedMphi(selectionCategory, 0.0);
    for (unsigned i = 0u; i < selectionCategory; i++)
    {
        currentStdDevSynthesisRate[i] = getStdDevSynthesisRate(i, false);
        currentMphi[i] = -((currentStdDevSynthesisRate[i] * currentStdDevSynthesisRate[i]) / 2);
        proposedStdDevSynthesisRate[i] = getStdDevSynthesisRate(i, true);
        proposedMphi[i] = -((proposedStdDevSynthesisRate[i] * proposedStdDevSynthesisRate[i]) / 2);
        // take the Jacobian into account for the non-linear transformation from logN to N distribution
        lpr -= (std::log(currentStdDevSynthesisRate[i]) - std::log(proposedStdDevSynthesisRate[i]));
        // take prior into account
        //TODO(Cedric): make sure you can control that prior from R
        //lpr -= Parameter::densityNorm(currentStdDevSynthesisRate[i], 1.0, 0.1, true) - Parameter::densityNorm(proposedStdDevSynthesisRate[i], 1.0, 0.1, true);
    }


    logProbabilityRatio.resize(2);
// #ifdef _OPENMP
//     //#ifndef __APPLE__
// #pragma omp parallel for reduction(+:lpr)
// #endif
    for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
    {
        unsigned mixture = getMixtureAssignment(i);
        mixture = getSynthesisRateCategory(mixture);
        double phi = getSynthesisRate(i, mixture, false);
        lpr += Parameter::densityLogNorm(phi, proposedMphi[mixture], proposedStdDevSynthesisRate[mixture], true) -
            Parameter::densityLogNorm(phi, currentMphi[mixture], currentStdDevSynthesisRate[mixture], true);
    }

    logProbabilityRatio[0] = lpr;

    Gene *gene;
    std::vector <double> currCodonSigmas;
    std::vector <double> propCodonSigmas;
    currCodonSigmas.resize(getGroupListSize(), INFINITY);
    propCodonSigmas.resize(getGroupListSize(), INFINITY);
    
    lpr = 0.0;
    double logLikelihood = 0.0;
    double logLikelihood_proposed = 0.0;
    for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
    {
        
        gene = &genome.getGene(i);
        double currSigma = 0.0;
        double propSigma = 0.0;
        double currTotalSigma = 0.0;
        double propTotalSigma = 0.0;
        // which mixture element does this gene belong to
        unsigned mixtureElement = parameter->getMixtureAssignment(i);
        double currU = getPartitionFunction(mixtureElement, false);
        double propU = getPartitionFunction(mixtureElement, true);
        // how is the mixture element defined. Which categories make it up
        unsigned alphaCategory = parameter->getMutationCategory(mixtureElement);
        unsigned lambdaPrimeCategory = parameter->getSelectionCategory(mixtureElement);
        unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(mixtureElement);
        // get non codon specific values, calculate likelihood conditional on these
        double phiValue = parameter->getSynthesisRate(i, synthesisRateCategory, false);

        std::vector <unsigned> positions = gene->geneData.getPositionCodonID();
        std::vector <int> rfpCounts = gene->geneData.getRFPCount(/*RFPCountColumn*/ 0);

        for (unsigned positionIndex = 0; positionIndex < positions.size(); positionIndex++)
        {

            int codonIndex = positions[positionIndex];
            std::string codon = gene->geneData.indexToCodon(codonIndex);
            double currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
            double currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, false);
            double currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);

            if(std::isinf(currCodonSigmas[codonIndex]))
            {
                currCodonSigmas[codonIndex] = elongationProbabilityLog(currAlpha, currLambdaPrime/currU, 1/currNSERate);
                propCodonSigmas[codonIndex] = elongationProbabilityLog(currAlpha, currLambdaPrime/propU, 1/currNSERate);
            }
            currTotalSigma += currCodonSigmas[codonIndex];
            propTotalSigma += propCodonSigmas[codonIndex];
        }


        for (unsigned positionIndex = 0; positionIndex < positions.size(); positionIndex++)
        {
            int positionalRFPCount = rfpCounts[positionIndex];
            int codonIndex = positions[positionIndex];
            std::string codon = gene->geneData.indexToCodon(codonIndex);
            double currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
            double currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, false);
            double currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);

            logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, positionalRFPCount, phiValue/std::exp(currTotalSigma), std::exp(currSigma));
           
            logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, positionalRFPCount, phiValue/std::exp(propTotalSigma), std::exp(propSigma));


            if(std::isinf(currCodonSigmas[codonIndex]))
            {
                currCodonSigmas[codonIndex] = elongationProbabilityLog(currAlpha, currLambdaPrime/currU, 1/currNSERate);
                propCodonSigmas[codonIndex] = elongationProbabilityLog(currAlpha, currLambdaPrime/propU, 1/currNSERate);
            }
            currSigma += currCodonSigmas[codonIndex];
            propSigma += propCodonSigmas[codonIndex];
        }

        
    }
    unsigned n = getNumMixtureElements();
    for (unsigned j = 0; j < n; j++)
    {
        lpr -= (std::log(getPartitionFunction(j, false)) - std::log(getPartitionFunction(j, true)));
    }
    lpr += logLikelihood_proposed - logLikelihood;
   // my_print("Hyperparameter: Proposed % Current %\n",logLikelihood_proposed,logLikelihood);
    logProbabilityRatio[1] = lpr;
}





//----------------------------------------------------------//
//---------- Initialization and Restart Functions ----------//
//----------------------------------------------------------//


void PANSEModel::initTraces(unsigned samples, unsigned num_genes, bool estimateSynthesisRate)
{
    parameter->initAllTraces(samples, num_genes, estimateSynthesisRate);
}


void PANSEModel::writeRestartFile(std::string filename)
{
    return parameter->writeEntireRestartFile(filename);
}





//----------------------------------------//
//---------- Category Functions ----------//
//----------------------------------------//


double PANSEModel::getCategoryProbability(unsigned i)
{
    return parameter->getCategoryProbability(i);
}


unsigned PANSEModel::getMutationCategory(unsigned mixture)
{
    return parameter->getMutationCategory(mixture);
}


unsigned PANSEModel::getSelectionCategory(unsigned mixture)
{
    return parameter->getSelectionCategory(mixture);
}


unsigned PANSEModel::getSynthesisRateCategory(unsigned mixture)
{
    return parameter->getSynthesisRateCategory(mixture);
}


std::vector<unsigned> PANSEModel::getMixtureElementsOfSelectionCategory(unsigned k)
{
    return parameter->getMixtureElementsOfSelectionCategory(k);
}





//------------------------------------------//
//---------- Group List Functions ----------//
//------------------------------------------//


unsigned PANSEModel::getGroupListSize()
{
    return parameter->getGroupListSize();
}


std::string PANSEModel::getGrouping(unsigned index)
{
    return parameter->getGrouping(index);
}





//---------------------------------------------------//
//---------- stdDevSynthesisRate Functions ----------//
//---------------------------------------------------//


double PANSEModel::getStdDevSynthesisRate(unsigned selectionCategory, bool proposed)
{
    return parameter->getStdDevSynthesisRate(selectionCategory, proposed);
}


double PANSEModel::getCurrentStdDevSynthesisRateProposalWidth()
{
    return parameter->getCurrentStdDevSynthesisRateProposalWidth();
}



void PANSEModel::updateStdDevSynthesisRate()
{
    parameter->updateStdDevSynthesisRate();
}


//---------------------------------------------------//
//----------- partitionFunction Functions -----------//
//---------------------------------------------------//


double PANSEModel::getPartitionFunction(unsigned selectionCategory, bool proposed)
{
    return parameter->getPartitionFunction(selectionCategory, proposed);
}


double PANSEModel::getCurrentPartitionFunctionProposalWidth()
{
    return parameter->getCurrentPartitionFunctionProposalWidth();
}



void PANSEModel::updatePartitionFunction()
{
    parameter->updatePartitionFunction();
}



//----------------------------------------------//
//---------- Synthesis Rate Functions ----------//
//----------------------------------------------//


double PANSEModel::getSynthesisRate(unsigned index, unsigned mixture, bool proposed)
{
    return parameter->getSynthesisRate(index, mixture, proposed);
}


void PANSEModel::updateSynthesisRate(unsigned i, unsigned k)
{
    parameter->updateSynthesisRate(i, k);
}





//-----------------------------------------//
//---------- Iteration Functions ----------//
//-----------------------------------------//


unsigned PANSEModel::getLastIteration()
{
    return parameter->getLastIteration();
}


void PANSEModel::setLastIteration(unsigned iteration)
{
    parameter->setLastIteration(iteration);
}





//-------------------------------------//
//---------- Trace Functions ----------//
//-------------------------------------//


void PANSEModel::updateStdDevSynthesisRateTrace(unsigned sample)
{
    parameter->updateStdDevSynthesisRateTrace(sample);
}


void PANSEModel::updatePartitionFunctionTrace(unsigned sample)
{
    parameter->updatePartitionFunctionTrace(sample);
}


void PANSEModel::updateSynthesisRateTrace(unsigned sample, unsigned i)
{
    parameter->updateSynthesisRateTrace(sample, i);
}


void PANSEModel::updateMixtureAssignmentTrace(unsigned sample, unsigned i)
{
    parameter->updateMixtureAssignmentTrace(sample, i);
}


void PANSEModel::updateMixtureProbabilitiesTrace(unsigned sample)
{
    parameter->updateMixtureProbabilitiesTrace(sample);
}


void PANSEModel::updateCodonSpecificParameterTrace(unsigned sample, std::string codon)
{
    parameter->updateCodonSpecificParameterTrace(sample, codon);
}


void PANSEModel::updateHyperParameterTraces(unsigned sample)
{
    updateStdDevSynthesisRateTrace(sample);
    updatePartitionFunctionTrace(sample);
}


void PANSEModel::updateTracesWithInitialValues(Genome & genome)
{
    std::vector <std::string> groupList = parameter->getGroupList();

    for (unsigned i = 0; i < genome.getGenomeSize(); i++)
    {
        parameter->updateSynthesisRateTrace(0, i);
        parameter->updateMixtureAssignmentTrace(0, i);
    }

    for (unsigned i = 0; i < groupList.size(); i++)
    {
        parameter->updateCodonSpecificParameterTrace(0, getGrouping(i));
    }
}





//----------------------------------------------//
//---------- Adaptive Width Functions ----------//
//----------------------------------------------//


void PANSEModel::adaptStdDevSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt)
{
    parameter->adaptStdDevSynthesisRateProposalWidth(adaptiveWidth, adapt);
}


void PANSEModel::adaptPartitionFunctionProposalWidth(unsigned adaptiveWidth, bool adapt)
{
    parameter->adaptPartitionFunctionProposalWidth(adaptiveWidth, adapt);
}


void PANSEModel::adaptSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt)
{
    parameter->adaptSynthesisRateProposalWidth(adaptiveWidth, adapt);
}


void PANSEModel::adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth, unsigned lastIteration, bool adapt)
{
    parameter->adaptCodonSpecificParameterProposalWidth(adaptiveWidth, lastIteration, adapt);
}


void PANSEModel::adaptHyperParameterProposalWidths(unsigned adaptiveWidth, bool adapt)
{
    adaptStdDevSynthesisRateProposalWidth(adaptiveWidth, adapt);
    adaptPartitionFunctionProposalWidth(adaptiveWidth, adapt);
}



//-------------------------------------//
//---------- Other Functions ----------//
//-------------------------------------//


void PANSEModel::proposeCodonSpecificParameter()
{
    parameter->proposeCodonSpecificParameter();
}


void PANSEModel::proposeHyperParameters()
{
    parameter->proposeStdDevSynthesisRate();
    parameter->proposePartitionFunction();
}


void PANSEModel::proposeSynthesisRateLevels()
{
    parameter->proposeSynthesisRateLevels();
}


unsigned PANSEModel::getNumPhiGroupings()
{
    return parameter->getNumObservedPhiSets();
}


unsigned PANSEModel::getMixtureAssignment(unsigned index)
{
    return parameter->getMixtureAssignment(index);
}


unsigned PANSEModel::getNumMixtureElements()
{
    return parameter->getNumMixtureElements();
}


unsigned PANSEModel::getNumSynthesisRateCategories()
{
    return parameter->getNumSynthesisRateCategories();
}


void PANSEModel::setNumPhiGroupings(unsigned value)
{
    parameter->setNumObservedPhiSets(value);
}


void PANSEModel::setMixtureAssignment(unsigned i, unsigned catOfGene)
{
    parameter->setMixtureAssignment(i, catOfGene);
}


void PANSEModel::setCategoryProbability(unsigned mixture, double value)
{
    parameter->setCategoryProbability(mixture, value);
}


void PANSEModel::updateCodonSpecificParameter(std::string aa)
{
    parameter->updateCodonSpecificParameter(aa);
}


void PANSEModel::updateGibbsSampledHyperParameters(Genome &genome)
{
    //TODO: fill in
}


void PANSEModel::updateAllHyperParameter()
{
    updateStdDevSynthesisRate();
    updatePartitionFunction();
}


void PANSEModel::updateHyperParameter(unsigned hp)
{
    // NOTE: when adding additional hyper parameter, also add to updateAllHyperParameter()
    switch (hp)
    {
        case 0:
            updateStdDevSynthesisRate();
            break;
        case 1:
            updatePartitionFunction();
            break;
        default:
            break;
    }
}

void PANSEModel::simulateGenome(Genome &genome)
{
    for (unsigned geneIndex = 0u; geneIndex < genome.getGenomeSize(); geneIndex++)
    {
        unsigned mixtureElement = getMixtureAssignment(geneIndex);
        Gene gene = genome.getGene(geneIndex);
        double phi = parameter->getSynthesisRate(geneIndex, mixtureElement, false);
        SequenceSummary sequence = gene.geneData;
        Gene tmpGene = gene;
        std::vector <unsigned> positions = sequence.getPositionCodonID();
        std::vector <int> rfpCount;
        unsigned alphaCat = parameter->getMutationCategory(mixtureElement);
        unsigned lambdaPrimeCat = parameter->getSelectionCategory(mixtureElement);
        double sigma =  1.0;
        double v;
        for (unsigned codonID : positions)
        {
            std::string codon = SequenceSummary::codonArray[codonID];
            double alpha = getParameterForCategory(alphaCat, PANSEParameter::alp, codon, false);
            double lambdaPrime = getParameterForCategory(lambdaPrimeCat, PANSEParameter::lmPri, codon, false);
            double NSERate = getParameterForCategory(alphaCat, PANSEParameter::nse, codon, false);

            if (NSERate == 0){v = 1000000000;}
            else {v = 1.0 / NSERate;}

#ifndef STANDALONE
            RNGScope scope;
            NumericVector xx(1);
            xx = rgamma(1, alpha, 1.0/lambdaPrime);
            sigma *= (v/(xx[0] + v));
            xx = rpois(1, xx[0] * phi * sigma);
            rfpCount.push_back(xx[0]);
#else
            std::gamma_distribution<double> GDistribution(alphaPrime, 1.0/lambdaPrime);
            double tmp = GDistribution(Parameter::generator);
            sigma *= (v/(tmp + v));
            std::poisson_distribution<unsigned> PDistribution(phi * tmp * sigma);
            unsigned simulatedValue = PDistribution(Parameter::generator);
            rfpCount.push_back(simulatedValue);
#endif
        }
        tmpGene.geneData.setRFPCount(rfpCount, RFPCountColumn);
        genome.addGene(tmpGene, true);
    }
}


void PANSEModel::printHyperParameters()
{
    for (unsigned i = 0u; i < getNumSynthesisRateCategories(); i++)
    {
        my_print("stdDevSynthesisRate posterior estimate for selection category %: %\n", i, parameter->getStdDevSynthesisRate(i));
    }
    my_print("\t current stdDevSynthesisRate proposal width: %\n", getCurrentStdDevSynthesisRateProposalWidth());
}


/* getParameter (RCPP EXPOSED)
 * Arguments: None
 *
 * Returns the PANSEParameter of the model.
 */
PANSEParameter* PANSEModel::getParameter()
{
    return parameter;
}


void PANSEModel::setParameter(PANSEParameter &_parameter)
{
    parameter = &_parameter;
}


double PANSEModel::calculateAllPriors()
{
    return 0.0; //TODO(Cedric): implement me, see ROCModel
}


double PANSEModel::getParameterForCategory(unsigned category, unsigned param, std::string codon, bool proposal)
{
    return parameter->getParameterForCategory(category, param, codon, proposal);
}

//Continued fractions helper function for upper incomplete gamma
double PANSEModel::UpperIncompleteGammaHelper(double s, double x)
{
    double rv;
    int i;

    rv = 1000.0 / x;
    for(i = 1000; i > 0; i--){
        if(i % 2 == 0) rv = (double) ((int) (i / 2)) / (x + rv);
        if(i % 2 != 0) rv = ((double) (((int) (i / 2) + 1) - s)) / (1 + rv);
    }

    return x + rv;
}

//Upper incomplete gamma function
double PANSEModel::UpperIncompleteGamma(double s, double x)
{
    double d, rv;
    rv = std::pow(x,s) * std::exp(0 - x);
    d = UpperIncompleteGammaHelper(s, x);
    return rv/d;

}

//log of upper incomplete gamma function
double PANSEModel::UpperIncompleteGammaLog(double s, double x)
{
    double rv, d;
    rv = s * std::log(x) - x;
    d = std::log(PANSEModel::UpperIncompleteGammaHelper(s, x));
    return rv - d;

}

//Generalized integral function
double PANSEModel::GeneralizedIntegral(double p, double z){
    return std::pow(z, p - 1.0) * UpperIncompleteGamma(1.0 - p, z);
}

//Log of generalized integral function
double PANSEModel::GeneralizedIntegralLog(double p, double z){
    return (p - 1.0) * std::log(z) + UpperIncompleteGammaLog(1.0 - p, z);
}

//Calculation of the probability of elongation at current codon
double PANSEModel::elongationProbability(double currAlpha, double currLambda, double currNSE)
{

    return std::pow(currLambda * currNSE, currAlpha) * std::exp(currLambda * currNSE) * UpperIncompleteGamma(1 - currAlpha, currLambda * currNSE);

}

//Log probability of elongation at current codon
double PANSEModel::elongationProbabilityLog(double currAlpha, double currLambda, double currNSE){
    return (currLambda * currNSE) + currAlpha * (std::log(currLambda) + std::log(currNSE)) + UpperIncompleteGammaLog(1- currAlpha,currLambda * currNSE);
}

double PANSEModel::elongationUnitilIndexProbability(int index, std::vector <double> lambdas, std::vector <double> NSERates){
    return 0;
}

double PANSEModel::elongationUnitilIndexProbabilityLog(int index, std::vector <double> lambdas, std::vector <double> NSERates){
    return 0;
}

double PANSEModel::prob_Y_g(double curralpha, int sample_size, double lambda_prime, double psi, double prevdelta){
    double term1, term2, term3;

    term1 = std::tgamma(curralpha + sample_size) / std::tgamma(curralpha);
    term2 = psi * prevdelta / (lambda_prime + (psi * prevdelta));
    term3 = lambda_prime / (lambda_prime + (psi * prevdelta));

    term2 = pow(term2, sample_size);
    term3 = pow(term3, curralpha);

    return term2 * term3 * term1;
}

double PANSEModel::prob_Y_g_log(double curralpha, int sample_size, double lambda_prime, double psi, double prevdelta){
    double term1, term2, term3;

    term1 = std::lgamma(curralpha + sample_size) - std::lgamma(curralpha);
    term2 = std::log(psi) + std::log(prevdelta) - std::log(lambda_prime + (psi * prevdelta));
    term3 = std::log(lambda_prime) - std::log(lambda_prime + (psi * prevdelta));

    term2 *= sample_size;
    term3 *= curralpha;

    return term1 + term2 + term3;
}

//TODO: Add sigma to parameter object to keep it from being calulcated and initialize as need using MCMC adjust functions to reflect this
double psi2phi(double psi, double sigma){
    return sigma * psi;
}
double phi2psi(double phi, double sigma){
    return phi / sigma;
}

std::vector <double> readAlphaValues(std::string filename){
    std::size_t pos;
    std::ifstream currentFile;
    std::string tmpString;
    std::vector <double> rv;

    rv.resize(64);

    currentFile.open(filename);
    if (currentFile.fail())
        my_printError("Error opening file %\n", filename.c_str());
    else
    {
        currentFile >> tmpString;
        while (currentFile >> tmpString){
            pos = tmpString.find(',');
            if (pos != std::string::npos)
            {
                std::string codon = tmpString.substr(0,3);
                std::string val = tmpString.substr(pos + 1, std::string::npos);
                rv[SequenceSummary::codonToIndex(codon, true)] = std::atof(val.c_str());
            }
        }
    }

    return rv;

}
std::vector <double> readLambdaValues(std::string filename){
    std::size_t pos;
    std::ifstream currentFile;
    std::string tmpString;
    std::vector <double> rv;

    rv.resize(64);

    currentFile.open(filename);
    if (currentFile.fail())
        my_printError("Error opening file %\n", filename.c_str());
    else
    {
        currentFile >> tmpString;
        while (currentFile >> tmpString){
            pos = tmpString.find(',');
            if (pos != std::string::npos)
            {
                std::string codon = tmpString.substr(0,3);
                std::string val = tmpString.substr(pos + 1, std::string::npos);
                rv[SequenceSummary::codonToIndex(codon, true)] = std::atof(val.c_str());
            }
        }
    }

    return rv;
}
std::vector <double> readNSERateValues(std::string filename){
    std::size_t pos;
    std::ifstream currentFile;
    std::string tmpString;
    std::vector <double> rv;

    rv.resize(64);

    currentFile.open(filename);
    if (currentFile.fail())
        my_printError("Error opening file %\n", filename.c_str());
    else
    {
        currentFile >> tmpString;
        while (currentFile >> tmpString){
            pos = tmpString.find(',');
            if (pos != std::string::npos)
            {
                std::string codon = tmpString.substr(0,3);
                std::string val = tmpString.substr(pos + 1, std::string::npos);
                rv[SequenceSummary::codonToIndex(codon, true)] = std::atof(val.c_str());
            }
        }
    }

    return rv;
}
