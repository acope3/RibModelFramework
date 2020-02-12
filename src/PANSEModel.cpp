#include "include/PANSE/PANSEModel.h"

//R runs only
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

//--------------------------------------------------//
//----------- Constructors & Destructors ---------- //
//--------------------------------------------------//

PANSEModel::PANSEModel(unsigned _RFPCountColumn, bool _withPhi, bool _fix_sEpsilon) : Model()
{
    parameter = NULL;
    RFPCountColumn = _RFPCountColumn - 1;
    withPhi = _withPhi;
    fix_sEpsilon = _fix_sEpsilon;
    //ctor
}


PANSEModel::~PANSEModel()
{
    //dtor
    //TODO: call Parent's deconstructor
    //delete parameter;
}


double PANSEModel::calculateLogLikelihoodPerCodonPerGene(double currAlpha, double currLambdaPrime,
        unsigned currRFPObserved, double phiValue, double prevSigma, double lgamma_currAlpha, double log_currLambdaPrime, double log_phi,double lgamma_rfp_alpha)
{
    double term1 = lgamma_rfp_alpha - lgamma_currAlpha;//std::lgamma(currAlpha);
    double term2 = log_phi + std::log(prevSigma) - std::log(currLambdaPrime + (phiValue * prevSigma));
    double term3 = log_currLambdaPrime - std::log(currLambdaPrime + (phiValue * prevSigma));

    term2 *= currRFPObserved;
    term3 *= currAlpha;


    double rv = term1 + term2 + term3;
    return rv;
}





//------------------------------------------------//
//---------- Likelihood Ratio Functions ----------//
//------------------------------------------------//


void PANSEModel::calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k, double* logProbabilityRatio)
{
    double currAlpha,currLambdaPrime,currNSERate;
    std::string codon;
    double logLikelihood = 0.0;
    double logLikelihood_proposed = 0.0;

    std::vector <unsigned> positions = gene.geneData.getPositionCodonID();
    std::vector <unsigned> rfpCounts = gene.geneData.getRFPCount(0);

    std::vector<double> sigma(positions.size()+1,0);
    std::vector<double> alpha(positions.size(),0);
    std::vector<double> lambda_prime(positions.size(),0);
    std::vector<double> local_lgamma_currentAlpha(getGroupListSize(), -1000);
    std::vector<double> local_log_currentLambdaPrime(getGroupListSize(), -1000);


    unsigned alphaCategory = parameter->getMutationCategory(k);
    unsigned lambdaPrimeCategory = parameter->getSelectionCategory(k);
    unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(k);
    unsigned mixtureElement = parameter->getMixtureAssignment(geneIndex);
    unsigned Y = parameter->getTotalRFPCount();
    double U = getPartitionFunction(mixtureElement, false)/Y;

    double phiValue = parameter->getSynthesisRate(geneIndex, synthesisRateCategory, false);
    double phiValue_proposed = parameter->getSynthesisRate(geneIndex, synthesisRateCategory, true);
    double logPhi = std::log(phiValue); 
    double logPhi_proposed = std::log(phiValue_proposed);

    for (unsigned index = 0; index < positions.size();index++)
    {
        codon = gene.geneData.indexToCodon(positions[index]);
        currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
        currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, false);
        currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);
        sigma[index+1] = elongationUntilIndexApproximation2ProbabilityLog(currAlpha, currLambdaPrime/U,1/currNSERate,sigma[index]);
        alpha[index] = currAlpha;
        lambda_prime[index] = currLambdaPrime;
    }
    // double log_final_sigma = sigma[sigma.size()-1];
    // double final_sigma = std::exp(log_final_sigma);
    // double initiation_rate = phiValue/final_sigma;
    // double initiation_rate_proposed = phiValue_proposed/final_sigma;
//Loop through the positions twice, but the first time, which is relatively cheap operations,
//sets it up so we can do the expensive operations in parallel
#ifdef _OPENMP
    //#ifndef __APPLE__
#pragma omp parallel for reduction(+:logLikelihood,logLikelihood_proposed)
#endif
    for (unsigned index = 0; index < positions.size(); index++)
    {
      
        unsigned codonIndex = positions[index];
        unsigned currRFPObserved = rfpCounts[index];
        double currSigma = std::exp(sigma[index]);
        if (local_lgamma_currentAlpha[codonIndex] < -5)
        {
            local_lgamma_currentAlpha[codonIndex] = std::lgamma(alpha[index]);
        }
        if (local_log_currentLambdaPrime[codonIndex] > 500)
        {
            local_log_currentLambdaPrime[codonIndex] = std::log(lambda_prime[index]);
        } 
        double currLgammaRFPAlpha = std::lgamma(alpha[index]+currRFPObserved);
        logLikelihood += calculateLogLikelihoodPerCodonPerGene(alpha[index], lambda_prime[index], currRFPObserved, phiValue, currSigma, local_lgamma_currentAlpha[codonIndex],
                                local_log_currentLambdaPrime[codonIndex], logPhi, currLgammaRFPAlpha);
        logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(alpha[index], lambda_prime[index], currRFPObserved, phiValue_proposed, currSigma, local_lgamma_currentAlpha[codonIndex],
                                local_log_currentLambdaPrime[codonIndex], logPhi_proposed, currLgammaRFPAlpha);
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
}



void PANSEModel::fillMatrices(Genome& genome,bool init_sigma_vector)
{
    unsigned n = getNumMixtureElements();
    //These vectors should be visible to all threads. If one thread tries to overwrite another one, should only result in it replacing the same value
    for (unsigned j = 0; j < n; j++)
    {
        std::vector<double> tmp = std::vector<double> (getGroupListSize(), -1000);
        std::vector<double> tmp_2 = std::vector<double> (getGroupListSize(),1000);
        lgamma_currentAlpha.push_back(tmp);
        log_currentLambdaPrime.push_back(tmp_2);
    }
    lgamma_rfp_alpha.resize(50);
    for (unsigned k = 0; k < 50;k++)
    {
        lgamma_rfp_alpha[k].resize(n);
        for (unsigned j = 0; j < n; j++)
        {
            std::vector<double> tmp = std::vector<double> (getGroupListSize(), -1000);
            lgamma_rfp_alpha[k][j]=tmp;
        }
    }
    if (init_sigma_vector)
    {
        unsigned numGenes = genome.getGenomeSize();
        current_sigma.resize(numGenes);
        for (unsigned k = 0; k < numGenes;k++)
        {
            Gene *gene = &genome.getGene(k);
            std::vector <unsigned> positions = gene->geneData.getPositionCodonID();
            // sigma represents probability, should be between 0 and 1 on non-log scale, so will range from negative infinity to 0 on log scale
            std::vector<double> tmp = std::vector<double> (positions.size()+1, 100);
            current_sigma[k] = tmp;
            current_sigma[k][0] = 0;
        }
    }

}

void PANSEModel::clearMatrices()
{
    lgamma_currentAlpha.clear();
    log_currentLambdaPrime.clear();
    lgamma_rfp_alpha.clear();
    current_sigma.clear();
}


void PANSEModel::calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome, std::vector<double> &logAcceptanceRatioForAllMixtures)
{
    std::vector<std::string> groups = parameter -> getGroupList();
    if (groups[0] == grouping)
    {
        fillMatrices(genome);
    }
    double logLikelihood = 0.0;
    double logLikelihood_proposed = 0.0;
    double logLikelihood_adjusted = 0.0;
    double logLikelihood_proposed_adjusted = 0.0;
    double propAlpha, propLambdaPrime, propNSERate;
    double currAlpha, currLambdaPrime, currNSERate;
    double currAdjustmentTerm = 0;
    double propAdjustmentTerm = 0;
    Gene *gene;
    unsigned index = SequenceSummary::codonToIndex(grouping);
    std::vector<double> codonSigmas;
    unsigned n = getNumMixtureElements();
    unsigned Y = genome.getSumRFP();
#ifdef _OPENMP
    //#ifndef __APPLE__
#pragma omp parallel for private(gene,currAlpha,currLambdaPrime,currNSERate,propAlpha,propLambdaPrime,propNSERate) reduction(+:logLikelihood,logLikelihood_proposed)
#endif
    for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
    {
        gene = &genome.getGene(i);
        unsigned currNumCodonsInMRNA = gene->geneData.getCodonCountForCodon(grouping);
        if (currNumCodonsInMRNA == 0) continue;
        unsigned positionalRFPCount;
        unsigned codonIndex;
        std::string codon;
        double currLgammaRFPAlpha;
        
        unsigned mixtureElement = parameter->getMixtureAssignment(i);
        double U = getPartitionFunction(mixtureElement, false)/Y;
        // how is the mixture element defined. Which categories make it up
        unsigned alphaCategory = parameter->getMutationCategory(mixtureElement);
        unsigned lambdaPrimeCategory = parameter->getSelectionCategory(mixtureElement);
        unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(mixtureElement);
        

        double phiValue = parameter->getSynthesisRate(i, synthesisRateCategory, false);
        double log_phi = std::log(phiValue);
        std::vector <unsigned> positions = gene->geneData.getPositionCodonID();
        std::vector <unsigned> rfpCounts = gene->geneData.getRFPCount(/*RFPCountColumn*/ 0);
        
        // std::vector<double> alpha(positions.size(),0);
        // std::vector<double> lambda_prime(positions.size(),0); 
        // std::vector<double> proposed_sigma(positions.size()+1,0);
        // for (unsigned positionIndex = 0; positionIndex < positions.size();positionIndex++)
        // {
        //     codon = gene->geneData.indexToCodon(positions[positionIndex]);
        //     codonIndex = positions[positionIndex];
        //     currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
        //     if (lgamma_currentAlpha[alphaCategory][codonIndex] < -5)
        //     {
        //         lgamma_currentAlpha[alphaCategory][codonIndex] = std::lgamma(currAlpha);
        //     }
        //     currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, false);
        //     if (log_currentLambdaPrime[lambdaPrimeCategory][codonIndex] > 500)
        //     {
        //         log_currentLambdaPrime[lambdaPrimeCategory][codonIndex] = std::log(currLambdaPrime);
        //     }
        //     currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);
        //     if (current_sigma[i][positionIndex+1] > 0)
        //     {
        //         current_sigma[i][positionIndex+1] = elongationUntilIndexApproximation2ProbabilityLog(currAlpha, currLambdaPrime/U,1/currNSERate,current_sigma[i][positionIndex]);
        //     }
        //     if(codon == grouping)
        //     {
        //         propAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, true);
        //         propLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, true);
        //         propNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, true);
        //         proposed_sigma[positionIndex+1]= elongationUntilIndexApproximation2ProbabilityLog(propAlpha, propLambdaPrime/U,1/propNSERate,proposed_sigma[positionIndex]);
        //     }
        //     else
        //     {
        //         proposed_sigma[positionIndex+1] = elongationUntilIndexApproximation2ProbabilityLog(currAlpha, currLambdaPrime/U,1/currNSERate,proposed_sigma[positionIndex]);
        //     }
        //     alpha[positionIndex] = currAlpha;
        //     lambda_prime[positionIndex] = currLambdaPrime;
        // }
        
        // double log_final_sigma = current_sigma[i][positions.size()];
        // double log_final_sigma_proposed = proposed_sigma[proposed_sigma.size()-1];
        // double final_sigma = std::exp(log_final_sigma);
        // double final_sigma_proposed = std::exp(log_final_sigma_proposed);
        // double initiation_rate = phiValue/final_sigma;
        // double initiation_rate_proposed = phiValue/final_sigma_proposed;
        double propSigma = 0;
        for (unsigned positionIndex = 0; positionIndex < positions.size(); positionIndex++)
        {
            positionalRFPCount = rfpCounts[positionIndex];
            codonIndex = positions[positionIndex];
            codon = gene->geneData.indexToCodon(codonIndex);
            currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
            if (lgamma_currentAlpha[alphaCategory][codonIndex] < -5)
            {
                lgamma_currentAlpha[alphaCategory][codonIndex] = std::lgamma(currAlpha);
            }
            currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, false);
            if (log_currentLambdaPrime[lambdaPrimeCategory][codonIndex] > 500)
            {
                log_currentLambdaPrime[lambdaPrimeCategory][codonIndex] = std::log(currLambdaPrime);
            }
            currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);
            //Check if lgamma(alpha + rfp_count) already counted for this codon 
            if (positionalRFPCount < 150)
            {
                if (lgamma_rfp_alpha[positionalRFPCount][alphaCategory][codonIndex] < -5)
                {
                    lgamma_rfp_alpha[positionalRFPCount][alphaCategory][codonIndex] = std::lgamma(currAlpha + positionalRFPCount);
                }
                currLgammaRFPAlpha = lgamma_rfp_alpha[positionalRFPCount][alphaCategory][codonIndex];
            }
            else
            {
                currLgammaRFPAlpha = std::lgamma(currAlpha + positionalRFPCount);
            }

            if(codon == grouping)
            {
                propAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, true);
                propLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, true);
                propNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, true);
                logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(propAlpha, propLambdaPrime, positionalRFPCount,
                                    phiValue,std::exp(propSigma),std::lgamma(propAlpha),std::log(propLambdaPrime),log_phi,std::lgamma(propAlpha+positionalRFPCount));
                propSigma = elongationUntilIndexApproximation2ProbabilityLog(propAlpha, propLambdaPrime/U,1/propNSERate,propSigma);
            }
            else
            {
                logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, positionalRFPCount,
                                   phiValue,std::exp(propSigma),lgamma_currentAlpha[alphaCategory][codonIndex],log_currentLambdaPrime[lambdaPrimeCategory][codonIndex],log_phi,currLgammaRFPAlpha);
                propSigma = elongationUntilIndexApproximation2ProbabilityLog(currAlpha, currLambdaPrime/U,1/currNSERate,propSigma);
            }

            logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, positionalRFPCount,
                                   phiValue, std::exp(current_sigma[i][positionIndex]),lgamma_currentAlpha[alphaCategory][codonIndex],log_currentLambdaPrime[lambdaPrimeCategory][codonIndex],log_phi,currLgammaRFPAlpha);  
        
            if (current_sigma[i][positionIndex+1] > 0)
            {
                current_sigma[i][positionIndex+1] = elongationUntilIndexApproximation2ProbabilityLog(currAlpha, currLambdaPrime/U,1/currNSERate,current_sigma[i][positionIndex]);
            }

        }
        //my_print("Done with gene %\n",index);
    }
    for (unsigned j = 0; j < n; j++)
    {
        unsigned alphaCategory = parameter->getMutationCategory(j);
        unsigned lambdaPrimeCategory = parameter->getSelectionCategory(j);
        currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, grouping, false);
        currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, grouping, false);
        currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, grouping, false);
        propAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, grouping, true);
        propLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, grouping, true);
        propNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, grouping, true);
        currAdjustmentTerm += std::log(currAlpha) + std::log(currLambdaPrime) + std::log(currNSERate);
        propAdjustmentTerm += std::log(propAlpha) + std::log(propLambdaPrime) + std::log(propNSERate);
    }
  
    logAcceptanceRatioForAllMixtures[0] = logLikelihood_proposed - logLikelihood - (currAdjustmentTerm - propAdjustmentTerm);
	logAcceptanceRatioForAllMixtures[1] = logLikelihood - propAdjustmentTerm;
	logAcceptanceRatioForAllMixtures[2] = logLikelihood_proposed - currAdjustmentTerm;
	logAcceptanceRatioForAllMixtures[3] = logLikelihood;
	logAcceptanceRatioForAllMixtures[4] = logLikelihood_proposed;

    if (groups[getGroupListSize()-1] == grouping)
    {
        clearMatrices();
    }
	//If this is the last codon continue with updating all accepted CSPs
    //Alex: This won't work because doesn't add last codon, adds it on next round
	// if (SequenceSummary::codonToIndex(grouping) == (getGroupListSize() - 1))
	// {
	//       parameter->completeUpdateCodonSpecificParameter();
	// }
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


    if (withPhi)
    {
        // one for each noiseOffset, and one for stdDevSynthesisRate
        logProbabilityRatio.resize(getNumPhiGroupings() + 2);
    }
    else
        logProbabilityRatio.resize(2);



#ifdef _OPENMP
    //#ifndef __APPLE__
#pragma omp parallel for reduction(+:lpr)
#endif
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
    double currAlpha, currLambdaPrime, currNSERate;
    double logLikelihood, logLikelihood_proposed;
    unsigned n = getNumMixtureElements();
    lpr = 0.0;

    fillMatrices(genome,false);
    unsigned Y = genome.getSumRFP();
#ifdef _OPENMP
//#ifndef __APPLE__
#pragma omp parallel for private(gene,currAlpha,currLambdaPrime,currNSERate) reduction(+:logLikelihood,logLikelihood_proposed)
#endif
    for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
    {
       
        unsigned positionalRFPCount;
        unsigned codonIndex;
        std::string codon;
        double currLgammaRFPAlpha;
        
        gene = &genome.getGene(i);
        unsigned mixtureElement = parameter->getMixtureAssignment(i);

        // how is the mixture element defined. Which categories make it up
        double currU = getPartitionFunction(mixtureElement, false)/Y;
        double propU = getPartitionFunction(mixtureElement, true)/Y;
        unsigned alphaCategory = parameter->getMutationCategory(mixtureElement);
        unsigned lambdaPrimeCategory = parameter->getSelectionCategory(mixtureElement);
        unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(mixtureElement);
        // get non codon specific values, calculate likelihood conditional on these
        double phiValue = parameter->getSynthesisRate(i, synthesisRateCategory, false);
        double log_phi = std::log(phiValue);
        
        std::vector <unsigned> positions = gene->geneData.getPositionCodonID();
        std::vector <unsigned> rfpCounts = gene->geneData.getRFPCount(/*RFPCountColumn*/ 0);

        // std::vector<double> alpha(positions.size(),0);
        // std::vector<double> lambda_prime(positions.size(),0);
        // std::vector<double> nse_rate(positions.size(),0); 

        // std::vector<double> proposed_sigma(positions.size()+1,0);
        // for (unsigned positionIndex = 0; positionIndex < positions.size();positionIndex++)
        // {
        //     codon = gene->geneData.indexToCodon(positions[positionIndex]);
        //     codonIndex = positions[positionIndex];
        //     currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
        //     if (lgamma_currentAlpha[alphaCategory][codonIndex] < -5)
        //     {
        //         lgamma_currentAlpha[alphaCategory][codonIndex] = std::lgamma(currAlpha);
        //     }
        //     currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, false);
        //     if (log_currentLambdaPrime[lambdaPrimeCategory][codonIndex] > 500)
        //     {
        //         log_currentLambdaPrime[lambdaPrimeCategory][codonIndex] = std::log(currLambdaPrime);
        //     }
        //     currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);
        //     current_sigma[i][positionIndex+1] = elongationUntilIndexApproximation2ProbabilityLog(currAlpha, currLambdaPrime/currU,1/currNSERate,current_sigma[i][positionIndex]);
        //     proposed_sigma[positionIndex+1] = elongationUntilIndexApproximation2ProbabilityLog(currAlpha, currLambdaPrime/propU,1/currNSERate,proposed_sigma[positionIndex]);
        //     alpha[positionIndex] = currAlpha;
        //     lambda_prime[positionIndex] = currLambdaPrime;
        //     nse_rate[positionIndex] = currNSERate; 
        // }

        // double log_final_sigma = current_sigma[i][positions.size()];
        // double log_final_sigma_proposed = proposed_sigma[proposed_sigma.size()-1];
        // double final_sigma = std::exp(log_final_sigma);
        // double final_sigma_proposed = std::exp(log_final_sigma_proposed);
        // double initiation_rate = phiValue/final_sigma;
        // double initiation_rate_proposed = phiValue/final_sigma_proposed;

        double currSigma = 0;
        double propSigma = 0;

        for (unsigned positionIndex = 0; positionIndex < positions.size(); positionIndex++)
        {
            positionalRFPCount = rfpCounts[positionIndex];
            codonIndex = positions[positionIndex];
            codon = gene->geneData.indexToCodon(codonIndex);
            currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
            if (lgamma_currentAlpha[alphaCategory][codonIndex] < -5)
            {
                lgamma_currentAlpha[alphaCategory][codonIndex] = std::lgamma(currAlpha);
            }
            currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PANSEParameter::lmPri, codon, false);
            if (log_currentLambdaPrime[lambdaPrimeCategory][codonIndex] > 500)
            {
                log_currentLambdaPrime[lambdaPrimeCategory][codonIndex] = std::log(currLambdaPrime);
            }
            currNSERate = getParameterForCategory(alphaCategory, PANSEParameter::nse, codon, false);           //Check if lgamma(alpha + rfp_count) already counted for this codon 
            if (positionalRFPCount < 50)
            {
                if (lgamma_rfp_alpha[positionalRFPCount][alphaCategory][codonIndex] < -5)
                {
                    lgamma_rfp_alpha[positionalRFPCount][alphaCategory][codonIndex] = std::lgamma(currAlpha + positionalRFPCount);
                }
                currLgammaRFPAlpha = lgamma_rfp_alpha[positionalRFPCount][alphaCategory][codonIndex];
            }
            else
            {
                currLgammaRFPAlpha = std::lgamma(currAlpha + positionalRFPCount);
            }
            logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, positionalRFPCount,
                                   phiValue,std::exp(propSigma),lgamma_currentAlpha[alphaCategory][codonIndex],log_currentLambdaPrime[lambdaPrimeCategory][codonIndex],log_phi,currLgammaRFPAlpha);
            logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, positionalRFPCount,
                                   phiValue, std::exp(currSigma),lgamma_currentAlpha[alphaCategory][codonIndex],log_currentLambdaPrime[lambdaPrimeCategory][codonIndex],log_phi,currLgammaRFPAlpha);   
       
           currSigma = elongationUntilIndexApproximation2ProbabilityLog(currAlpha, currLambdaPrime/currU,1/currNSERate,currSigma);
           propSigma = elongationUntilIndexApproximation2ProbabilityLog(currAlpha, currLambdaPrime/propU,1/currNSERate,propSigma);
            
        }
          
    }
    for (unsigned j = 0; j < n; j++)
    {
        lpr -= (std::log(getPartitionFunction(j, false)) - std::log(getPartitionFunction(j, true)));
    }
    lpr += logLikelihood_proposed - logLikelihood; 
    logProbabilityRatio[1] = lpr;
    clearMatrices();
    if (withPhi)
    {
        for (unsigned i = 0; i < parameter->getNumObservedPhiSets(); i++)
        {
            
            lpr = 0.0;
            double noiseOffset = getNoiseOffset(i, false);
            double noiseOffset_proposed = getNoiseOffset(i, true);
            double observedSynthesisNoise = getObservedSynthesisNoise(i);
#ifdef _OPENMP
//#ifndef __APPLE__
#pragma omp parallel for reduction(+:lpr)
#endif
            for (unsigned j = 0u; j < genome.getGenomeSize(); j++)
            {
                unsigned mixtureAssignment = getMixtureAssignment(j);
                mixtureAssignment = getSynthesisRateCategory(mixtureAssignment);
                double logPhi = std::log(getSynthesisRate(j, mixtureAssignment, false));
                double obsPhi = genome.getGene(j).getObservedSynthesisRate(i);
                if (obsPhi > -1.0)
                {
                    double logObsPhi = std::log(obsPhi);
                    double proposed = Parameter::densityNorm(logObsPhi, logPhi + noiseOffset_proposed, observedSynthesisNoise, true);
                    double current = Parameter::densityNorm(logObsPhi, logPhi + noiseOffset, observedSynthesisNoise, true);
                    lpr += proposed - current;
                }
            }
            logProbabilityRatio[i+2] = lpr;
        }
    }
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
    if (withPhi)
    {
        updateNoiseOffsetTrace(sample);
        updateObservedSynthesisNoiseTrace(sample);
    }
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
    if (withPhi)
        adaptNoiseOffsetProposalWidth(adaptiveWidth, adapt);
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
    if (withPhi)
    {
        parameter->proposeNoiseOffset();
    }
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


void PANSEModel::completeUpdateCodonSpecificParameter()
{
    parameter->completeUpdateCodonSpecificParameter();
}

//Noise offset functions

double PANSEModel::getNoiseOffset(unsigned index, bool proposed)
{
    return parameter->getNoiseOffset(index, proposed);
}


double PANSEModel::getObservedSynthesisNoise(unsigned index)
{
    return parameter->getObservedSynthesisNoise(index);
}


double PANSEModel::getCurrentNoiseOffsetProposalWidth(unsigned index)
{
    return parameter->getCurrentNoiseOffsetProposalWidth(index);
}


void PANSEModel::updateNoiseOffset(unsigned index)
{
    parameter->updateNoiseOffset(index);
}


void PANSEModel::updateNoiseOffsetTrace(unsigned sample)
{
    parameter->updateNoiseOffsetTraces(sample);
}


void PANSEModel::updateObservedSynthesisNoiseTrace(unsigned sample)
{
    parameter->updateObservedSynthesisNoiseTraces(sample);
}


void PANSEModel::adaptNoiseOffsetProposalWidth(unsigned adaptiveWidth, bool adapt)
{
    parameter->adaptNoiseOffsetProposalWidth(adaptiveWidth, adapt);
}



void PANSEModel::updateGibbsSampledHyperParameters(Genome &genome)
{
  // estimate s_epsilon by sampling from a gamma distribution and transforming it into an inverse gamma sample
    
    if (withPhi)
    {
        if(!fix_sEpsilon)
        {
            double shape = ((double)genome.getGenomeSize() - 1.0) / 2.0;
            for (unsigned i = 0; i < parameter->getNumObservedPhiSets(); i++)
            {
                double rate = 0.0; //Prior on s_epsilon goes here?
                unsigned mixtureAssignment;
                double noiseOffset = getNoiseOffset(i);
                for (unsigned j = 0; j < genome.getGenomeSize(); j++)
                {
                    mixtureAssignment = parameter->getMixtureAssignment(j);
                    double obsPhi = genome.getGene(j).getObservedSynthesisRate(i);
                    if (obsPhi > -1.0)
                    {
                        double sum = std::log(obsPhi) - noiseOffset - std::log(parameter->getSynthesisRate(j, mixtureAssignment, false));
                        rate += (sum * sum);
                    }else{
                        // missing observation.
                        shape -= 0.5;
                        //Reduce shape because initial estimate assumes there are no missing observations
                    }
                }
                rate /= 2.0;
                double rand = parameter->randGamma(shape, rate);

                // Below the gamma sample is transformed into an inverse gamma sample
                // According to Gilchrist et al (2015) Supporting Materials p. S6
                // The sample 1/T is supposed to be equal to $s_\epsilon^2$.
                double sepsilon = std::sqrt(1.0/rand);
                parameter->setObservedSynthesisNoise(i, sepsilon);
            }
        }
    }
}


void PANSEModel::updateAllHyperParameter()
{
    updateStdDevSynthesisRate();
    updatePartitionFunction();
    if (withPhi)
    {
        for (unsigned i =0; i < parameter->getNumObservedPhiSets(); i++)
        {
            updateNoiseOffset(i);
        }
    }
}


void PANSEModel::updateHyperParameter(unsigned hp)
{
    // NOTE: when adding additional hyper parameter, also add to updateAllHyperParameter()
    if (hp == 0)
    {
        updateStdDevSynthesisRate();
    }
    else if (hp == 1)
    {
        updatePartitionFunction();
    }       
    else if (hp > 1 and withPhi)
    {   
        //subtract off 2 because the first two parameters withh be the updateStdDevSynthesisRate
        updateNoiseOffset(hp - 2);
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
        std::vector <unsigned> rfpCount;
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
            xx = rgamma(1, alpha,1/lambdaPrime);
            sigma *= (v/(xx[0] + v));
            xx = rpois(1, xx[0] * phi * sigma);
            rfpCount.push_back(xx[0]);
#else
            std::gamma_distribution<double> GDistribution(alpha, 1.0/lambdaPrime);
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

    rv = 10000.0 / x;
    for(i = 10000; i > 0; i--){
        if(i % 2 == 0) rv = (double) ((int) (i / 2)) / (x + rv);
        if(i % 2 != 0) rv = ((double) (((int) (i / 2) + 1) - s)) / (1 + rv);
    }

    return x + rv;
}

//Upper incomplete gamma function
double PANSEModel::UpperIncompleteGamma(double s, double x)
{
    double d, rv;

    rv = pow(x, s) * exp(0 - x);

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
double PANSEModel::elongationProbability(double currAlpha, double currLambda, double currNSE){
    return std::pow(currLambda * currNSE, currAlpha) * std::exp(currLambda * currNSE) * UpperIncompleteGamma(currAlpha, currLambda * currNSE);
}

//Log probability of elongation at current codon
double PANSEModel::elongationProbabilityLog(double currAlpha, double currLambda, double currNSE){
    return (currLambda * currNSE) + currAlpha * (std::log(currLambda) + std::log(currNSE)) + UpperIncompleteGammaLog(1- currAlpha, currLambda * currNSE);
}

double PANSEModel::elongationUntilIndexProbability(int index, std::vector <double> lambdas, std::vector <double> NSERates){
    return 0;
}

double PANSEModel::elongationUntilIndexProbabilityLog(int index, std::vector <double> lambdas, std::vector <double> NSERates){
    return 0;
}


double PANSEModel::elongationUntilIndexApproximation1Probability(double alpha, double lambda, double v, double current)
{
	//my_print("Alpha % Lambda % v %\n",alpha,lambda,v);
//    if (proposed)
//    {
//        propSigmaCalculationSummationFor1 += (alpha/(lambda * v));
//        return 1 - propSigmaCalculationSummationFor1;
//    }
//    else
//    {
//        currSigmaCalculationSummationFor1 += (alpha/(lambda * v));
//        return 1 - currSigmaCalculationSummationFor1;
//    }
	return (current+(alpha/(lambda * v)));
}

double PANSEModel::elongationUntilIndexApproximation2Probability(double alpha, double lambda, double v, bool proposed)
{
    if (proposed)
    {
        propSigmaCalculationSummationFor2 += (alpha/(lambda * v)) *
                (-1*elongationUntilIndexApproximation1Probability(alpha, lambda, v, proposed))
                + (alpha/(lambda * lambda * v * v));
        return 1 + propSigmaCalculationSummationFor2;
    }
    else
    {
        currSigmaCalculationSummationFor2 += (alpha/(lambda * v)) *
                (-1*elongationUntilIndexApproximation1Probability(alpha, lambda, v, proposed))
                + (alpha/(lambda * lambda * v * v));
        return 1 + currSigmaCalculationSummationFor2;
    }
}

double PANSEModel::elongationUntilIndexApproximation1ProbabilityLog(double alpha, double lambda, double v, bool proposed)
{
    if (proposed)
    {
        propSigmaCalculationSummationFor1 -= (alpha/(lambda * v));
        return propSigmaCalculationSummationFor1;
    }
    else
    {
    	currSigmaCalculationSummationFor1 -= (alpha/(lambda * v));
    	return currSigmaCalculationSummationFor1;
    }
}
double PANSEModel::elongationUntilIndexApproximation2ProbabilityLog(double alpha, double lambda, double v, double current)
{
//    if (proposed)
//       {
//        propSigmaCalculationSummationFor2 += -(alpha/(lambda * v)) + (alpha/(lambda * lambda * v * v))
//                   + (alpha/(lambda * v)) * (alpha/(lambda * v)) / 2;
//           return propSigmaCalculationSummationFor2;
//       }
//       currSigmaCalculationSummationFor2 += -(alpha/(lambda * v)) + (alpha/(lambda * lambda * v * v))
//                           + (alpha/(lambda * v)) * (alpha/(lambda * v)) / 2;
//       return currSigmaCalculationSummationFor2;
		return current + (-(alpha/(lambda * v)) + (alpha/(lambda * lambda * v * v))
                          + ((alpha/(lambda * v)) * (alpha/(lambda * v))) / 2);

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

    term1 = std::lgamma(curralpha + sample_size) - lgamma(curralpha);
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
