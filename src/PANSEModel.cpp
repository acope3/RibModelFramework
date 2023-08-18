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
    parameter_types = {"Elongation","NSERate"};
    //ctor
}


PANSEModel::~PANSEModel()
{
    //dtor
    //TODO: call Parent's deconstructor
    //delete parameter;
}

std::vector<std::vector<unsigned>> PANSEModel::getElongationMixtureCategories()
{
  unsigned alphaCategory, lambdaCategory, nseCategory;
  std::vector<std::vector<unsigned>> mixture_to_category; 
  unsigned n = getNumElongationMixtureElements(); //Want to use elongation mixtures, not phi mixtures (the number of mixtures specificed by numMixtures)
  mixture_to_category.resize(n);
  for (unsigned i = 0; i < n; i++)
  {
    mixture_to_category[i].resize(3); //for alpha, lambda, and NSE categories
    alphaCategory = parameter->getMutationCategory(i);
    lambdaCategory = parameter->getSelectionCategory(i);
    nseCategory = parameter->getNSECategory(i);
    mixture_to_category[i][0] = alphaCategory;
    mixture_to_category[i][1] = lambdaCategory;
    mixture_to_category[i][2] = nseCategory;
  }
  return(mixture_to_category);
}


void PANSEModel::fillMatrices(Genome& genome)
{
    std::string codon;
    unsigned alphaCategory, lambdaCategory, nseCategory;
    double currAlpha,currLambda,currNSERate;
    double U;
    unsigned n = getNumElongationMixtureElements(); //Want to use elongation mixtures, not phi mixtures (the number of mixtures specificed by numMixtures)
    unsigned mixtures = getNumMixtureElements();
    unsigned num_codons = getGroupListSize();

    //These vectors should be visible to all threads. If one thread tries to overwrite another one, should only result in it replacing the same value
    unsigned long Y = parameter->getTotalRFPCount();
    lgamma_rfp_alpha.resize(50);
    for (unsigned i=0; i < 50; i++)
    {
        lgamma_rfp_alpha[i].resize(n);
        for (unsigned j = 0; j < n; j ++)
        {
            lgamma_rfp_alpha[i][j].resize(getGroupListSize());
        }
    }
    //TO DO: Make prob_successful matrix instead of vector
    prob_successful.resize(n);
    for (unsigned j = 0; j < n; j++)
    {
        prob_successful[j].resize(getGroupListSize(),0);

        alphaCategory = parameter->getMutationCategory(j);
        lambdaCategory = parameter->getSelectionCategory(j);
        nseCategory = parameter->getNSECategory(j);

        std::vector<double> tmp = std::vector<double> (getGroupListSize(), 0);

        for (unsigned k = 0; k < num_codons; k++)
        {
            codon = getGrouping(k);
            currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
            currLambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, codon, false);
            currNSERate = getParameterForCategory(nseCategory, PANSEParameter::nse, codon, false);

            tmp[k] = std::lgamma(currAlpha);
            prob_successful[j][k] = elongationUntilIndexApproximation2ProbabilityLog(currAlpha, currLambda, 1/currNSERate);

            if (prob_successful[j][k] > 0.0)
            {
            	//prob_successful[j][k] = std::numeric_limits<double>::quiet_NaN();
              prob_successful[j][k] = 0.0;
            }
            for (unsigned l=0; l < 50; l++)
            {
                // l: RFP count at position, j is mixture, k is codon
                lgamma_rfp_alpha[l][j][k] = std::lgamma(currAlpha + l);
            }
        }
        lgamma_currentAlpha.push_back(tmp);
        //log_currentLambda.push_back(tmp_2);
    }
    log_currentLambda.resize(mixtures);
    for (unsigned i = 0; i < mixtures; i++)
    {
    	U = getPartitionFunction(i, false)/Y;
    	std::vector<std::vector<double>> tmp = std::vector<std::vector<double>>(n);
    	for (unsigned j = 0; j < n; j++)
    	{
    		std::vector<double> tmp_2 = std::vector<double> (num_codons,0);
            lambdaCategory = parameter->getSelectionCategory(j);
            for (unsigned k = 0; k < num_codons; k++)
            {
            	codon = getGrouping(k);
            	currLambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, codon, false);
            	tmp_2[k] = std::log(currLambda) + std::log(U);
            }
            tmp[j] = tmp_2;
    	}
    	log_currentLambda[i] = tmp;
    }
    mixture_to_category = getElongationMixtureCategories();
}

void PANSEModel::clearMatrices()
{
    lgamma_currentAlpha.clear();
    log_currentLambda.clear();
    lgamma_rfp_alpha.clear();
    prob_successful.clear();
    mixture_to_category.clear();
}


double PANSEModel::calculateLogLikelihoodPerCodonPerGeneByPosition(double currAlpha, double currLambdaPrime,
                                                                unsigned currRFPObserved, double phiValue, double prevSigma)
{
  
  double term1 = std::lgamma(currAlpha + currRFPObserved) - std::lgamma(currAlpha);
  double term2 = std::log(phiValue) + std::log(prevSigma) - std::log(currLambdaPrime + phiValue * prevSigma);
  double term3 = std::log(currLambdaPrime) - std::log(currLambdaPrime + phiValue * prevSigma);
  
  term2 *= currRFPObserved;
  term3 *= currAlpha;
  
  double rv = term1 + term2 + term3;
  return rv;
}



double PANSEModel::calculateLogLikelihoodPerCodonPerGene(double currAlpha, double currLambdaPrime,
        unsigned currRFPObserved, double phiValue, double prevSigma, double lgamma_currAlpha, double log_currLambdaPrime, double log_phi,double lgamma_rfp_alpha)
{

    double term1 = lgamma_rfp_alpha - (lgamma_currAlpha);//std::lgamma(currAlpha);
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
    double currAlpha,currLambda,currNSERate;
    std::string codon;
    double currLgammaRFPAlpha;
    double logLikelihood = 0.0;
    double logLikelihood_proposed = 0.0;
    unsigned codonIndex;
    unsigned long positionalRFPCount;
    unsigned alphaCategory, lambdaCategory, nseCategory;
    int codonMixture, codonMixture_w_flag;

    std::vector<int> positionMixture = gene.geneData.getPositionMixture();
    std::vector <unsigned> positions = gene.geneData.getPositionCodonID();
    std::vector <unsigned long> rfpCounts = gene.geneData.getRFPCount(0);


    unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(k);

    unsigned mixtureElement = parameter->getMixtureAssignment(geneIndex);
    
    unsigned long Y = parameter->getTotalRFPCount();
    double U = getPartitionFunction(k, false)/Y;
    

    double phiValue = parameter->getSynthesisRate(geneIndex, synthesisRateCategory, false);
    double phiValue_proposed = parameter->getSynthesisRate(geneIndex, synthesisRateCategory, true);
    

    double logPhi = std::log(phiValue); 
    double logPhi_proposed = std::log(phiValue_proposed);
    
    double currSigma = 0;
    
    
    for (unsigned positionIndex = 0; positionIndex < positions.size();positionIndex++)
    {
        codonMixture_w_flag = positionMixture[positionIndex] + 1; // put back on 1-indexed scale to check if original value was negative or not
        if (codonMixture_w_flag < 0)
        {
          codonMixture = -1 * (codonMixture_w_flag) - 1; //if negative, get codonMixture if were not ignoring
        }
        else if (codonMixture_w_flag > 0)
        {
          codonMixture = codonMixture_w_flag - 1; //if positive, get codonMixture
        } else{
          my_print("ERROR: Your column indicating elongation mixtures contains 0. Should be a non-zero positive (include position in likelihood) or negative (use only for sigma) number. Exiting program.\n");
          exit(1);
        }
        codonIndex = positions[positionIndex];
        if (codonMixture_w_flag > 0)
        {
          positionalRFPCount = rfpCounts[positionIndex];
          codon = gene.geneData.indexToCodon(codonIndex);
          alphaCategory = mixture_to_category[codonMixture][0];
          lambdaCategory = mixture_to_category[codonMixture][1];
          nseCategory = mixture_to_category[codonMixture][2];
  
  
          currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
          currLambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, codon, false);
          currNSERate = getParameterForCategory(nseCategory, PANSEParameter::nse, codon, false);
  
          if (positionalRFPCount < 50)
          {
              
              currLgammaRFPAlpha = lgamma_rfp_alpha[positionalRFPCount][alphaCategory][codonIndex];
          }
          else
          {
              currLgammaRFPAlpha = std::lgamma(currAlpha + positionalRFPCount);
          }
          logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambda * U, positionalRFPCount, phiValue, std::exp(currSigma), 
                                  lgamma_currentAlpha[alphaCategory][codonIndex],log_currentLambda[synthesisRateCategory][lambdaCategory][codonIndex], logPhi, currLgammaRFPAlpha);
          logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambda * U, positionalRFPCount, phiValue_proposed, std::exp(currSigma), 
                                  lgamma_currentAlpha[alphaCategory][codonIndex],log_currentLambda[synthesisRateCategory][lambdaCategory][codonIndex], logPhi_proposed, currLgammaRFPAlpha);
          
        }
        currSigma = currSigma + prob_successful[codonMixture][codonIndex];
    }


    double stdDevSynthesisRate = parameter->getStdDevSynthesisRate(synthesisRateCategory, false);
    double logPhiProbability = Parameter::densityLogNorm(phiValue, (-(stdDevSynthesisRate * stdDevSynthesisRate) * 0.5), stdDevSynthesisRate, true);
    double logPhiProbability_proposed = Parameter::densityLogNorm(phiValue_proposed, (-(stdDevSynthesisRate * stdDevSynthesisRate) * 0.5), stdDevSynthesisRate, true);
	
    if (withPhi)
    {
        for (unsigned i = 0; i < parameter->getNumObservedPhiSets(); i++)
        {
            double obsPhi = gene.getObservedSynthesisRate(i);
            if (obsPhi > -1.0)
            {
                double logObsPhi = std::log(obsPhi);
                logPhiProbability += Parameter::densityNorm(logObsPhi, std::log(phiValue) + getNoiseOffset(i), getObservedSynthesisNoise(i), true);
                logPhiProbability_proposed += Parameter::densityNorm(logObsPhi, std::log(phiValue_proposed) + getNoiseOffset(i), getObservedSynthesisNoise(i), true);
            }
        }
    }


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





void PANSEModel::calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome, std::vector<double> &logAcceptanceRatioForAllMixtures,std::string param)
{
    std::vector<std::string> groups = parameter -> getGroupList();
    Gene *gene;
    double propAlpha, propLambda, propNSERate;
    double currAlpha, currLambda, currNSERate;
    unsigned alphaCategory, lambdaCategory, nseCategory;

    double logLikelihood, logPosterior = 0.0;
    double logLikelihood_proposed, logPosterior_proposed = 0.0;
    double currAdjustmentTerm = 0;
    double propAdjustmentTerm = 0;
    
    unsigned n = getNumElongationMixtureElements();
    unsigned long Y = genome.getSumRFP();

    bool share_nse = shareNSE();

    fillMatrices(genome);
    //my_print("start!\n");
#ifdef _OPENMP
//#ifndef __APPLE__
#pragma omp parallel for private(gene,currAlpha,currLambda,currNSERate,propAlpha,propLambda,propNSERate,alphaCategory, lambdaCategory, nseCategory) reduction(+:logLikelihood,logLikelihood_proposed)
#endif
    for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
    {
    	unsigned long positionalRFPCount;
      int codonMixture, codonMixture_w_flag;
    	unsigned codonIndex;
    	std::string codon;
    	double currLgammaRFPAlpha;
        std::vector<std::vector<double>> prop_prob_successful(n);
        for (unsigned j = 0; j < n; j++ )
        {
        	prop_prob_successful[j] = std::vector<double>(getGroupListSize(),1000);
        }
        gene = &genome.getGene(i);

        unsigned mixtureElement = parameter->getMixtureAssignment(i);
        
        double U = getPartitionFunction(mixtureElement, false)/Y;
        
        std::vector <unsigned> positions = gene->geneData.getPositionCodonID();
        std::vector <unsigned long> rfpCounts = gene->geneData.getRFPCount(0);
        std::vector<int> positionMixture = gene->geneData.getPositionMixture();

        
        // how is the mixture element defined. Which categories make it up
        unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(mixtureElement);

        double phiValue = parameter->getSynthesisRate(i, synthesisRateCategory, false);

        double logPhi = std::log(phiValue); 
        
        double propSigma = 0;
        double currSigma = 0;
        
        for (unsigned positionIndex = 0; positionIndex < positions.size(); positionIndex++)
        {
        	  codonMixture_w_flag = positionMixture[positionIndex] + 1; // put back on 1-indexed scale to check if original value was negative or not
            if (codonMixture_w_flag < 0)
            {
              codonMixture = -1 * (codonMixture_w_flag) - 1; //if negative, get codonMixture if were not ignoring
            }
            else if (codonMixture_w_flag > 0)
            {
              codonMixture = codonMixture_w_flag - 1; //if positive, get codonMixture
            } else{
              my_print("ERROR: Your column indicating elongation mixtures contains 0. Should be a non-zero positive (include position in likelihood) or negative (use only for sigma) number. Exiting program.\n");
              exit(1);
            }
            codonIndex = positions[positionIndex];
            positionalRFPCount = rfpCounts[positionIndex];
            
            codon = gene->geneData.indexToCodon(codonIndex);

            alphaCategory = mixture_to_category[codonMixture][0];
            lambdaCategory = mixture_to_category[codonMixture][1];
            nseCategory = mixture_to_category[codonMixture][2];
            
            currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
            currLambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, codon, false);
            currNSERate = getParameterForCategory(nseCategory, PANSEParameter::nse, codon, false);

            if (positionalRFPCount < 50)
            {
                currLgammaRFPAlpha = lgamma_rfp_alpha[positionalRFPCount][alphaCategory][codonIndex];
            }
            else
            {
                currLgammaRFPAlpha = std::lgamma(currAlpha + positionalRFPCount);
            }
           
            if (share_nse && param == "NSERate")
            {

                propNSERate = getParameterForCategory(nseCategory, PANSEParameter::nse, codon, true);
                
                if (codonMixture_w_flag > 0)
                {
                  logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambda * U, positionalRFPCount,
                                phiValue,std::exp(propSigma),lgamma_currentAlpha[alphaCategory][codonIndex],log_currentLambda[synthesisRateCategory][lambdaCategory][codonIndex],logPhi,currLgammaRFPAlpha);
                }
                if (prop_prob_successful[codonMixture][codonIndex] > 500)
                {
                    prop_prob_successful[codonMixture][codonIndex] = elongationUntilIndexApproximation2ProbabilityLog(currAlpha, currLambda,1/propNSERate);
                    if (prop_prob_successful[codonMixture][codonIndex] > 0.0)
                    {
                        //prop_prob_successful[codonMixture][codonIndex] = std::numeric_limits<double>::quiet_NaN();
                        prop_prob_successful[codonMixture][codonIndex] = 0.0;
                    }
                }
                propSigma = propSigma + prop_prob_successful[codonMixture][codonIndex];
            }
            else if (codon == grouping)
            {
                if (param == "Elongation")
                {
                    propAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, true);
                    propLambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, codon, true);
                    if (codonMixture_w_flag > 0)
                    {
                      logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(propAlpha, propLambda * U, positionalRFPCount,
                                    phiValue,std::exp(propSigma),std::lgamma(propAlpha),std::log(propLambda) + std::log(U),logPhi,std::lgamma(propAlpha+positionalRFPCount));
                    }
                    if (prop_prob_successful[codonMixture][codonIndex] > 500)
                    {
                        prop_prob_successful[codonMixture][codonIndex] = elongationUntilIndexApproximation2ProbabilityLog(propAlpha, propLambda,1/currNSERate);
                        if (prop_prob_successful[codonMixture][codonIndex] > 0.0)
                        {
                            //prop_prob_successful[codonMixture][codonIndex] = std::numeric_limits<double>::quiet_NaN();
                            prop_prob_successful[codonMixture][codonIndex] = 0.0;
                        }
               
                    }
   
                }
                else if (param == "NSERate")
                {
                    propNSERate = getParameterForCategory(nseCategory, PANSEParameter::nse, codon, true);
                    if (codonMixture_w_flag > 0)
                    {
                      logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambda * U, positionalRFPCount,
                                    phiValue,std::exp(propSigma),lgamma_currentAlpha[alphaCategory][codonIndex],log_currentLambda[synthesisRateCategory][lambdaCategory][codonIndex],logPhi,currLgammaRFPAlpha);
                    }
                    if (prop_prob_successful[codonMixture][codonIndex] > 500)
                    {
                        prop_prob_successful[codonMixture][codonIndex] = elongationUntilIndexApproximation2ProbabilityLog(currAlpha, currLambda,1/propNSERate);
                        if (prop_prob_successful[codonMixture][codonIndex] > 0.0)
                        {
                        	//prop_prob_successful[codonMixture][codonIndex] = std::numeric_limits<double>::quiet_NaN();
                        	prop_prob_successful[codonMixture][codonIndex] = 0.0;
                        }
                    }
                }
                propSigma = propSigma + prop_prob_successful[codonMixture][codonIndex];
           }
           else
           {
                if (codonMixture_w_flag > 0)
                {
                  logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambda * U, positionalRFPCount,
                                  phiValue,std::exp(propSigma),lgamma_currentAlpha[alphaCategory][codonIndex],log_currentLambda[synthesisRateCategory][lambdaCategory][codonIndex],logPhi,currLgammaRFPAlpha);
                }
                propSigma = propSigma + prob_successful[codonMixture][codonIndex];
       
           }
           if (codonMixture_w_flag > 0)
           {
              logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambda * U, positionalRFPCount,
                                    phiValue, std::exp(currSigma),lgamma_currentAlpha[alphaCategory][codonIndex],log_currentLambda[synthesisRateCategory][lambdaCategory][codonIndex],logPhi,currLgammaRFPAlpha);
           }
            currSigma = currSigma + prob_successful[codonMixture][codonIndex];
        }
    }
	    

    for (unsigned j = 0; j < n; j++)
    {
        unsigned alphaCategory = parameter->getMutationCategory(j);
        unsigned lambdaCategory = parameter->getSelectionCategory(j);
        unsigned nseCategory = parameter->getNSECategory(j);
        if (param == "Elongation")
        {
            currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, grouping, false);
            currLambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, grouping, false);
            propAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, grouping, true);
            propLambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, grouping, true);
            if (std::isnan(logLikelihood_proposed))
            {
              my_print("WARNING: proposed logLikelihood for % is NaN\n",grouping);
              my_print("\tProposed alpha: % \n",getParameterForCategory(alphaCategory, PANSEParameter::alp, grouping, true));
              my_print("\tProposed lambda: %\n",getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, grouping, true));
            }
            currAdjustmentTerm += std::log(currAlpha) + std::log(currLambda);
            propAdjustmentTerm += std::log(propAlpha) + std::log(propLambda);
        }
        else
        {
            currNSERate = getParameterForCategory(nseCategory, PANSEParameter::nse, grouping, false);
            propNSERate = getParameterForCategory(nseCategory, PANSEParameter::nse, grouping, true);
            if (std::isnan(logLikelihood_proposed))
            {
               my_print("WARNING: proposed logLikelihood for % is NaN\n",grouping);
               my_print("\tProposed NSE Rate: %\n",getParameterForCategory(nseCategory, PANSEParameter::nse, grouping, true));
            }
            currAdjustmentTerm += std::log(currNSERate);
            propAdjustmentTerm += std::log(propNSERate);
        }
    }

    logPosterior_proposed = logLikelihood_proposed + calculateNSERatePrior(grouping,true) + calculateAlphaPrior(grouping,true) + calculateLambdaPrior(grouping,true);
    logPosterior = logLikelihood + calculateNSERatePrior(grouping,false) + calculateAlphaPrior(grouping,false) + calculateLambdaPrior(grouping,false);
    
    //Should never accept parameters that give NaN, so just check proposed parameters
    
  logAcceptanceRatioForAllMixtures[0] = logPosterior_proposed - logPosterior - (currAdjustmentTerm - propAdjustmentTerm);
	logAcceptanceRatioForAllMixtures[1] = logLikelihood;
	logAcceptanceRatioForAllMixtures[2] = logLikelihood_proposed;
	logAcceptanceRatioForAllMixtures[3] = logPosterior;
	logAcceptanceRatioForAllMixtures[4] = logPosterior_proposed;
	

    clearMatrices();
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
    double currAlpha, currLambda, currNSERate;
    double logLikelihood = 0;
    double logLikelihood_proposed = 0;
    unsigned alphaCategory, lambdaCategory, nseCategory;

    unsigned n = getNumMixtureElements(); //This is phi mixture (numMixtures), not elongation mixtures (numElongationMixtures)
    lpr = 0.0;

    unsigned long Y = genome.getSumRFP();
    fillMatrices(genome);
  
#ifdef _OPENMP
//#ifndef __APPLE__
#pragma omp parallel for private(gene,currAlpha,currLambda,currNSERate,alphaCategory,lambdaCategory,nseCategory) reduction(+:logLikelihood,logLikelihood_proposed)
#endif
    for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
    {

        unsigned long positionalRFPCount;
        unsigned codonIndex;
        int codonMixture, codonMixture_w_flag;
        std::string codon;
        double currLgammaRFPAlpha;
        
        gene = &genome.getGene(i);
        unsigned mixtureElement = parameter->getMixtureAssignment(i);

        // how is the mixture element defined. Which categories make it up
        double currU = getPartitionFunction(mixtureElement, false)/Y;
        double propU = getPartitionFunction(mixtureElement, true)/Y;


        unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(mixtureElement);
        // get non codon specific values, calculate likelihood conditional on these
        double phiValue = parameter->getSynthesisRate(i, synthesisRateCategory, false);

        std::vector <unsigned> positions = gene->geneData.getPositionCodonID();
        std::vector <unsigned long> rfpCounts = gene->geneData.getRFPCount(0);
        std::vector<int> positionMixture = gene->geneData.getPositionMixture();
    
        double logPhi = std::log(phiValue);
        
        double currSigma = 0;
        double propSigma = 0;

        for (unsigned positionIndex = 0; positionIndex < positions.size(); positionIndex++)
        {
            codonMixture_w_flag = positionMixture[positionIndex] + 1; // put back on 1-indexed scale to check if original value was negative or not
            if (codonMixture_w_flag < 0)
            {
              codonMixture = -1 * (codonMixture_w_flag) - 1; //if negative, get codonMixture if were not ignoring
            }
            else if (codonMixture_w_flag > 0)
            {
              codonMixture = codonMixture_w_flag - 1; //if positive, get codonMixture
            } else{
              my_print("ERROR: Your column indicating elongation mixtures contains 0. Should be a non-zero positive (include position in likelihood) or negative (use only for sigma) number. Exiting program.\n");
              exit(1);
            }
            positionalRFPCount = rfpCounts[positionIndex];
            codonIndex = positions[positionIndex];
            alphaCategory = mixture_to_category[codonMixture][0];
            lambdaCategory = mixture_to_category[codonMixture][1];
            nseCategory = mixture_to_category[codonMixture][2];

            codon = gene->geneData.indexToCodon(codonIndex);
            currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
            currLambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, codon, false);
            currNSERate = getParameterForCategory(nseCategory, PANSEParameter::nse, codon, false);
            
             //Check if lgamma(alpha + rfp_count) already counted for this codon 
            if (positionalRFPCount < 50)
            {
                currLgammaRFPAlpha = lgamma_rfp_alpha[positionalRFPCount][alphaCategory][codonIndex];
            }
            else
            {
                currLgammaRFPAlpha = std::lgamma(currAlpha + positionalRFPCount);
            }


            if (codonMixture_w_flag > 0)
            {
              logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, (currLambda * propU), positionalRFPCount,
                                   phiValue,std::exp(propSigma),lgamma_currentAlpha[alphaCategory][codonIndex],std::log(currLambda)+ std::log(propU),logPhi,currLgammaRFPAlpha);
              logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, (currLambda * currU), positionalRFPCount,
                                   phiValue,std::exp(currSigma),lgamma_currentAlpha[alphaCategory][codonIndex],log_currentLambda[synthesisRateCategory][lambdaCategory][codonIndex],logPhi,currLgammaRFPAlpha);
            }
            currSigma = currSigma + prob_successful[codonMixture][codonIndex];
            propSigma = propSigma + prob_successful[codonMixture][codonIndex];
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



bool PANSEModel::fixedAlpha()
{
    return(parameter->isAlphaFixed());
}

bool PANSEModel::fixedLambda()
{
    return(parameter->isLambdaFixed());
}

bool PANSEModel::fixedNSE()
{
    return(parameter->isNSEFixed());
}

bool PANSEModel::shareNSE()
{
    return(parameter->isNSEShared());
}

bool PANSEModel::getParameterTypeFixed(std::string csp_parameter)
{
    bool fixed = false;
    if (csp_parameter == parameter_types[0]) // == Elongation
    {
        bool alpha_fixed = fixedAlpha();
        bool lambda_fixed = fixedLambda();
        fixed = (alpha_fixed && lambda_fixed);
    }
    else if (csp_parameter == parameter_types[1]) // == NSERate
    {
        fixed = fixedNSE();
    }
    return(fixed);
}



bool PANSEModel::isShared(std::string csp_parameter)
{
    bool shared = false;
    if (csp_parameter == parameter_types[1])
    {
        shared = parameter->isNSEShared();
    }
    return(shared);
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


unsigned PANSEModel::getNumElongationMixtureElements()
{
    return parameter->getNumElongationMixtureElements();
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


void PANSEModel::updateCodonSpecificParameter(std::string codon)
{
    parameter->updateCodonSpecificParameter(codon,"Elongation");
    parameter->updateCodonSpecificParameter(codon,"NSERate");
}


void PANSEModel::updateCodonSpecificParameter(std::string codon,std::string param)
{
    parameter->updateCodonSpecificParameter(codon,param);
}


void PANSEModel::completeUpdateCodonSpecificParameter()
{
    parameter->completeUpdateCodonSpecificParameter();
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
    unsigned alphaCategory, lambdaCategory, nseCategory;
    unsigned long Y = genome.getSumRFP();
    my_print("##Total Number of Counts in Simulation: %\n",Y);
    double Z = 0;
    std::vector<std::vector<double>> wait_times;
    wait_times.resize(genome.getGenomeSize());
    
    for (unsigned geneIndex = 0u; geneIndex < genome.getGenomeSize(); geneIndex++)
    {

        unsigned mixtureElement = getMixtureAssignment(geneIndex);
        Gene gene = genome.getGene(geneIndex);
        double phi = parameter->getSynthesisRate(geneIndex, mixtureElement, false);
        SequenceSummary sequence = gene.geneData;
        Gene tmpGene = gene;
        std::vector <unsigned> positions = sequence.getPositionCodonID();
        std::vector<int> positionMixture = gene.geneData.getPositionMixture();
        
        wait_times[geneIndex].resize(positions.size());
        std::vector <unsigned> rfpCount;
        
        double sigma =  1.0;
        double v;
        for (unsigned positionIndex = 0; positionIndex < positions.size(); positionIndex++)
        {
            int codonMixture = positionMixture[positionIndex] + 1;
            if (codonMixture < 0)
            {
              codonMixture = -1 * (codonMixture) - 1;
            } 
            else
            {
              codonMixture = codonMixture - 1;
            }
            alphaCategory = parameter->getMutationCategory(codonMixture);
            lambdaCategory = parameter->getSelectionCategory(codonMixture);
            nseCategory = parameter->getNSECategory(codonMixture);
          
            unsigned codonIndex = positions[positionIndex];
            std::string codon = sequence.indexToCodon(codonIndex);
            if (codon == "TAG" || codon == "TGA" || codon == "TAA")
            {
                my_print("Stop codon being used during simulations\n");
            }
            double alpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, codon, false);
            double lambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, codon, false);
            double NSERate = getParameterForCategory(nseCategory, PANSEParameter::nse, codon, false);
            v = 1.0 / NSERate;
#ifndef STANDALONE
            RNGScope scope;
            NumericVector wt(1);
            wt = rgamma(1, alpha,1/(lambda));
            wait_times[geneIndex][positionIndex] = wt[0];
            Z += phi * wait_times[geneIndex][positionIndex] * sigma;
            sigma *= (v/(wait_times[geneIndex][positionIndex] + v));
#else
            std::gamma_distribution<double> GDistribution(alpha, 1.0/(lambda));
            wait_times[geneIndex][positionIndex] = GDistribution(Parameter::generator);
            Z += phi * wait_times[geneIndex][positionIndex] * sigma;
            sigma *= (v/(wait_times[geneIndex][positionIndex] + v));
#endif
            
        }
    }
    double U = Z/Y;
    my_print("##True Z value based on provided phi and CSPs:%\n",Z);
    my_print("Gene\tLog.Sigma\tSigma\n");
    for (unsigned geneIndex = 0u; geneIndex < genome.getGenomeSize(); geneIndex++)
    {

        unsigned mixtureElement = getMixtureAssignment(geneIndex);
        Gene gene = genome.getGene(geneIndex);
        double phi = parameter->getSynthesisRate(geneIndex, mixtureElement, false);
        SequenceSummary sequence = gene.geneData;
        Gene tmpGene = gene;
        std::vector <unsigned> positions = sequence.getPositionCodonID();
        std::vector<int> positionMixture = gene.geneData.getPositionMixture();
        std::vector <unsigned long> rfpCount;
        
        double sigma =  1.0;
        double v;
        
        for (unsigned positionIndex = 0; positionIndex < positions.size(); positionIndex++)
        {
            int codonMixture = positionMixture[positionIndex] + 1;
            if (codonMixture < 0)
            {
              codonMixture = -1 * (codonMixture) - 1;
            } 
            else
            {
              codonMixture = codonMixture - 1;
            }
            alphaCategory = parameter->getMutationCategory(codonMixture);
            lambdaCategory = parameter->getSelectionCategory(codonMixture);
            nseCategory = parameter->getNSECategory(codonMixture);
          
            unsigned codonIndex = positions[positionIndex];
            std::string codon = sequence.indexToCodon(codonIndex);
            if (codon == "TAG" || codon == "TGA" || codon == "TAA")
            {
                my_print("Stop codon being used during simulations\n");
            }
            double NSERate = getParameterForCategory(nseCategory, PANSEParameter::nse, codon, false);
            v = 1.0 / NSERate;
#ifndef STANDALONE
            RNGScope scope;
            NumericVector ribo_count(1);
            ribo_count = rpois(1, wait_times[geneIndex][positionIndex] * (1.0/U) * phi * sigma);
            sigma *= (v/(wait_times[geneIndex][positionIndex] + v));
            rfpCount.push_back(ribo_count[0]);
#else
            std::poisson_distribution<unsigned> PDistribution(phi * wait_times[geneIndex][positionIndex] * (1.0/U) * sigma);
            unsigned long simulatedValue = PDistribution(Parameter::generator);
            sigma *= (v/(wait_times[geneIndex][positionIndex] + v));
            rfpCount.push_back(simulatedValue);
#endif
        }
        std::cout << std::setprecision(15) << std::fixed;
        std::cout << gene.getId() << "\t" << std::log(sigma) << "\t" << sigma << "\n";
        tmpGene.geneData.setRFPCount(rfpCount, RFPCountColumn);
        genome.addGene(tmpGene, true);
    }
}


void PANSEModel::printHyperParameters()
{
    for (unsigned i = 0u; i < getNumSynthesisRateCategories(); i++)
    {
        my_print("stdDevSynthesisRate posterior estimate for selection category %: %\n", i, parameter->getStdDevSynthesisRate(i));
        my_print("partition function posterior estimate for selection category %: %\n", i, parameter->getPartitionFunction(i,false));
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

double PANSEModel::calculateAlphaPrior(std::string grouping,bool proposed)
{
    double priorValue = 0.0;

    unsigned numMutCat = parameter->getNumMutationCategories();
    for (unsigned i = 0u; i < numMutCat; i++)
    {
        double alpha = parameter->getParameterForCategory(i, PANSEParameter::alp, grouping, proposed);
        if (alpha < 0 || alpha > 100)
        {
            priorValue = std::log(0);
        }
        else
        {
            priorValue = std::log(1);
        }
    }
    return priorValue;
}

double PANSEModel::calculateLambdaPrior(std::string grouping,bool proposed)
{
    double priorValue = 0.0;

    unsigned numSelCat = parameter->getNumSelectionCategories();
    for (unsigned i = 0u; i < numSelCat; i++)
    {
        double lambda = parameter->getParameterForCategory(i, PANSEParameter::lmPri, grouping, proposed);
        if (lambda < 0 || lambda > 1000000000)
        {
            priorValue = std::log(0);
        }
        else
        {
            priorValue = std::log(1);
        }
    }
    return priorValue;
}


double PANSEModel::calculateNSERatePrior(std::string grouping,bool proposed)
{
  double NSERate,logNSERate;
  double priorValue = 0.0;
  double dist_mean = 25000;
  unsigned numMutCat = parameter->getNumNSECategories();
  for (unsigned i = 0u; i < numMutCat; i++)
  {
    NSERate = parameter->getParameterForCategory(i, PANSEParameter::nse, grouping, proposed);
    priorValue += (std::log(dist_mean) - dist_mean * NSERate);  
    
  }
  return priorValue;
}


double PANSEModel::calculateAllPriors(bool proposed)
{
  double prior = 0.0;
	unsigned size = getGroupListSize();
	bool share_nse = shareNSE();
	for (unsigned i = 0; i < size; i++)
	{
  	std::string grouping = getGrouping(i);
	  if (share_nse && i == 0)
	  {
  	  prior += calculateNSERatePrior(grouping, proposed);
	  } 
	  else if (share_nse && i > 0)
	  {
	    prior += 0;
	  }
	  else if (!share_nse)
	  {
	    prior += calculateNSERatePrior(grouping, proposed);
	  }
    prior += calculateAlphaPrior(grouping, proposed);
    prior += calculateLambdaPrior(grouping, proposed);
	}
	return prior;
}


bool PANSEModel::checkValues(bool proposed)
{
	unsigned alphaCategory,lambdaCategory,nseCategory;
	double currAlpha, currLambda, currNSERate;

   	double prior = 0.0;
   	double prob_success;
	unsigned size = getGroupListSize();
	mixture_to_category = getElongationMixtureCategories();
	bool good_values = true;
	for (unsigned i = 0; i < size; i++)
	{
		std::string grouping = getGrouping(i);
		prior += calculateNSERatePrior(grouping, proposed);
        prior += calculateAlphaPrior(grouping, proposed);
        prior += calculateLambdaPrior(grouping, proposed);
        if (std::isnan(prior) || !std::isfinite(prior))
        {
        	good_values = false;
        }
        for (unsigned j=0; j < mixture_to_category.size(); j++)
        {
        	alphaCategory = mixture_to_category[j][0];
        	lambdaCategory = mixture_to_category[j][1];
			nseCategory = mixture_to_category[j][2];

			currAlpha = getParameterForCategory(alphaCategory, PANSEParameter::alp, grouping, proposed);
			currLambda = getParameterForCategory(lambdaCategory, PANSEParameter::lmPri, grouping, proposed);
			currNSERate = getParameterForCategory(nseCategory, PANSEParameter::nse, grouping, proposed);
			prob_success = elongationUntilIndexApproximation2ProbabilityLog(currAlpha, currLambda, 1/currNSERate);
			if (prob_success > 0.0)
			{
				good_values = false;
			}
        }
	}

	return good_values;
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
        if(i % 2 == 0) rv = ((i / 2)) / (x + rv);
        if(i % 2 != 0) rv = ((((i / 2) + 1) - s)) / (1 + rv);
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
    d = std::log(UpperIncompleteGammaHelper(s, x));

    return rv - d;

}

//Log probability of elongation at current codon
double PANSEModel::elongationProbabilityLog(double currAlpha, double currLambda, double currNSE)
{ 
    double tmp = std::exp(std::log(currLambda) + std::log(currNSE));
    double x = tmp + currAlpha * std::log(currLambda)  + currAlpha * std::log(currNSE) + UpperIncompleteGammaLog(1.0-currAlpha,tmp);
    return x;
}


//Calculation of the probability of elongation at current codon
double PANSEModel::elongationProbability(double currAlpha, double currLambda, double currNSE)
{
    return std::pow(currLambda * currNSE, currAlpha) * std::exp(currLambda * currNSE) * UpperIncompleteGamma(1-currAlpha, currLambda * currNSE);
}







double PANSEModel::elongationUntilIndexApproximation1Probability(double alpha, double lambda, double v, double current)
{
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

double PANSEModel::elongationUntilIndexApproximation1ProbabilityLog(double alpha, double lambda, double v)
{
    return (-1*(alpha/(lambda * v)));   	   
}
double PANSEModel::elongationUntilIndexApproximation2ProbabilityLog(double alpha, double lambda, double v)
{
	return (-(alpha/(lambda * v)) + (alpha/(lambda * lambda * v * v))
                      + ((alpha/(lambda * v)) * (alpha/(lambda * v))) / 2);
}

double PANSEModel::calculateLogLikelihood(Genome &genome, std::vector<std::vector<double>> alpha, std::vector<std::vector<double>> lambda, std::vector<std::vector<double>> NSERate, std::vector<double> phi, double Z)
{
  double currAlpha, currLambda, currNSERate;
  unsigned alphaCategory, lambdaCategory, nseCategory;
  Gene *gene;
  
  std::vector<std::string> groups = parameter -> getGroupList();
  double logLikelihood = 0.0;
  unsigned n = getNumElongationMixtureElements();
  unsigned long Y = genome.getSumRFP();
  
  double U = Z/Y;
  for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
  {
    unsigned long positionalRFPCount;
    unsigned codonMixture;
    unsigned codonIndex;
    std::string codon;
    double currLgammaRFPAlpha;
    std::vector<std::vector<double>> prob_successful(n);
    for (unsigned j = 0; j < n; j++ )
    {
      prob_successful[j] = std::vector<double>(getGroupListSize(),1000);
    }
    gene = &genome.getGene(i);
    
    
    std::vector <unsigned> positions = gene->geneData.getPositionCodonID();
    std::vector <unsigned long> rfpCounts = gene->geneData.getRFPCount(0);
    std::vector<unsigned> positionMixture = gene->geneData.getPositionMixture();
    
    double phiValue = phi[i];
    double currSigma = 0;
    
    for (unsigned positionIndex = 0; positionIndex < positions.size(); positionIndex++)
    {
      codonMixture = positionMixture[positionIndex];
      
      positionalRFPCount = rfpCounts[positionIndex];
      codonIndex = positions[positionIndex];
      codon = gene->geneData.indexToCodon(codonIndex);
      
      currAlpha = alpha[codonMixture][codonIndex];
      currLambda = lambda[codonMixture][codonIndex];
      currNSERate = NSERate[codonMixture][codonIndex];
      logLikelihood += calculateLogLikelihoodPerCodonPerGeneByPosition(currAlpha, currLambda * U, positionalRFPCount,
                                                             phiValue, std::exp(currSigma));
      if (prob_successful[codonMixture][codonIndex] > 500)
      {
        prob_successful[codonMixture][codonIndex] = elongationUntilIndexApproximation2ProbabilityLog(currAlpha, currLambda,1/currNSERate);
        if (prob_successful[codonMixture][codonIndex] > 0.0)
        {
          prob_successful[codonMixture][codonIndex] = std::numeric_limits<double>::quiet_NaN();
          my_print("Warning: greater than 1. Prob now %\n",prob_successful[codonMixture][codonIndex]);
        }
      }
      currSigma = currSigma + prob_successful[codonMixture][codonIndex];
      
    }
  }
  
  return(logLikelihood);
}


#ifndef STANDALONE

double PANSEModel::calculateLogLikelihoodR(Genome &genome, std::vector<double> _alpha, std::vector<double> _lambda, std::vector<double> _NSERate, std::vector<double> _phi, double _Z)
{
  std::vector<std::vector<double>> alpha,lambda,NSERate;
  double loglikelihood;
  unsigned numElongationMixtures;
  std::vector<std::string> groups = parameter -> getGroupList();
  unsigned numParam = groups.size();
  if ((_alpha.size() % numParam != 0) || (_lambda.size() % numParam != 0) || (_NSERate.size() % numParam != 0))
  {
    my_print("ERROR: alpha, lambda, and NSERate must be a multiple of the number of codons (61).\n");
    std::exit(1);
  }
  if (_alpha.size() != _lambda.size() || _alpha.size() != _NSERate.size())
  {
    my_print("ERROR: alpha, lambda, and NSERate must be the same size. If you are trying to share the same parameter values across elongation mixture, then these values should be included twice.\n");
    std::exit(1);
  }
  numElongationMixtures = _alpha.size()/numParam;
  alpha.resize(numElongationMixtures);
  lambda.resize(numElongationMixtures);
  NSERate.resize(numElongationMixtures);
  unsigned index = 0;
  for (unsigned i = 0; i < numElongationMixtures; i++)
  {
    alpha[i].resize(numParam);
    lambda[i].resize(numParam);
    NSERate[i].resize(numParam);
    for (unsigned j = 0; j < numParam; j++,index++)
    {
      alpha[i][j] = _alpha[index];
      lambda[i][j] = _lambda[index];
      NSERate[i][j] = _NSERate[index];
    }
  }
  loglikelihood = calculateLogLikelihood(genome,alpha,lambda,NSERate,_phi,_Z);
  return(loglikelihood);
}

#endif //STANDALONE		  

