#include "include/PA/PAModel.h"

//R runs only
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

//--------------------------------------------------//
//----------- Constructors & Destructors ---------- //
//--------------------------------------------------//

PAModel::PAModel(unsigned _RFPCountColumn,bool _withPhi, bool _fix_sEpsilon) : Model()
{
	parameter = NULL;
	RFPCountColumn = _RFPCountColumn - 1;
	withPhi = _withPhi;
	fix_sEpsilon = _fix_sEpsilon;
	parameter_types = {"Elongation"};
}


PAModel::~PAModel()
{
	//dtor
	//TODO: call Parent's deconstructor
	//delete parameter;
}


double PAModel::calculateAlphaPrior(std::string grouping,bool proposed)
{
    double priorValue = 0.0;

    unsigned numMutCat = parameter->getNumMutationCategories();
    for (unsigned i = 0u; i < numMutCat; i++)
    {
        double alpha = parameter->getParameterForCategory(i, PAParameter::alp, grouping, proposed);
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

double PAModel::calculateLambdaPrior(std::string grouping,bool proposed)
{
    double priorValue = 0.0;

    unsigned numSelCat = parameter->getNumSelectionCategories();
    for (unsigned i = 0u; i < numSelCat; i++)
    {
        double lambda = parameter->getParameterForCategory(i, PAParameter::lmPri, grouping, proposed);
        if (lambda < 0 || lambda > 100)
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



double PAModel::calculateLogLikelihoodPerCodonPerGene(double currAlpha, double currLambdaPrime, unsigned currRFPValue,
                                                      unsigned currNumCodonsInMRNA, double phiValue)
{
	
	double logLikelihood = ((std::lgamma((currNumCodonsInMRNA * currAlpha) + currRFPValue)) - (std::lgamma(currNumCodonsInMRNA * currAlpha)))
						   + (currRFPValue * (std::log(phiValue) - std::log(currLambdaPrime + phiValue)))
						   + ((currNumCodonsInMRNA * currAlpha) * (std::log(currLambdaPrime) - std::log(currLambdaPrime + phiValue)));

	return logLikelihood;
}





//------------------------------------------------//
//---------- Likelihood Ratio Functions ----------//
//------------------------------------------------//


void PAModel::calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k, double* logProbabilityRatio)
{
	double logLikelihood = 0.0;
	double logLikelihood_proposed = 0.0;

	unsigned alphaCategory = parameter->getMutationCategory(k);
	unsigned lambdaPrimeCategory = parameter->getSelectionCategory(k);
	unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(k);

	double phiValue = parameter->getSynthesisRate(geneIndex, synthesisRateCategory, false);
	double phiValue_proposed = parameter->getSynthesisRate(geneIndex, synthesisRateCategory, true);

// #ifdef _OPENMP
// //#ifndef __APPLE__
// #pragma omp parallel for reduction(+:logLikelihood,logLikelihood_proposed)
// #endif
	for (unsigned index = 0; index < getGroupListSize(); index++) //number of codons, without the stop codons
	{
		std::string codon = getGrouping(index);
		double currAlpha = getParameterForCategory(alphaCategory, PAParameter::alp, codon, false);
		double currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PAParameter::lmPri, codon, false);
		unsigned currRFPValue = gene.geneData.getCodonSpecificSumRFPCount(index,0 /*RFPCountColumn*/);
		unsigned currNumCodonsInMRNA = gene.geneData.getCodonCountForCodon(index);
		if (currNumCodonsInMRNA == 0) 
		{	
			continue;
		}
		logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, currRFPValue, currNumCodonsInMRNA, phiValue);
		logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, currRFPValue, currNumCodonsInMRNA, phiValue_proposed);
	}

	unsigned mixture = getMixtureAssignment(geneIndex);
	mixture = getSynthesisRateCategory(mixture);
	double stdDevSynthesisRate = parameter->getStdDevSynthesisRate(mixture, false);
	double mPhi = (-(stdDevSynthesisRate * stdDevSynthesisRate) * 0.5); // X * 0.5 = X / 2
	double logPhiProbability = Parameter::densityLogNorm(phiValue, mPhi, stdDevSynthesisRate, true);
	double logPhiProbability_proposed = Parameter::densityLogNorm(phiValue_proposed, mPhi, stdDevSynthesisRate, true);


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


void PAModel::calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome,
                                                                std::vector<double> &logAcceptanceRatioForAllMixtures,std::string param)
{
	double currAlpha,currLambdaPrime,propAlpha,propLambdaPrime;
	double logLikelihood, logPosterior = 0.0;
    double logLikelihood_proposed, logPosterior_proposed = 0.0;
	Gene *gene;
	unsigned index = SequenceSummary::codonToIndex(grouping);
	unsigned n = getNumMixtureElements();
	double currAdjustmentTerm = 0;
    double propAdjustmentTerm = 0;
#ifdef _OPENMP
//#ifndef __APPLE__
#pragma omp parallel for private(gene,currAlpha,currLambdaPrime,propAlpha,propLambdaPrime) reduction(+:logLikelihood,logLikelihood_proposed)
#endif
    for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		gene = &genome.getGene(i);
		
		unsigned mixtureElement = parameter->getMixtureAssignment(i);
		// how is the mixture element defined. Which categories make it up
		unsigned alphaCategory = parameter->getMutationCategory(mixtureElement);
		unsigned lambdaPrimeCategory = parameter->getSelectionCategory(mixtureElement);
		unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(mixtureElement);
		
		double phiValue = parameter->getSynthesisRate(i, synthesisRateCategory, false);
		
		unsigned currRFPValue = gene->geneData.getCodonSpecificSumRFPCount(index, RFPCountColumn);
		unsigned currNumCodonsInMRNA = gene->geneData.getCodonCountForCodon(index);
		if (currNumCodonsInMRNA == 0) continue;

		currAlpha = getParameterForCategory(alphaCategory, PAParameter::alp, grouping, false);
		currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PAParameter::lmPri, grouping, false);

		propAlpha = getParameterForCategory(alphaCategory, PAParameter::alp, grouping, true);
		propLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PAParameter::lmPri, grouping, true);

		logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, currRFPValue, currNumCodonsInMRNA, phiValue);
		logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(propAlpha, propLambdaPrime, currRFPValue, currNumCodonsInMRNA, phiValue);

	}
	for (unsigned j = 0; j < n; j++)
    {
        unsigned alphaCategory = parameter->getMutationCategory(j);
        unsigned lambdaPrimeCategory = parameter->getSelectionCategory(j);
        currAlpha = getParameterForCategory(alphaCategory, PAParameter::alp, grouping, false);
        currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PAParameter::lmPri, grouping, false);
        propAlpha = getParameterForCategory(alphaCategory, PAParameter::alp, grouping, true);
        propLambdaPrime = getParameterForCategory(lambdaPrimeCategory, PAParameter::lmPri, grouping, true);

        currAdjustmentTerm += std::log(currAlpha) + std::log(currLambdaPrime);
        propAdjustmentTerm += std::log(propAlpha) + std::log(propLambdaPrime);
    }

    logPosterior_proposed = logLikelihood_proposed + calculateAlphaPrior(grouping,true) + calculateLambdaPrior(grouping,true);
    logPosterior = logLikelihood + calculateAlphaPrior(grouping,false) + calculateLambdaPrior(grouping,false);
    
	
	logAcceptanceRatioForAllMixtures[0] = logPosterior_proposed - logPosterior - (currAdjustmentTerm - propAdjustmentTerm);
	logAcceptanceRatioForAllMixtures[3] = logLikelihood;
	logAcceptanceRatioForAllMixtures[4] = logLikelihood_proposed;
	logAcceptanceRatioForAllMixtures[1] = logPosterior;
	logAcceptanceRatioForAllMixtures[2] = logPosterior_proposed;
}


void PAModel::calculateLogLikelihoodRatioForHyperParameters(Genome &genome, unsigned iteration, std::vector <double> & logProbabilityRatio)
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
		currentMphi[i] = -((currentStdDevSynthesisRate[i] * currentStdDevSynthesisRate[i]) * 0.5);
		proposedStdDevSynthesisRate[i] = getStdDevSynthesisRate(i, true);
		proposedMphi[i] = -((proposedStdDevSynthesisRate[i] * proposedStdDevSynthesisRate[i]) * 0.5);
		// take the Jacobian into account for the non-linear transformation from logN to N distribution
		lpr -= (std::log(currentStdDevSynthesisRate[i]) - std::log(proposedStdDevSynthesisRate[i]));
		// take prior into account
		//TODO(Cedric): make sure you can control that prior from R
		//lpr -= Parameter::densityNorm(currentStdDevSynthesisRate[i], 1.0, 0.1, true) - Parameter::densityNorm(proposedStdDevSynthesisRate[i], 1.0, 0.1, true);
	}

	if (withPhi)
    {
        // one for each noiseOffset, and one for stdDevSynthesisRate
        logProbabilityRatio.resize(getNumPhiGroupings() + 1);
    }
    else
        logProbabilityRatio.resize(1);

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
            logProbabilityRatio[i+1] = lpr;
            my_print("this should not be here\n");
        }
    }
}


//----------------------------------------------------------//
//---------- Initialization and Restart Functions ----------//
//----------------------------------------------------------//


void PAModel::initTraces(unsigned samples, unsigned num_genes, bool estimateSynthesisRate)
{
	parameter->initAllTraces(samples, num_genes, estimateSynthesisRate);
}


void PAModel::writeRestartFile(std::string filename)
{
	return parameter->writeEntireRestartFile(filename);
}





//----------------------------------------//
//---------- Category Functions ----------//
//----------------------------------------//


double PAModel::getCategoryProbability(unsigned i)
{
	return parameter->getCategoryProbability(i);
}


unsigned PAModel::getMutationCategory(unsigned mixture)
{
	return parameter->getMutationCategory(mixture);
}


unsigned PAModel::getSelectionCategory(unsigned mixture)
{
	return parameter->getSelectionCategory(mixture);
}


unsigned PAModel::getSynthesisRateCategory(unsigned mixture)
{
	return parameter->getSynthesisRateCategory(mixture);
}


std::vector<unsigned> PAModel::getMixtureElementsOfSelectionCategory(unsigned k)
{
	return parameter->getMixtureElementsOfSelectionCategory(k);
}





//------------------------------------------//
//---------- Group List Functions ----------//
//------------------------------------------//


unsigned PAModel::getGroupListSize()
{
	return parameter->getGroupListSize();
}


std::string PAModel::getGrouping(unsigned index)
{
	return parameter->getGrouping(index);
}





//---------------------------------------------------//
//---------- stdDevSynthesisRate Functions ----------//
//---------------------------------------------------//


double PAModel::getStdDevSynthesisRate(unsigned selectionCategory, bool proposed)
{
	return parameter->getStdDevSynthesisRate(selectionCategory, proposed);
}


double PAModel::getCurrentStdDevSynthesisRateProposalWidth()
{
	return parameter->getCurrentStdDevSynthesisRateProposalWidth();
}


void PAModel::updateStdDevSynthesisRate()
{
	parameter->updateStdDevSynthesisRate();
}





//----------------------------------------------//
//---------- Synthesis Rate Functions ----------//
//----------------------------------------------//


double PAModel::getSynthesisRate(unsigned index, unsigned mixture, bool proposed)
{
	return parameter->getSynthesisRate(index, mixture, proposed);
}


void PAModel::updateSynthesisRate(unsigned i, unsigned k)
{
	parameter->updateSynthesisRate(i, k);
}





//-----------------------------------------//
//---------- Iteration Functions ----------//
//-----------------------------------------//


unsigned PAModel::getLastIteration()
{
	return parameter->getLastIteration();
}


void PAModel::setLastIteration(unsigned iteration)
{
	parameter->setLastIteration(iteration);
}





//-------------------------------------//
//---------- Trace Functions ----------//
//-------------------------------------//


void PAModel::updateStdDevSynthesisRateTrace(unsigned sample)
{
	parameter->updateStdDevSynthesisRateTrace(sample);
}


void PAModel::updateSynthesisRateTrace(unsigned sample, unsigned i)
{
	parameter->updateSynthesisRateTrace(sample, i);
}


void PAModel::updateMixtureAssignmentTrace(unsigned sample, unsigned i)
{
	parameter->updateMixtureAssignmentTrace(sample, i);
}


void PAModel::updateMixtureProbabilitiesTrace(unsigned sample)
{
	parameter->updateMixtureProbabilitiesTrace(sample);
}


void PAModel::updateCodonSpecificParameterTrace(unsigned sample, std::string codon)
{
	parameter->updateCodonSpecificParameterTrace(sample, codon);
}


void PAModel::updateHyperParameterTraces(unsigned sample)
{
	updateStdDevSynthesisRateTrace(sample);
	if(withPhi)
	{
		updateNoiseOffsetTrace(sample);
    	updateObservedSynthesisNoiseTrace(sample);
    }
}


void PAModel::updateTracesWithInitialValues(Genome & genome)
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


void PAModel::adaptStdDevSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt)
{
	parameter->adaptStdDevSynthesisRateProposalWidth(adaptiveWidth, adapt);
}


void PAModel::adaptSynthesisRateProposalWidth(unsigned adaptiveWidth, bool adapt)
{
	parameter->adaptSynthesisRateProposalWidth(adaptiveWidth, adapt);
}


void PAModel::adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth, unsigned lastIteration, bool adapt)
{
	parameter->adaptCodonSpecificParameterProposalWidth(adaptiveWidth, lastIteration, adapt);
}


void PAModel::adaptHyperParameterProposalWidths(unsigned adaptiveWidth, bool adapt)
{
	adaptStdDevSynthesisRateProposalWidth(adaptiveWidth, adapt);
	if (withPhi)
		adaptNoiseOffsetProposalWidth(adaptiveWidth, adapt);
}





//-------------------------------------//
//---------- Other Functions ----------//
//-------------------------------------//


void PAModel::proposeCodonSpecificParameter()
{
	parameter->proposeCodonSpecificParameter();
}


void PAModel::proposeHyperParameters()
{
	parameter->proposeStdDevSynthesisRate();
	if (withPhi)
	{
		parameter->proposeNoiseOffset();
	}
}


void PAModel::proposeSynthesisRateLevels()
{
	parameter->proposeSynthesisRateLevels();
}


unsigned PAModel::getNumPhiGroupings()
{
	return parameter->getNumObservedPhiSets();
}


unsigned PAModel::getMixtureAssignment(unsigned index)
{
	return parameter->getMixtureAssignment(index);
}


unsigned PAModel::getNumMixtureElements()
{
	return parameter->getNumMixtureElements();
}


unsigned PAModel::getNumSynthesisRateCategories()
{
	return parameter->getNumSynthesisRateCategories();
}


void PAModel::setNumPhiGroupings(unsigned value)
{
	parameter->setNumObservedPhiSets(value);
}


void PAModel::setMixtureAssignment(unsigned i, unsigned catOfGene)
{
	parameter->setMixtureAssignment(i, catOfGene);
}


void PAModel::setCategoryProbability(unsigned mixture, double value)
{
	parameter->setCategoryProbability(mixture, value);
}


void PAModel::updateCodonSpecificParameter(std::string codon)
{
	parameter->updateCodonSpecificParameter(codon);
}

void PAModel::updateCodonSpecificParameter(std::string codon, std::string param)
{
	parameter->updateCodonSpecificParameter(codon);
}

void PAModel::completeUpdateCodonSpecificParameter()
{
    parameter->completeUpdateCodonSpecificParameter();
}


//Noise offset functions

// double PAModel::getNoiseOffset(unsigned index, bool proposed)
// {
// 	return parameter->getNoiseOffset(index, proposed);
// }


// double PAModel::getObservedSynthesisNoise(unsigned index)
// {
// 	return parameter->getObservedSynthesisNoise(index);
// }


// double PAModel::getCurrentNoiseOffsetProposalWidth(unsigned index)
// {
// 	return parameter->getCurrentNoiseOffsetProposalWidth(index);
// }


// void PAModel::updateNoiseOffset(unsigned index)
// {
// 	parameter->updateNoiseOffset(index);
// }


// void PAModel::updateNoiseOffsetTrace(unsigned sample)
// {
// 	parameter->updateNoiseOffsetTraces(sample);
// }


// void PAModel::updateObservedSynthesisNoiseTrace(unsigned sample)
// {
// 	parameter->updateObservedSynthesisNoiseTraces(sample);
// }


// void PAModel::adaptNoiseOffsetProposalWidth(unsigned adaptiveWidth, bool adapt)
// {
// 	parameter->adaptNoiseOffsetProposalWidth(adaptiveWidth, adapt);
// }



// void PAModel::updateGibbsSampledHyperParameters(Genome &genome)
// {
//   // estimate s_epsilon by sampling from a gamma distribution and transforming it into an inverse gamma sample
	
// 	if (withPhi)
// 	{
// 		if(!fix_sEpsilon)
// 		{
// 			double shape = ((double)genome.getGenomeSize() - 1.0) / 2.0;
// 			for (unsigned i = 0; i < parameter->getNumObservedPhiSets(); i++)
// 			{
// 				double rate = 0.0; //Prior on s_epsilon goes here?
// 				unsigned mixtureAssignment;
// 				double noiseOffset = getNoiseOffset(i);
// 				for (unsigned j = 0; j < genome.getGenomeSize(); j++)
// 				{
// 					mixtureAssignment = parameter->getMixtureAssignment(j);
// 					double obsPhi = genome.getGene(j).getObservedSynthesisRate(i);
// 					if (obsPhi > -1.0)
// 					{
// 						double sum = std::log(obsPhi) - noiseOffset - std::log(parameter->getSynthesisRate(j, mixtureAssignment, false));
// 						rate += (sum * sum);
// 					}else{
// 						// missing observation.
// 						shape -= 0.5;
// 						//Reduce shape because initial estimate assumes there are no missing observations
// 					}
// 				}
// 				rate /= 2.0;
// 				double rand = parameter->randGamma(shape, rate);

// 				// Below the gamma sample is transformed into an inverse gamma sample
// 				// According to Gilchrist et al (2015) Supporting Materials p. S6
// 				// The sample 1/T is supposed to be equal to $s_\epsilon^2$.
// 				double sepsilon = std::sqrt(1.0/rand);
// 				parameter->setObservedSynthesisNoise(i, sepsilon);
// 			}
// 		}
// 	}
// }


void PAModel::updateAllHyperParameter()
{
	updateStdDevSynthesisRate();
	if (withPhi)
	{
		for (unsigned i =0; i < parameter->getNumObservedPhiSets(); i++)
		{
			updateNoiseOffset(i);
		}
	}
}


void PAModel::updateHyperParameter(unsigned hp)
{
	// NOTE: when adding additional hyper parameter, also add to updateAllHyperParameter()
	if (hp == 0)
	{
		updateStdDevSynthesisRate();
	}
	else if (hp > 0 && withPhi)
	{
		updateNoiseOffset(hp-1);
	}		
}


void PAModel::simulateGenome(Genome &genome)
{
	  unsigned long Y = genome.getSumRFP();
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
            unsigned alphaCategory = parameter->getMutationCategory(codonMixture);
            unsigned lambdaCategory = parameter->getSelectionCategory(codonMixture);
            unsigned codonIndex = positions[positionIndex];
            std::string codon = sequence.indexToCodon(codonIndex);
            if (codon == "TAG" || codon == "TGA" || codon == "TAA")
            {
                my_print("Stop codon being used during simulations\n");
            }
            double alpha = getParameterForCategory(alphaCategory, PAParameter::alp, codon, false);
            double lambda = getParameterForCategory(lambdaCategory, PAParameter::lmPri, codon, false);
#ifndef STANDALONE
            RNGScope scope;
            NumericVector wt(1);
            wt = rgamma(1, alpha,1/(lambda));
            wait_times[geneIndex][positionIndex] = wt[0];
            Z += phi * wait_times[geneIndex][positionIndex];
#else
            std::gamma_distribution<double> GDistribution(alpha, 1.0/(lambda));
            wait_times[geneIndex][positionIndex] = GDistribution(Parameter::generator);
            Z += phi * wait_times[geneIndex][positionIndex];
#endif
            
        }
    }
    double U = Z/Y;
    my_print("##True Z value based on provided phi and CSPs:%\n",Z);
    for (unsigned geneIndex = 0u; geneIndex < genome.getGenomeSize(); geneIndex++)
    {
        unsigned mixtureElement = getMixtureAssignment(geneIndex);
        Gene gene = genome.getGene(geneIndex);
        double phi = parameter->getSynthesisRate(geneIndex, mixtureElement, false);
        SequenceSummary sequence = gene.geneData;
        Gene tmpGene = gene;
        std::vector <unsigned> positions = sequence.getPositionCodonID();
        std::vector <unsigned long> rfpCount;
        std::vector<int> positionMixture = gene.geneData.getPositionMixture();
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
            unsigned alphaCategory = parameter->getMutationCategory(codonMixture);
            unsigned lambdaCategory = parameter->getSelectionCategory(codonMixture);
            unsigned codonIndex = positions[positionIndex];
            std::string codon = sequence.indexToCodon(codonIndex);
            if (codon == "TAG" || codon == "TGA" || codon == "TAA")
            {
                my_print("Stop codon being used during simulations\n");
            }
#ifndef STANDALONE
            RNGScope scope;
            NumericVector ribo_count(1);
            ribo_count = rpois(1, wait_times[geneIndex][positionIndex] * (1.0/U) * phi);
            rfpCount.push_back(ribo_count[0]);
#else
            std::poisson_distribution<unsigned> PDistribution(phi * wait_times[geneIndex][positionIndex] * (1.0/U));
            unsigned long simulatedValue = PDistribution(Parameter::generator);
            rfpCount.push_back(simulatedValue);
#endif
        }
        tmpGene.geneData.setRFPCount(rfpCount, RFPCountColumn);
        genome.addGene(tmpGene, true);
    }
}


void PAModel::printHyperParameters()
{
	for (unsigned i = 0u; i < getNumSynthesisRateCategories(); i++)
	{
		my_print("stdDevSynthesisRate posterior estimate for selection category %: %\n", i, parameter->getStdDevSynthesisRate(i));
	}
	my_print("\t current stdDevSynthesisRate proposal width: %\n", getCurrentStdDevSynthesisRateProposalWidth());
    //printCodonSpecificParameters(); //TODO put this in MCMC instead
}

////TODO: Assumed single mixture correct this and label values
//void PAModel::printCodonSpecificParameters()
//{
//    std::vector<std::vector<double>> alphas = parameter-> getCurrentAlphaParameter();
//    std::vector<std::vector<double>> lambdaPrimes = parameter-> getCurrentLambdaPrimeParameter();
//
//	for (unsigned i = 0u; i < alphas.size(); i++)
//	{
//        for (unsigned j = 0u; j < alphas[i].size(); j++)
//        {
//    		my_print("Alpha estimate for selection category %: %\n", i, alphas[i][j]);
//	    	my_print("Lambda Prime  estimate for selection category %: %\n", i, lambdaPrimes[i][j]);
//        }
//	}
//}

/* getParameter (RCPP EXPOSED)
 * Arguments: None
 *
 * Returns the PAParameter of the model.
*/
PAParameter* PAModel::getParameter()
{
	return parameter;
}


/* setParameter (RCPP EXPOSED)
 * Arguments: The PAParameter object to be referenced and stored.
 * Sets the private parameter object to the parameter object specified.
*/
void PAModel::setParameter(PAParameter &_parameter)
{
	parameter = &_parameter;
}


double PAModel::calculateAllPriors(bool proposed)
{
	double prior = 0.0;
	unsigned size = getGroupListSize();

	for (unsigned i = 0; i < size; i++)
	{
		std::string grouping = getGrouping(i);
        prior += calculateAlphaPrior(grouping, proposed);
        prior += calculateLambdaPrior(grouping, proposed);
	}

	// add more priors if necessary.

	return prior;
}


double PAModel::getParameterForCategory(unsigned category, unsigned param, std::string codon, bool proposal)
{
	return parameter->getParameterForCategory(category, param, codon, proposal);
}


bool PAModel::getParameterTypeFixed(std::string csp_parameter)
{
    bool fixed = false;
    if (csp_parameter == parameter_types[0]) // == Evolutionary
    {
        bool alpha_fixed = parameter -> isAlphaFixed();
        bool lambda_fixed = parameter -> isLambdaFixed();
        fixed == alpha_fixed && lambda_fixed;
    }
    return(fixed);
}

bool PAModel::isShared(std::string csp_parameters)
{
	return false;
}

