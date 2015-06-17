#include "include/MCMCAlgorithm.h"
#include "include/CovarianceMatrix.h"

#include <random>
#include <cstdlib>
#include <thread>
#include <sstream>

#include <iostream>
#include <fstream>
#include <stdlib.h> //can be removed later

//#include <omp.h>



MCMCAlgorithm::MCMCAlgorithm() : samples(1000), thining(1), adaptiveWidth(100 * thining), estimateExpression(true),
estimateCodonSpecificParameter(true), estimateHyperParameter(true)
{
    likelihoodTrace.resize(samples);
}

MCMCAlgorithm::MCMCAlgorithm(int _samples, int _thining, int _adaptiveWidth, bool _estimateExpression, bool _estimateCodonSpecificParameter, bool _estimateHyperParameter)
    : samples(_samples), thining(_thining), adaptiveWidth(_adaptiveWidth * thining), estimateExpression(_estimateExpression), estimateCodonSpecificParameter(_estimateCodonSpecificParameter),
        estimateHyperParameter(_estimateHyperParameter)
{
    likelihoodTrace.resize(samples);
}

MCMCAlgorithm::~MCMCAlgorithm()
{
    //dtor
}

MCMCAlgorithm::MCMCAlgorithm(const MCMCAlgorithm& other)
{
        samples = other.samples;
        thining = other.thining;
        adaptiveWidth = other.adaptiveWidth;

        estimateExpression = other.estimateExpression;
        estimateCodonSpecificParameter = other.estimateCodonSpecificParameter;
        estimateHyperParameter = other.estimateHyperParameter;

        likelihoodTrace = other.likelihoodTrace;
}

double MCMCAlgorithm::acceptRejectExpressionLevelForAllGenes(Genome& genome, ROCParameter& parameter, ROCModel& model, int iteration)
{
    // TODO move the likelihood calculation out off here. make it a void function again.

    double logLikelihood = 0.0;
    int numGenes = genome.getGenomeSize();

    // just for testing
    //unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
    //omp_set_num_threads(concurentThreadsSupported);
    //#pragma omp parallel for //shared(parameter)
    // testing end

    unsigned numExpressionCategories = parameter.getNumExpressionCategories();
    unsigned numMixtures = parameter.getNumMixtureElements();
    double* dirichletParameters = new double[numMixtures]();
		//initialize parameter's size
    for(int i = 0; i < numGenes; i++)
    {
        Gene gene = genome.getGene(i);
        /*
            Since some values returned by calculateLogLiklihoodRatioPerGene are veyr small (~ -1100), exponantiation leads to 0.
            To solve this problem, we adjust the value by a constant c. I choose to use the average value across all mixtures.
            We justify this by
            P = Sum(p_i*f(...))
                => f' = c*f
                => ln(f') = ln(c) + ln(f)
                => ln(P) = ln( Sum(p_i*f'(...)) )
                => ln(P) = ln(P') - ln(c)
            Note that we use invere sign because our values of ln(f) and ln(f') are negative.
        */
        double* unscaledLogProb_curr = new double[numExpressionCategories]();
        double* unscaledLogProb_prop = new double[numExpressionCategories]();

        double* unscaledLogProb_curr_singleMixture = new double[numMixtures]();
        double maxValue = -1000000.0;
        unsigned mixtureIndex = 0u;
        for(unsigned k = 0u; k < numExpressionCategories; k++)
        {
            // logProbabilityRatio contains the logProbabilityRatio in element 0,
            // the current unscaled probability in elemen 1 and the proposed unscaled probability in element 2
            std::vector<unsigned> mixtureElements = parameter.getMixtureElementsOfSelectionCategory(k);
            for(unsigned n = 0u; n < mixtureElements.size(); n++)
            {
                unsigned mixtureElement = mixtureElements[n];
                double logProbabilityRatio[3];
                model.calculateLogLiklihoodRatioPerGene(gene, i, parameter, mixtureElement, logProbabilityRatio);

                // store values so they can be processed
                unscaledLogProb_curr[k] += logProbabilityRatio[1];
                unscaledLogProb_prop[k] += logProbabilityRatio[2];

                unscaledLogProb_curr_singleMixture[mixtureIndex] = logProbabilityRatio[1];
                maxValue = unscaledLogProb_curr_singleMixture[mixtureIndex] > maxValue ? unscaledLogProb_curr_singleMixture[mixtureIndex] : maxValue;
                mixtureIndex++;
            }
        }
        unsigned mixtureAssignmentOfGene = parameter.getMixtureAssignment(i);
        for(unsigned k = 0u; k < numExpressionCategories; k++)
        {
            // We do not need to add std::log(parameter.getCategoryProbability(k)) since it will cancel in the ratio!
            double currLogLike = unscaledLogProb_curr[k];
            double propLogLike = unscaledLogProb_prop[k];
            if( -ROCParameter::randExp(1) < (propLogLike - currLogLike) )
            {
                parameter.updateExpression(i, k);
                // only count each gene once, not numExpressionCategories times
                //#pragma omp critical
                if(mixtureAssignmentOfGene == k) logLikelihood += propLogLike;
            }else{
            	// only count each gene once, not numExpressionCategories times
                //#pragma omp critical
            	if(mixtureAssignmentOfGene == k) logLikelihood += currLogLike;
            }
        }

        // adjust the the unscaled probabilities by the constant c
        // ln(f') = ln(c) + ln(f)
        // calculate ln(P) = ln( Sum(p_i*f'(...)) ) and obtain normalizing constant for new p_i
        double normalizingProbabilityConstant = 0.0;
        double* probabilities = new double[numMixtures]();
        for(unsigned k = 0u; k < numMixtures; k++)
        {
            unscaledLogProb_curr_singleMixture[k] -= maxValue;
            //TODO compare log vs non log calculation!
            probabilities[k] = std::log(parameter.getCategoryProbability(k)) + unscaledLogProb_curr_singleMixture[k];
            probabilities[k] = std::exp(probabilities[k]);
            normalizingProbabilityConstant += probabilities[k];
        }
        // normalize probabilities
        for (unsigned k = 0u; k < numMixtures; k++)
        {
            probabilities[k] = probabilities[k] / normalizingProbabilityConstant;
        }
        // Get category in which the gene is placed in.
        // If we use multiple sequence observation (like different mutatnts) randMultinom needs an parameter N to place N observations in numMixture buckets
        unsigned categoryOfGene = ROCParameter::randMultinom(probabilities, numMixtures);
        parameter.setMixtureAssignment(i, categoryOfGene);
        dirichletParameters[categoryOfGene] += 1;
        if((iteration % thining) == 0)
        {
            parameter.updateExpressionTrace(iteration/thining, i);
            parameter.updateMixtureAssignmentTrace(iteration/thining, i);
        }
			delete [] probabilities;
    		delete [] unscaledLogProb_curr_singleMixture;
			delete [] unscaledLogProb_prop;
			delete [] unscaledLogProb_curr;
		}

    double* newMixtureProbabilities = new double[numMixtures]();
    ROCParameter::randDirichlet(dirichletParameters, numMixtures, newMixtureProbabilities);
		for(unsigned k = 0u; k < numMixtures; k++)
    {
      parameter.setCategoryProbability(k, newMixtureProbabilities[k]);
    }
    if((iteration % thining) == 0)
    {
        parameter.updateMixtureProbabilitiesTrace(iteration/thining);
    }
    delete [] dirichletParameters;
	delete [] newMixtureProbabilities;
    return logLikelihood;
}

void MCMCAlgorithm::acceptRejectHyperParameter(int numGenes, ROCParameter& parameter, ROCModel& model, int iteration)
{
    double logProbabilityRatio = 0.0;

    double currentSphi = parameter.getSphi(false);
    double currentMPhi = -(currentSphi * currentSphi) / 2;

    double proposedSphi = parameter.getSphi(true);
    double proposedMPhi = -(proposedSphi * proposedSphi) / 2;

    for(int i = 0; i < numGenes; i++)
    {
        unsigned mixture = parameter.getMixtureAssignment(i);
        mixture = parameter.getExpressionCategory(mixture);
        double phi = parameter.getExpression(i, mixture, false);
        logProbabilityRatio += std::log(ROCParameter::densityLogNorm(phi, proposedMPhi, proposedSphi)) - std::log(ROCParameter::densityLogNorm(phi, currentMPhi, currentSphi));
    }

    // Take reverse jumb probability into account if phi proposal width is not identical.
    // If phi proposal width is identical, the term cancels and does not have to be calculated.
    double curr_std_sphi = parameter.getCurrentSphiProposalWidth();
    double prev_std_sphi = parameter.getPreviousSphiProposalWidth();
	double revJump_proposed = 0.0;
	double revJump = 0.0;
	if(curr_std_sphi != prev_std_sphi)
	{
		revJump_proposed = std::log(ROCParameter::densityNorm(proposedSphi, currentSphi, prev_std_sphi));
		revJump = std::log(ROCParameter::densityNorm(currentSphi, proposedSphi, curr_std_sphi));
	}
    logProbabilityRatio -= (std::log(currentSphi) - std::log(proposedSphi)) + (revJump_proposed - revJump);
    if(!std::isfinite(logProbabilityRatio))
    {
    	std::cout << "logProbabilityRatio not finite!\n";
    }

    if( -ROCParameter::randExp(1) < logProbabilityRatio )
    {
        // moves proposed Sphi to current Sphi
        parameter.updateSphi();
    }
    if((iteration % thining) == 0)
    {
        parameter.updateSphiTrace(iteration/thining);
    }
}

void MCMCAlgorithm::acceptRejectCodonSpecificParameter(Genome& genome, ROCParameter& parameter, ROCModel& model, int iteration)
{
    double acceptanceRatioForAllMixtures = 0.0;
    for(unsigned i = 0; i < 22; i++)
    {
    	char curAA = SequenceSummary::AminoAcidArray[i];
        // skip amino acids with only one codon or stop codons
        if(curAA == 'X' || curAA == 'M' || curAA == 'W') continue;
        // calculate likelihood ratio for every Category for current AA
        model.calculateLogLikelihoodRatioPerAAPerCategory(curAA, genome, parameter, acceptanceRatioForAllMixtures);
        if( -ROCParameter::randExp(1) < acceptanceRatioForAllMixtures )
        {
            // moves proposed codon specific parameters to current codon specific parameters
            parameter.updateCodonSpecificParameter(curAA);
		}
        if((iteration % thining) == 0)
        {
        	parameter.updateCodonSpecificParameterTrace(iteration/thining, curAA);
        }
    }
}

void MCMCAlgorithm::run(Genome& genome, ROCModel& model, ROCParameter& parameter)
{
    unsigned maximumIterations = samples * thining;
    // initialize everything

    parameter.initAllTraces(samples, genome.getGenomeSize(), maximumIterations/adaptiveWidth); 
    // starting the MCMC

    std::cout << "entering MCMC loop" << std::endl;
    std::cout << "\tEstimate Codon Specific Parameters? " << (estimateCodonSpecificParameter ? "TRUE" : "FALSE") << std::endl;
    std::cout << "\tEstimate Hyper Parameters? " << (estimateHyperParameter ? "TRUE" : "FALSE") << std::endl;
    std::cout << "\tEstimate Expression Parameters? " << (estimateExpression ? "TRUE" : "FALSE") << std::endl;


    std::cout << "\tStarting MCMC with " << maximumIterations << " iterations\n";
    std::cout << "\tAdaptive width is  " << adaptiveWidth << " iterations\n";


    for(unsigned iteration = 0; iteration < maximumIterations; iteration++)
    {

        if( (iteration + 1) % 100 == 0)
        {
            std::cout << "Status at iteration " << (iteration+1) << std::endl;
            std::cout << "\t current logLikelihood: " << likelihoodTrace[(iteration/thining) - 1] << std::endl;
            std::cout << "\t current Sphi estimate: " << parameter.getSphi() << std::endl;
            std::cout << "\t current Sphi proposal width: " << parameter.getSphiProposalWidth() << std::endl;
            for(unsigned i = 0u; i < parameter.getNumMixtureElements(); i++)
            {
                std::cout << "\t current Mixture element probability for element " << i << ": " << parameter.getCategoryProbability(i) << std::endl;
            }

        }
        if(estimateCodonSpecificParameter) //should the "is" functions be used here instead?
        {
            parameter.proposeCodonSpecificParameter();
            acceptRejectCodonSpecificParameter(genome, parameter, model, iteration);
            if( ( (iteration + 1u) % adaptiveWidth) == 0u)
            {
                parameter.adaptCodonSpecificParameterProposalWidth(adaptiveWidth);
            }
        }
        // update hyper parameter
        if(estimateHyperParameter)
        {
            parameter.proposeSPhi();
            acceptRejectHyperParameter(genome.getGenomeSize(), parameter, model, iteration);
            if( ( (iteration + 1u) % adaptiveWidth) == 0u)
            {
               parameter.adaptSphiProposalWidth(adaptiveWidth);
            }
        }
        // update expression level values
        if(estimateExpression)
        {
            parameter.proposeExpressionLevels();
            double logLike = acceptRejectExpressionLevelForAllGenes(genome, parameter, model, iteration);
            if((iteration % thining) == 0)
            {
                likelihoodTrace[iteration/thining] = logLike;
            }
            if( ( (iteration + 1u) % adaptiveWidth) == 0u)
            {
               parameter.adaptExpressionProposalWidth(adaptiveWidth);
            }
        }
    } // end MCMC loop
    std::cout << "leaving MCMC loop" << std::endl;

//NOTE: The following files used to be written here:
		//selectionParamTrace_#.csv
		//phiTrace_nmix_#.csv
		//mutationParamTrace_#.csv
    //liklihoodTrace.csv
		//sphiTrace.csv
		//expressionLevelTrace.csv

}

double MCMCAlgorithm::getLogLikelihoodPosteriorMean(unsigned samples)
{
	double posteriorMean = 0.0;
	unsigned traceLength = likelihoodTrace.size();


	if(samples > traceLength)
	{
		std::cerr << "Warning in MCMCAlgorithm::getLogLikelihoodPosteriorMean throws: Number of anticipated samples (" <<
			samples << ") is greater than the length of the available trace (" << traceLength << ")." << "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}
	unsigned start = traceLength - samples;
	for(unsigned i = start; i < traceLength; i++)
	{
		posteriorMean += likelihoodTrace[i];
	}

	return posteriorMean / (double)samples;
}
// ---------------------------------------------------------------------------
// ----------------------------- RCPP STUFF ----------------------------------
// ---------------------------------------------------------------------------
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;

RCPP_EXPOSED_CLASS(Genome)
RCPP_EXPOSED_CLASS(ROCParameter)
RCPP_EXPOSED_CLASS(ROCModel)

RCPP_MODULE(MCMCAlgorithm_mod)
{
	class_<MCMCAlgorithm>( "MCMCAlgorithm" )
		.constructor("empty constructor")
		.constructor <int, int, int, bool, bool, bool>()
		.method("run", &MCMCAlgorithm::run)
		.method("getLogLikelihoodTrace", &MCMCAlgorithm::getLogLikelihoodTrace)
		.method("getLogLikelihoodPosteriorMean", &MCMCAlgorithm::getLogLikelihoodPosteriorMean)
	;
}
#endif
