#include "../include/MCMCAlgorithm.h"

#include <random>
#include <cstdlib>
#include <thread>

#include <iostream>
#include <fstream>


#include <omp.h>



MCMCAlgorithm::MCMCAlgorithm()
{
    MCMCAlgorithm(1000, 1, true, true, true);
    adaptiveWidth = 100;
}

MCMCAlgorithm::MCMCAlgorithm(int _samples, int _thining, bool _estimateExpression, bool _estimateCodonSpecificParameter, bool _estimateHyperParameter)
    : samples(_samples), thining(_thining), estimateExpression(_estimateExpression), estimateCodonSpecificParameter(_estimateCodonSpecificParameter),
        estimateHyperParameter(_estimateHyperParameter)
{
    likelihoodTrace.resize(samples + 1);
    adaptiveWidth = 100;
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
    unsigned numMixtures = parameter.getNumMixtureElements();

		//std::cout <<"Entered acceptRejectExpressionLevelForAllGenes\n";
		//std::cout <<"numGenes = " << numGenes <<"\n";
		//std::cout <<"numMixtures = " << numMixtures <<"\n";

    double dirichletParameters[numMixtures];
		//initialize parameter's size
    for(int i = 0; i < numGenes; i++)
    {
        Gene gene = genome.getGene(i);
        double currLike = 0.0;
        double propLike = 0.0;

        double probabilities[numMixtures];

        for(unsigned k = 0; k < numMixtures; k++)
        {
            double* logProbabilityRatio = model.calculateLogLiklihoodRatioPerGene(gene, i, parameter, k);
            //double check probabilites
            probabilities[k] = parameter.getCategoryProbability(k, iteration/thining) * std::exp(logProbabilityRatio[1]);
            currLike += probabilities[k];
            propLike += parameter.getCategoryProbability(k, i) * std::exp(logProbabilityRatio[2]);
        }

        //	std::cout <<"out of initial for loop in acceptRejectExpressionLevelForAllGenes\n";

        // Get category in which the gene is placed in.
        // If we use multiple sequence observrvation (like different mutatnts) randMultinom needs an parameter N to place N observations in numMixture buckets
        unsigned categoryOfGene = ROCParameter::randMultinom(probabilities, numMixtures);
        parameter.setMixtureAssignment(i, categoryOfGene);
        dirichletParameters[categoryOfGene] += 1;

        // accept/reject proposed phi values
        if( ( (double)std::rand() / (double)RAND_MAX ) < (propLike / currLike) )
        {
            // moves proposed phi to current phi
            //std::cout << "Update expression for Gene i = " << i << " in iteration " << iteration << std::endl;
            parameter.updateExpression(i);
            //#pragma omp atomic
            logLikelihood += propLike;
        }else{
            //#pragma omp atomic
            logLikelihood += currLike;
        }
        if((iteration % thining) == 0)
        {
            parameter.updateExpressionTrace(iteration/thining, i);
        }
    }
    double newMixtureProbabilities[numMixtures];
    ROCParameter::randDirichlet(dirichletParameters, numMixtures, newMixtureProbabilities);
    for(unsigned k = 0u; k < numMixtures; k++)
    {
      parameter.setCategoryProbability(k, iteration/thining, newMixtureProbabilities[k]);
    }
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
        double phi = parameter.getExpression(i, false);
        logProbabilityRatio += std::log(ROCParameter::densityLogNorm(phi, proposedMPhi, proposedSphi)) - std::log(ROCParameter::densityLogNorm(phi, currentMPhi, currentSphi));
    }

    logProbabilityRatio -= (std::log(currentSphi) - std::log(proposedSphi));


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
    unsigned numMutationCategories = parameter.getNumMutationCategories();
    unsigned numSelectionCategories = parameter.getNumSelectionCategories();
    unsigned totalNumCategories = numMutationCategories + numSelectionCategories;
    double logAcceptanceRatioPerCategory[totalNumCategories];

    for(unsigned i = 0; i < 22; i++)
    {
        char curAA = SequenceSummary::AminoAcidArray[i];
        // skip amino acids with only one codon or stop codons
        if(curAA == 'X' || curAA == 'M' || curAA == 'W') continue;
        // calculate likelihood ratio for every Category for current AA
        model.calculateLogLikelihoodRatioPerAAPerCategory(curAA, genome, parameter, 1, logAcceptanceRatioPerCategory);

        for(unsigned i = 0; i <  numMutationCategories; i++)
        {
            if( -ROCParameter::randExp(1) < logAcceptanceRatioPerCategory[i] )
            {
                // moves proposed Sphi to current Sphi
                parameter.updateCodonSpecificParameter(i, ROCParameter::dM, curAA);
            }
            if((iteration % thining) == 0)
            {
                //parameter.updateSphiTrace(iteration/thining);
            }
        }
        for(unsigned i = numMutationCategories; i < totalNumCategories; i++)
        {
            if( -ROCParameter::randExp(1) < logAcceptanceRatioPerCategory[i] )
            {
                // moves proposed Sphi to current Sphi
                parameter.updateCodonSpecificParameter(i - numMutationCategories, ROCParameter::dEta, curAA);
            }
            if((iteration % thining) == 0)
            {
                //parameter.updateSphiTrace(iteration/thining);
            }
        }
    }
}

void MCMCAlgorithm::run(Genome& genome, ROCModel& model, ROCParameter& parameter)
{
    unsigned maximumIterations = samples * thining;
    // initialize everything

    parameter.initAllTraces(samples, genome.getGenomeSize());
		parameter.resizeCategoryProbabilites(samples);

    // starting the MCMC
    std::cout << "entering MCMC loop" << std::endl;
		std::cout <<maximumIterations <<"\n";
    for(unsigned iteration = 0; iteration < maximumIterations; iteration++)
    {

        if(iteration % 100 == 0) {std::cout << iteration << std::endl;}
        //std::cout << ".";
        // update codon specific parameter
        if(estimateCodonSpecificParameter) //should the "is" functions be used here instead?
        {
					//	std::cout <<"Doing codon specific parameter update\n";
            parameter.proposeCodonSpecificParameter();
            acceptRejectCodonSpecificParameter(genome, parameter, model, iteration);
        }
        // update hyper parameter
        if(estimateHyperParameter)
        {
						//std::cout <<"Doing hyper parameter update\n";
            parameter.proposeSPhi();
            acceptRejectHyperParameter(genome.getGenomeSize(), parameter, model, iteration);
            if( ( (iteration + 1) % adaptiveWidth) == 0)
            {
               parameter.adaptSphiProposalWidth(adaptiveWidth);
            }
        }
        // update expression level values
        if(estimateExpression)
        {
						//std::cout <<"Doing expression level update\n";
            parameter.proposeExpressionLevels();
            double logLike = acceptRejectExpressionLevelForAllGenes(genome, parameter, model, iteration);
            if((iteration % thining) == 0)
            {
							//	std::cout <<"Thining\n";
                likelihoodTrace[iteration/thining] = logLike;
            }
            if((iteration % adaptiveWidth) == 0)
            {
								//std::cout <<"would call adaptExpressionPro.....\n";
               //parameter.adaptExpressionProposalWidth(adaptiveWidth);
            }
						//std::cout <<"finished expression level update\n";
        }
    } // end MCMC loop
    std::cout << "leaving MCMC loop" << std::endl;

    // development output
    std::ofstream likout("../results/test.lik", std::ofstream::out);
    likout.close();
    std::ofstream sphiout("../results/test.sphi", std::ofstream::out);
    std::ofstream phitraceout("../results/test.phiTrace", std::ofstream::out);
    std::vector<double> sphiTrace = parameter.getSPhiTrace();
    std::vector<std::vector<std::vector<double>>> expressionTrace = parameter.getExpressionTrace();

		for(unsigned iteration = 0; iteration < samples; iteration++)
    {
        likout << likelihoodTrace[iteration] << std::endl;
        sphiout << sphiTrace[iteration] << std::endl;
        for(int i = 0; i < genome.getGenomeSize(); i++)
        {
            phitraceout << expressionTrace[0][iteration][i] << ",";
        }
        phitraceout << std::endl;
    }
    likout.close();
    sphiout.close();
    phitraceout.close();
}

