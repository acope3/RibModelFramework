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
}

MCMCAlgorithm::MCMCAlgorithm(int _samples, int _thining, bool _estimateExpression, bool _estimateCodonSpecificParameter, bool _estimateHyperParameter)
    : samples(_samples), thining(_thining), estimateExpression(_estimateExpression), estimateCodonSpecificParameter(_estimateCodonSpecificParameter),
        estimateHyperParameter(_estimateHyperParameter)
{
    likelihoodTrace.resize(samples + 1);
}

MCMCAlgorithm::~MCMCAlgorithm()
{
    //dtor
}

MCMCAlgorithm::MCMCAlgorithm(const MCMCAlgorithm& other)
{
    //copy ctor
}

double MCMCAlgorithm::acceptRejectExpressionLevelForAllGenes(Genome& genome, ROCParameter& parameter, ROCModel& model, int iteration)
{
    // TODO move the likelihood calculation out off here. make it a void function again.

    double logLikelihood = 0.0;
    int numGenes = genome.getGenomeSize();

    // just for testing
    //unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
    //omp_set_num_threads(concurentThreadsSupported);
    //#pragma omp parallel for shared(parameter)
    // testing end
    for(int i = 0; i < numGenes; i++)
    {
        Gene gene = genome.getGene(i);
        double currentLogLikelihood = model.calculateLogLiklihoodPerGene(gene, i, parameter, false);
        double proposedLogLikelihood = model.calculateLogLiklihoodPerGene(gene, i, parameter, true);

        if(currentLogLikelihood != currentLogLikelihood)
        {
            double bla = 0.0;
        }
        if(proposedLogLikelihood != proposedLogLikelihood)
        {
            double blub = 0.0;
        }
        // TODO move next line into parameter object
        //double curr = parameter.getExpression(i, false);
        //double propo = parameter.getExpression(i, true);
        //double logImportanceRatio = curr - propo; // we need that if we do a change of variables from phi to log(phi). ONLY IF: log(phi)~N(m,s)... NOT IF phi~logN(m,s)
        // accept/reject proposed phi values
        if( ( (double)std::rand() / (double)RAND_MAX ) < std::exp(proposedLogLikelihood - currentLogLikelihood) )
        {
            // moves proposed phi to current phi
            parameter.updateExpression(i);
            //#pragma omp atomic
            logLikelihood += proposedLogLikelihood;
        }else{
            //#pragma omp atomic
            logLikelihood += currentLogLikelihood;
        }
        if((iteration % thining) == 0)
        {
            parameter.updateExpressionTrace(iteration/thining, i);
        }

    }
    if(logLikelihood != logLikelihood)
    {
        double foo = 0.0;
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

    //double logImportanceRatio = currentSphi - proposedSphi; Sphi is proposed using a logN dist. So thi line is not necessary!!

    for(int i = 0; i < numGenes; i++)
    {
        double phi = parameter.getExpression(i, false);
        double a = std::log(ROCParameter::densityLogNorm(phi, proposedMPhi, proposedSphi));
        double b = std::log(ROCParameter::densityLogNorm(phi, currentMPhi, currentSphi));
        logProbabilityRatio += (a - b);
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


}

void MCMCAlgorithm::run(Genome& genome, ROCModel& model, ROCParameter& parameter)
{
    int maximumIterations = samples * thining;
    // initialize everything
    parameter.initAllTraces(samples, genome.getGenomeSize());
    //parameter.initExpressionTrace(samples, genome.getGenomeSize()); // not necessary if initAllTraces is implemented!

    //double likelihoodTrace[maximumIterations];
    // starting the MCMC
    std::cout << "entering MCMC loop" << std::endl;
    for(int iteration = 0; iteration < maximumIterations; iteration++)
    {

        if(iteration % 100 == 0) {std::cout << iteration << std::endl;}
        std::cout << ".";
        // update codon specific parameter
        if(estimateCodonSpecificParameter)
        {
            //parameter.proposeCodonSpecificParameter();
            //acceptRejectCodonSpecificParameter(genome, parameter, model);
        }
        // update hyper parameter
        if(estimateHyperParameter)
        {
            parameter.proposeSPhi();
            acceptRejectHyperParameter(genome.getGenomeSize(), parameter, model, iteration);
            // TODO here goes the estimation of S_phi and so on...
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

        }

    } // end MCMC loop
    std::cout << "leaving MCMC loop" << std::endl;

    std::ofstream likout("/home/clandere/CodonUsageBias/organisms/yeast/results/test.lik");
    std::ofstream sphiout("/home/clandere/CodonUsageBias/organisms/yeast/results/test.sphi");
    std::vector<double> sphiTrace = parameter.getSPhiTrace();
    for(int iteration = 0; iteration < samples; iteration++)
    {
        likout << likelihoodTrace[iteration] << std::endl;
        sphiout << sphiTrace[iteration] << std::endl;
    }


}


