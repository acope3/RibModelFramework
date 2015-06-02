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



MCMCAlgorithm::MCMCAlgorithm() : samples(1000), thining(1), estimateExpression(true), estimateCodonSpecificParameter(true),
estimateHyperParameter(true)
{
    MCMCAlgorithm(1000, 1, true, true, true);
    likelihoodTrace.resize(samples + 1);
    adaptiveWidth = 100;
}

MCMCAlgorithm::MCMCAlgorithm(int _samples, int _thining, bool _estimateExpression, bool _estimateCodonSpecificParameter, bool _estimateHyperParameter)
    : samples(_samples), thining(_thining), estimateExpression(_estimateExpression), estimateCodonSpecificParameter(_estimateCodonSpecificParameter),
        estimateHyperParameter(_estimateHyperParameter)
{
    // TODO add adaptiveWidth to constructor
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

    unsigned numExpressionCategories = parameter.getNumExpressionCategories();
    unsigned numMixtures = parameter.getNumMixtureElements();
    double dirichletParameters[numMixtures];
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
        double unscaledLogProb_curr[numExpressionCategories];
        double unscaledLogProb_prop[numExpressionCategories];

        double unscaledLogProb_curr_singleMixture[numMixtures];
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
                //#pragma omp atomic
                // only count each gene once, not numExpressionCategories times
                if(mixtureAssignmentOfGene == k) logLikelihood += propLogLike;
            }else{
                //#pragma omp atomic
            	// only count each gene once, not numExpressionCategories times
            	if(mixtureAssignmentOfGene == k) logLikelihood += currLogLike;
            }
        }

        // adjust the the unscaled probabilities by the constant c
        // ln(f') = ln(c) + ln(f)
        // calculate ln(P) = ln( Sum(p_i*f'(...)) ) and obtain normalizing constant for new p_i
        double normalizingProbabilityConstant = 0.0;
        double probabilities[numMixtures];
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
    }

    double newMixtureProbabilities[numMixtures];
    ROCParameter::randDirichlet(dirichletParameters, numMixtures, newMixtureProbabilities);
    for(unsigned k = 0u; k < numMixtures; k++)
    {
      parameter.setCategoryProbability(k, newMixtureProbabilities[k]);
    }
    if((iteration % thining) == 0)
    {
        parameter.updateCategoryProbabilitiesTrace(iteration/thining);
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
        unsigned mixture = parameter.getMixtureAssignment(i);
        mixture = parameter.getExpressionCategory(mixture);
        double phi = parameter.getExpression(i, mixture, false);
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
            //parameter.updateSphiTrace(iteration/thining);
        }
    }
}

void MCMCAlgorithm::run(Genome& genome, ROCModel& model, ROCParameter& parameter)
{
    unsigned maximumIterations = samples * thining;
    // initialize everything

    parameter.initAllTraces(samples, genome.getGenomeSize());
    CovarianceMatrix covmat = CovarianceMatrix(2);
    // starting the MCMC

    std::cout << "entering MCMC loop" << std::endl;
    std::cout << "\tEstimate Codon Specific Parameters? " << (estimateCodonSpecificParameter ? "TRUE" : "FALSE") << std::endl;
    std::cout << "\tEstimate Hyper Parameters? " << (estimateHyperParameter ? "TRUE" : "FALSE") << std::endl;
    std::cout << "\tEstimate Expression Parameters? " << (estimateExpression ? "TRUE" : "FALSE") << std::endl;


    std::cout << "\tStarting MCMC with " << maximumIterations << " iterations\n";
    for(unsigned iteration = 0; iteration < maximumIterations; iteration++)
    {

        if( (iteration + 1) % 100 == 0)
        {
            std::cout << (iteration+1) << std::endl;
            std::cout << "\t current logLikelihood: " << likelihoodTrace[iteration/thining] << std::endl;
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
            if( ( (iteration + 1) % adaptiveWidth) == 0)
            {
                parameter.adaptCodonSpecificParameterProposalWidth(adaptiveWidth);
            }
        }
        // update hyper parameter
        if(estimateHyperParameter)
        {
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
            parameter.proposeExpressionLevels();
            double logLike = acceptRejectExpressionLevelForAllGenes(genome, parameter, model, iteration);
            if((iteration % thining) == 0)
            {
                likelihoodTrace[iteration/thining] = logLike;
            }
            if( ( (iteration + 1) % adaptiveWidth) == 0)
            {
               parameter.adaptExpressionProposalWidth(adaptiveWidth);
            }
        }
    } // end MCMC loop
    std::cout << "leaving MCMC loop" << std::endl;

    // development output
    std::vector<std::vector<std::vector<double>>> expressionTrace = parameter.getExpressionTrace();
    unsigned numParam = parameter.getNumParam();
    for(unsigned nm = 0u; nm < parameter.getNumSelectionCategories(); nm++)
    {
        std::vector<std::vector<std::vector<double>>> selectionParameterTrace = parameter.getSelectionParameterTrace();
        std::stringstream oss;
        std::stringstream oss2;
        oss2 << "results/selectionParamTrace_" << nm<< ".csv";
        oss << "results/phiTrace_nmix_" << nm << ".csv";
        std::ofstream phitraceout(oss.str(), std::ofstream::out);
        std::ofstream selectTraceOut(oss2.str(), std::ofstream::out);
        for (int i = 0; i < 22; i++)
        {
            char aa = SequenceSummary::AminoAcidArray[i];
            if (aa == 'X' || aa == 'W' || aa == 'M') continue;
            unsigned aaRange[2];
            SequenceSummary::AAToCodonRange(aa, false, aaRange);
            for (unsigned j = aaRange[0]; j < aaRange[1] - 1; j++)
            {
                selectTraceOut << aa <<"." << SequenceSummary::codonArray[j] <<",";
            }
        }
        selectTraceOut <<"\n";
        for(unsigned iteration = 0; iteration < samples; iteration++)
        {
            for(unsigned i = 0u; i < genome.getGenomeSize(); i++)
            {
                phitraceout << expressionTrace[nm][iteration][i] << ",";
            }
            phitraceout << std::endl;
            for(unsigned i = 0; i < numParam; i++)
            {
                selectTraceOut << selectionParameterTrace[nm][iteration][i] <<",";
            }
            selectTraceOut << std::endl;
        }
        phitraceout.close();
        selectTraceOut.close();
    }

    for (unsigned nm = 0u; nm < parameter.getNumMutationCategories(); nm++)
    {
        std::vector<std::vector<std::vector<double>>> mutationParameterTrace = parameter.getMutationParameterTrace();
        std::stringstream oss;
        oss << "results/mutationParamTrace_" << nm << ".csv";
        std::ofstream mutateTraceOut(oss.str(), std::ofstream::out);
        for (int i = 0; i < 22; i++)
        {
            char aa = SequenceSummary::AminoAcidArray[i];
            if (aa == 'X' || aa == 'W' || aa == 'M') continue;
            unsigned aaRange[2];
            SequenceSummary::AAToCodonRange(aa, false, aaRange);
            for (unsigned j = aaRange[0]; j < aaRange[1] - 1; j++)
            {
                mutateTraceOut << aa <<"." << SequenceSummary::codonArray[j] <<",";
            }
        }
        mutateTraceOut <<"\n";
        for(unsigned iteration = 0; iteration < samples; iteration++)
        {
            for(unsigned i = 0u; i < numParam; i++)
            {
                mutateTraceOut << mutationParameterTrace[nm][iteration][i] <<",";
            }
            mutateTraceOut << std::endl;
        }
        mutateTraceOut.close();
    }
    std::ofstream likout("results/liklihoodTrace.csv", std::ofstream::out);
    std::ofstream sphiout("results/sphiTrace.csv", std::ofstream::out);
    std::vector<double> sphiTrace = parameter.getSPhiTrace();
    for(unsigned iteration = 0u; iteration < samples; iteration++)
    {
        likout << likelihoodTrace[iteration] << std::endl;
        sphiout << sphiTrace[iteration] << std::endl;
    }
    likout.close();
    sphiout.close();


    std::vector<std::vector<double>> phiTraces(genome.getGenomeSize());
    for(unsigned i = 0u; i < genome.getGenomeSize(); i++)
    {
        phiTraces[i] = parameter.getExpressionTrace(i);
    }
    std::ofstream phitraceout("results/expressionLevelTrace.csv", std::ofstream::out);
    for(unsigned iteration = 0u; iteration < samples; iteration++)
    {
        for(unsigned i = 0u; i < genome.getGenomeSize(); i++)
        {
            phitraceout << phiTraces[i][iteration] << ",";
        }
        phitraceout << std::endl;
    }

}

