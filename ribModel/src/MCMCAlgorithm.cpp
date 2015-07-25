#include "include/MCMCAlgorithm.h"
#include "include/CovarianceMatrix.h"

#ifdef STANDALONE
#include <random>
#endif

#include <cstdlib>
#include <sstream>

#include <iostream>
#include <fstream>
#include <stdlib.h> //can be removed later

#ifndef __APPLE__
#include <omp.h>
#include <thread>
#endif



MCMCAlgorithm::MCMCAlgorithm() : samples(1000), thining(1), adaptiveWidth(100 * thining), estimateSynthesisRate(true),
	estimateCodonSpecificParameter(true), estimateHyperParameter(true)
{
	MCMCAlgorithm(1000, 1, true, true, true);
	likelihoodTrace.resize(samples);
	writeRestartFile = false;
	multipleFiles = false;
	fileWriteInterval = 1u;
	numCores = 1u;
}

	MCMCAlgorithm::MCMCAlgorithm(unsigned _samples, unsigned _thining, unsigned _adaptiveWidth, bool _estimateSynthesisRate, bool _estimateCodonSpecificParameter, bool _estimateHyperParameter)
: samples(_samples), thining(_thining), adaptiveWidth(_adaptiveWidth * thining), estimateSynthesisRate(_estimateSynthesisRate), estimateCodonSpecificParameter(_estimateCodonSpecificParameter),
	estimateHyperParameter(_estimateHyperParameter)
{
	likelihoodTrace.resize(samples);
	writeRestartFile = false;
	multipleFiles = false;
	fileWriteInterval = 1u;
	numCores = 1u;
}

MCMCAlgorithm::~MCMCAlgorithm()
{
	//dtor
}

double MCMCAlgorithm::acceptRejectSynthesisRateLevelForAllGenes(Genome& genome, Model& model, int iteration)
{
	// TODO move the likelihood calculation out off here. make it a void function again.

	double logLikelihood = 0.0;
	int numGenes = genome.getGenomeSize();

	// just for testing
	//unsigned concurrentThreadsSupported = std::thread::hardware_concurrency();
	//omp_set_num_threads(concurrentThreadsSupported);
	
	// testing end
	unsigned numSynthesisRateCategories = model.getNumSynthesisRateCategories();
	unsigned numMixtures = model.getNumMixtureElements();
	double* dirichletParameters = new double[numMixtures]();
	//initialize parameter's size
	#pragma omp parallel for //shared(parameter)
	for(int i = 0; i < numGenes; i++)
	{
		Gene gene = genome.getGene(i);

		/*
			 Since some values returned by calculateLogLikelihoodRatioPerGene are veyr small (~ -1100), exponantiation leads to 0.
			 To solve this problem, we adjust the value by a constant c. I choose to use the average value across all mixtures.
			 We justify this by
			 P = Sum(p_i*f(...))
			 => f' = c*f
			 => ln(f') = ln(c) + ln(f)
			 => ln(P) = ln( Sum(p_i*f'(...)) )
			 => ln(P) = ln(P') - ln(c)
			 Note that we use invere sign because our values of ln(f) and ln(f') are negative.
		 */        
		double* unscaledLogProb_curr = new double[numSynthesisRateCategories]();
		double* unscaledLogProb_prop = new double[numSynthesisRateCategories]();

		double* unscaledLogProb_curr_singleMixture = new double[numMixtures]();
		double maxValue = -1000000.0;
		unsigned mixtureIndex = 0u;
		for(unsigned k = 0u; k < numSynthesisRateCategories; k++)
		{
			// logProbabilityRatio contains the logProbabilityRatio in element 0,
			// the current unscaled probability in elemen 1 and the proposed unscaled probability in element 2
			std::vector<unsigned> mixtureElements = model.getMixtureElementsOfSelectionCategory(k);
			for(unsigned n = 0u; n < mixtureElements.size(); n++)
			{
				unsigned mixtureElement = mixtureElements[n];
				double logProbabilityRatio[3];
				model.calculateLogLikelihoodRatioPerGene(gene, i, mixtureElement, logProbabilityRatio);

				// store values so they can be processed
				unscaledLogProb_curr[k] += logProbabilityRatio[1];
				unscaledLogProb_prop[k] += logProbabilityRatio[2];

				unscaledLogProb_curr_singleMixture[mixtureIndex] = logProbabilityRatio[1];
				maxValue = unscaledLogProb_curr_singleMixture[mixtureIndex] > maxValue ? unscaledLogProb_curr_singleMixture[mixtureIndex] : maxValue;
				mixtureIndex++;
			}
		}
		unsigned mixtureAssignmentOfGene = model.getMixtureAssignment(i);
		for(unsigned k = 0u; k < numSynthesisRateCategories; k++)
		{
			// We do not need to add std::log(model.getCategoryProbability(k)) since it will cancel in the ratio!
			double currLogLike = unscaledLogProb_curr[k];
			double propLogLike = unscaledLogProb_prop[k];
			if( -Parameter::randExp(1) < (propLogLike - currLogLike) )
			{
				model.updateSynthesisRate(i, k);
				// only count each gene once, not numSynthesisRateCategories times
				#pragma omp critical
				logLikelihood += model.getCategoryProbability(k) * propLogLike;
			}else{
				// only count each gene once, not numSynthesisRateCategories times
				#pragma omp critical
				logLikelihood += model.getCategoryProbability(k) * currLogLike;
			}
		}

		if (std::isinf(logLikelihood)) std::cout <<"\tInfinity reached\n";
		// adjust the the unscaled probabilities by the constant c
		// ln(f') = ln(c) + ln(f)
		// calculate ln(P) = ln( Sum(p_i*f'(...)) ) and obtain normalizing constant for new p_i
		double normalizingProbabilityConstant = 0.0;
		double* probabilities = new double[numMixtures]();
		for(unsigned k = 0u; k < numMixtures; k++)
		{
			unscaledLogProb_curr_singleMixture[k] -= maxValue;
			//TODO compare log vs non log calculation!
			probabilities[k] = std::log(model.getCategoryProbability(k)) + unscaledLogProb_curr_singleMixture[k];
			probabilities[k] = std::exp(probabilities[k]);
			normalizingProbabilityConstant += probabilities[k];
		}
		// normalize probabilities
		for (unsigned k = 0u; k < numMixtures; k++)
		{
			probabilities[k] = probabilities[k] / normalizingProbabilityConstant;
		}
		// Get category in which the gene is placed in.
		// If we use multiple sequence observation (like different mutants) randMultinom needs an parameter N to place N observations in numMixture buckets
		unsigned categoryOfGene = Parameter::randMultinom(probabilities, numMixtures);
		model.setMixtureAssignment(i, categoryOfGene);
		dirichletParameters[categoryOfGene] += 1;
		if((iteration % thining) == 0)
		{
			model.updateSynthesisRateTrace(iteration/thining, i);
			model.updateMixtureAssignmentTrace(iteration/thining, i);
		}
		delete [] probabilities;
		delete [] unscaledLogProb_curr_singleMixture;
		delete [] unscaledLogProb_prop;
		delete [] unscaledLogProb_curr;
	}
	double* newMixtureProbabilities = new double[numMixtures]();
	Parameter::randDirichlet(dirichletParameters, numMixtures, newMixtureProbabilities);
	for(unsigned k = 0u; k < numMixtures; k++)
	{
		model.setCategoryProbability(k, newMixtureProbabilities[k]);
	}
	if((iteration % thining) == 0)
	{
		model.updateMixtureProbabilitiesTrace(iteration/thining);
	}
	delete [] dirichletParameters;
	delete [] newMixtureProbabilities;
	return logLikelihood;
}

void MCMCAlgorithm::acceptRejectHyperParameter(int numGenes, Model& model, int iteration)
{
	double logProbabilityRatio = 0.0;
	double currentSphi = model.getSphi(false);
	double currentMPhi = -(currentSphi * currentSphi) / 2;

	double proposedSphi = model.getSphi(true);
	double proposedMPhi = -(proposedSphi * proposedSphi) / 2;

#ifndef __APPLE__
#pragma omp parallel for reduction(+:logProbabilityRatio)
#endif
	for(int i = 0; i < numGenes; i++)
	{
		unsigned mixture = model.getMixtureAssignment(i);
		mixture = model.getSynthesisRateCategory(mixture);
		double phi = model.getSynthesisRate(i, mixture, false);
		logProbabilityRatio += std::log(Parameter::densityLogNorm(phi, proposedMPhi, proposedSphi)) - std::log(Parameter::densityLogNorm(phi, currentMPhi, currentSphi));
	}

	logProbabilityRatio -= (std::log(currentSphi) - std::log(proposedSphi));
	if(!std::isfinite(logProbabilityRatio))
	{
		std::cout << "logProbabilityRatio not finite!\n";
		std::cout <<"currentSphi = " << currentSphi << " proposedSphi = " << proposedSphi <<"\n";
	}

	if( -Parameter::randExp(1) < logProbabilityRatio )
	{
		// moves proposed Sphi to current Sphi
		model.updateSphi();
	}

	if((iteration % thining) == 0)
	{
		model.updateSphiTrace(iteration/thining);
	}
}

void MCMCAlgorithm::acceptRejectCodonSpecificParameter(Genome& genome, Model& model, int iteration)
{
	double acceptanceRatioForAllMixtures = 0.0;
	unsigned size = model.getGroupListSize();

	for(unsigned i = 0; i < size; i++)
	{
		std::string grouping = model.getGrouping(i);

		// calculate likelihood ratio for every Category for current AA
		model.calculateLogLikelihoodRatioPerGroupingPerCategory(grouping, genome, acceptanceRatioForAllMixtures);
		if( -Parameter::randExp(1) < acceptanceRatioForAllMixtures )
		{
			// moves proposed codon specific parameters to current codon specific parameters
			model.updateCodonSpecificParameter(grouping);
		}
		if((iteration % thining) == 0)
		{
			model.updateCodonSpecificParameterTrace(iteration/thining, grouping);
		}
	}
}

void MCMCAlgorithm::run(Genome& genome, Model& model, unsigned numCores)
{
#ifndef __APPLE__
	omp_set_num_threads(numCores);
#endif

	unsigned maximumIterations = samples * thining;
	// initialize everything

	model.initTraces(samples, genome.getGenomeSize()); 
	// starting the MCMC

	std::cout << "entering MCMC loop" << std::endl;
	std::cout << "\tEstimate Codon Specific Parameters? " << (estimateCodonSpecificParameter ? "TRUE" : "FALSE") << std::endl;
	std::cout << "\tEstimate Hyper Parameters? " << (estimateHyperParameter ? "TRUE" : "FALSE") << std::endl;
	std::cout << "\tEstimate SynthesisRate Parameters? " << (estimateSynthesisRate ? "TRUE" : "FALSE") << std::endl;


	std::cout << "\tStarting MCMC with " << maximumIterations << " iterations\n";
	for(unsigned iteration = 0u; iteration < maximumIterations; iteration++)
	{
		if (writeRestartFile)
		{
			if ((iteration + 1u) % fileWriteInterval  == 0u)
			{
				std::cout <<"Writing restart file!\n";
				if (multipleFiles)
				{
					std::ostringstream oss;
					oss << (iteration + 1) / thining << file;
					std::string tmp = oss.str();
					model.writeRestartFile(tmp);
				}
				else
				{
					model.writeRestartFile(file);
				}
			}
		}
		if( (iteration + 1u) % 100u == 0u)
		{
			std::cout << "Status at iteration " << (iteration+1) << std::endl;
			std::cout << "\t current logLikelihood: " << likelihoodTrace[(iteration/thining) - 1] << std::endl;
			std::cout << "\t current Sphi estimate: " << model.getSphi() << std::endl;
			std::cout << "\t current Sphi proposal width: " << model.getSphiProposalWidth() << std::endl;
			for(unsigned i = 0u; i < model.getNumMixtureElements(); i++)
			{
				std::cout << "\t current Mixture element probability for element " << i << ": " << model.getCategoryProbability(i) << std::endl;
			}
		}
		if(estimateCodonSpecificParameter) 
		{
			model.proposeCodonSpecificParameter();
			acceptRejectCodonSpecificParameter(genome, model, iteration);
			if( ( (iteration + 1u) % adaptiveWidth) == 0u)
			{
				model.adaptCodonSpecificParameterProposalWidth(adaptiveWidth);
			}
		}
		// update hyper parameter
		if(estimateHyperParameter)
		{
			model.proposeSPhi();
			acceptRejectHyperParameter(genome.getGenomeSize(), model, iteration);
			if( ( (iteration + 1u) % adaptiveWidth) == 0u)
			{
				model.adaptSphiProposalWidth(adaptiveWidth);
			}
		}
		// update expression level values
		if(estimateSynthesisRate)
		{
			model.proposeSynthesisRateLevels();
			double logLike = acceptRejectSynthesisRateLevelForAllGenes(genome, model, iteration);
			if((iteration % thining) == 0u)
			{
				likelihoodTrace[iteration/thining] = logLike;
			}
			if( ( (iteration + 1u) % adaptiveWidth) == 0u)
			{
				model.adaptSynthesisRateProposalWidth(adaptiveWidth);
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

void MCMCAlgorithm::setRestartFileSettings(std::string filename, unsigned interval, bool multiple)
{
	file = filename;
	fileWriteInterval = interval * thining;
	multipleFiles = multiple;
	writeRestartFile = true;
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
RCPP_EXPOSED_CLASS(Model)

RCPP_MODULE(MCMCAlgorithm_mod)
{
	class_<MCMCAlgorithm>( "MCMCAlgorithm" )
		.constructor("empty constructor")
		.constructor <unsigned, unsigned, unsigned, bool, bool, bool>()
		.method("run", &MCMCAlgorithm::run)
		.method("getLogLikelihoodTrace", &MCMCAlgorithm::getLogLikelihoodTrace)
		.method("getLogLikelihoodPosteriorMean", &MCMCAlgorithm::getLogLikelihoodPosteriorMean)
		.method("setRestartFileSettings", &MCMCAlgorithm::setRestartFileSettings)
		;
}
#endif
