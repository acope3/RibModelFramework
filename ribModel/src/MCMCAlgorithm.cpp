#include "include/MCMCAlgorithm.h"
#include "include/CovarianceMatrix.h"

#ifdef STANDALONE
#include <random>
#endif

#include <cstdlib>
#include <sstream>
#include <chrono>
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
	lastConvergenceTest = 0u;

	estimateMixtureAssignment = true;
}

MCMCAlgorithm::MCMCAlgorithm(unsigned _samples, unsigned _thining, unsigned _adaptiveWidth, bool _estimateSynthesisRate, bool _estimateCodonSpecificParameter, bool _estimateHyperParameter)
: samples(_samples), thining(_thining), adaptiveWidth(_adaptiveWidth * thining), estimateSynthesisRate(_estimateSynthesisRate), estimateCodonSpecificParameter(_estimateCodonSpecificParameter),
	estimateHyperParameter(_estimateHyperParameter)
{
	likelihoodTrace.resize(samples);
	writeRestartFile = false;
	multipleFiles = false;
	fileWriteInterval = 1u;
	lastConvergenceTest = 0u;

	estimateMixtureAssignment = true;
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

	unsigned numSynthesisRateCategories = model.getNumSynthesisRateCategories();
	unsigned numMixtures = model.getNumMixtureElements();
	double* dirichletParameters = new double[numMixtures]();

	
	for (unsigned i = 0u; i < numMixtures; i++) {
		dirichletParameters[i] = 0.0;
	}

	//initialize parameter's size
	for(int i = 0; i < numGenes; i++)
	{
		Gene gene = genome.getGene(i);

		/*
			 Since some values returned by calculateLogLikelihoodRatioPerGene are veyr small (~ -1100), exponentiation leads to 0.
			 To solve this problem, we adjust the value by a constant c. I choose to use the average value across all mixtures.
			 We justify this by
			 P = Sum(p_i*f(...))
			 => f' = c*f
			 => ln(f') = ln(c) + ln(f)
			 => ln(P) = ln( Sum(p_i*f'(...)) )
			 => ln(P) = ln(P') - ln(c)
			 Note that we use invere sign because our values of ln(f) and ln(f') are negative.
		 */        
		
		double maxValue = -1000000.0;
		unsigned mixtureIndex = 0u;
		
		double* unscaledLogProb_curr = new double[numSynthesisRateCategories]();
		double* unscaledLogProb_prop = new double[numSynthesisRateCategories]();

		double* unscaledLogProb_curr_singleMixture = new double[numMixtures]();
		double* probabilities = new double[numMixtures]();

		for (unsigned j = 0u; j < numMixtures; j++) 
		{
			probabilities[j] = 0.0;
			unscaledLogProb_curr_singleMixture[j] = 0.0;
		}

		for (unsigned j = 0u; j < numSynthesisRateCategories; j++)
		{
			unscaledLogProb_curr[j] = 0.0;
			unscaledLogProb_prop[j] = 0.0;
		}
		
		for(unsigned k = 0u; k < numSynthesisRateCategories; k++)
		{
			// logProbabilityRatio contains the logProbabilityRatio in element 0,
			// the current unscaled probability in element 1 and the proposed unscaled probability in element 2
			std::vector<unsigned> mixtureElements = model.getMixtureElementsOfSelectionCategory(k);
			for(unsigned n = 0u; n < mixtureElements.size(); n++)
			{
				unsigned mixtureElement = mixtureElements[n];
				double logProbabilityRatio[3];
				model.calculateLogLikelihoodRatioPerGene(gene, i, mixtureElement, logProbabilityRatio);

/*
				if (std::isinf(logProbabilityRatio[1])) {
					std::cout << "logprob1 inf\n";
				}
				if (std::isinf(logProbabilityRatio[2])) {
					std::cout << "logprob2 inf\n";
				}
*/
				unscaledLogProb_curr[k] += logProbabilityRatio[1];
				unscaledLogProb_prop[k] += logProbabilityRatio[2];

				unscaledLogProb_curr_singleMixture[mixtureIndex] = logProbabilityRatio[1];
				maxValue = unscaledLogProb_curr_singleMixture[mixtureIndex] > maxValue ? unscaledLogProb_curr_singleMixture[mixtureIndex] : maxValue;
				mixtureIndex++;
			}
		}

		unsigned mixAssign = model.getMixtureAssignment(i);
		unsigned geneSynthCat = model.getSynthesisRateCategory(mixAssign);

		// adjust the the unscaled probabilities by the constant c
		// ln(f') = ln(c) + ln(f)
		// calculate ln(P) = ln( Sum(p_i*f'(...)) ) and obtain normalizing constant for new p_i
		double normalizingProbabilityConstant = 0.0;

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

//		unsigned mixtureAssignmentOfGene = model.getMixtureAssignment(i);
		for(unsigned k = 0u; k < numSynthesisRateCategories; k++)
		{
			// We do not need to add std::log(model.getCategoryProbability(k)) since it will cancel in the ratio!
			double currLogLike = unscaledLogProb_curr[k];
			double propLogLike = unscaledLogProb_prop[k];
			if( -Parameter::randExp(1) < (propLogLike - currLogLike) )
			{
				model.updateSynthesisRate(i, k);
				// only count each gene once, not numSynthesigeneIndexsRateCategories times
				//if(geneSynthCat == k)
					//logLikelihood += std::log(model.getCategoryProbability(k)) + propLogLike;
				logLikelihood += probabilities[k] * propLogLike;
			}else{
				// only count each gene once, not numSynthesisRateCategories times
				//if(geneSynthCat == k)
					//logLikelihood += std::log(model.getCategoryProbability(k)) + currLogLike;
				logLikelihood += probabilities[k] * currLogLike;
			}
		}

		if (std::isinf(logLikelihood))
		{
			std::cout << "\tInfinity reached (Gene: " << i << ")\n";
		}

		// Get category in which the gene is placed in.
		// If we use multiple sequence observation (like different mutants) randMultinom needs an parameter N to place N observations in numMixture buckets
		unsigned categoryOfGene = Parameter::randMultinom(probabilities, numMixtures);
		if(estimateMixtureAssignment)
		{
			model.setMixtureAssignment(i, categoryOfGene);
		}

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

	// take all priors into account
	logLikelihood += model.calculateAllPriors();

	double *newMixtureProbabilities = new double[numMixtures]();
	Parameter::randDirichlet(dirichletParameters, numMixtures, newMixtureProbabilities);
	for(unsigned k = 0u; k < numMixtures; k++)
	{
		model.setCategoryProbability(k, newMixtureProbabilities[k]);
	}
	if((iteration % thining) == 0)
	{
		model.updateMixtureProbabilitiesTrace(iteration/thining);
	}
	//delete [] dirichletParameters;
//	delete [] newMixtureProbabilities;
	return logLikelihood;
}

void MCMCAlgorithm::acceptRejectHyperParameter(Genome &genome, Model& model, int iteration)
{
	std::vector <double> logProbabilityRatios;

	model.calculateLogLikelihoodRatioForHyperParameters(genome, iteration, logProbabilityRatios);

	for (unsigned i = 0; i < logProbabilityRatios.size(); i++)
	{
		if (!std::isfinite(logProbabilityRatios[i]))
		{
			std::cout << "logProbabilityRatio " << i << " not finite!\n";
		}

		if (-Parameter::randExp(1) < logProbabilityRatios[i])
		{
			model.updateHyperParameter(i);
		}
	}

	if((iteration % thining) == 0)
	{
		model.updateHyperParameterTraces(iteration/thining);
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

// Allows to diverge from initial conditions (divergenceIterations controls the divergence).
// This allows for varying initial conditions for better exploration of the parameter space.
// This functions does not take the likelihood into account, the random walk is only constrained by the prior
// for each parameter. If no prior (flat prior) is present, an unconstrained random walk is performed.
void MCMCAlgorithm::varyInitialConditions(Genome& genome, Model& model, unsigned divergenceIterations)
{
	// TODO THIS FUNCTON THROWS THE ACCEPTANCE COUNTER OFF!


	// NOTE: IF PRIORS ARE ADDED, TAKE INTO ACCOUNT HERE!
	std::cout << "Allowing divergence from initial conditions for " << divergenceIterations << " iterations.\n" << std::endl;
	// divergence from initial conditions is not stored in trace

	// how many steps do you want to walk "away" from the initial conditions
	for(unsigned iteration = 0u; iteration < divergenceIterations; iteration++)
	{
		// propose all parameters
		model.proposeCodonSpecificParameter();
		model.proposeHyperParameters();
		model.proposeSynthesisRateLevels();

		// no prior on codon specific parameters -> just accept everything
		unsigned size = model.getGroupListSize();
		for(unsigned i = 0; i < size; i++)
		{
			std::string grouping = model.getGrouping(i);
			model.updateCodonSpecificParameter(grouping);
		}

		// no prior on hyper parameters -> just accept everything
		model.updateAllHyperParameter();

		// prior on phi values -> take prior into account, but only the prior no likelihood
		int numGenes = genome.getGenomeSize();
		unsigned numSynthesisRateCategories = model.getNumSynthesisRateCategories();
		for(int i = 0; i < numGenes; i++)
		{

			for(unsigned k = 0u; k < numSynthesisRateCategories; k++)
			{
				// map from mixture to category and obtain corresponding phi value
				unsigned expressionCategory = model.getSynthesisRateCategory(k);
				double phiValue = model.getSynthesisRate(i, expressionCategory, false);
				double phiValue_proposed = model.getSynthesisRate(i, expressionCategory, true);

				unsigned mixture = model.getMixtureAssignment(k);
				mixture = model.getSynthesisRateCategory(mixture);
				double sPhi = model.getSphi(mixture, false);
				double mPhi = (-(sPhi * sPhi) / 2);

				// accept/ reject based on prior ratio
				double logPhiProbability = Parameter::densityLogNorm(phiValue, mPhi, sPhi, true);
				double logPhiProbability_proposed = Parameter::densityLogNorm(phiValue_proposed, mPhi, sPhi, true);
				if( -Parameter::randExp(1) < (logPhiProbability_proposed - logPhiProbability) )
				{
					model.updateSynthesisRate(i, k);
				}
			}
		}

		// update Gibbs sampled parameter.
		model.updateGibbsSampledHyperParameters(genome);
	}
}

void MCMCAlgorithm::run(Genome& genome, Model& model, unsigned numCores, unsigned divergenceIterations)
{
#ifndef __APPLE__
	omp_set_num_threads(numCores);
#endif

	// Allows to diverge from initial conditions (divergenceIterations controls the divergence).
	// This allows for varying initial conditions for better exploration of the parameter space.
	varyInitialConditions(genome, model, divergenceIterations);

	unsigned maximumIterations = samples * thining;
	// initialize everything

	model.setNumPhiGroupings(genome.getGene(0).getObservedPhiValues().size());
	model.initTraces(samples + 1, genome.getGenomeSize()); //Samples + 2 so we can store the starting and ending values.
	// starting the MCMC

	model.updateTracesWithInitialValues(genome);


	std::cout << "entering MCMC loop" << std::endl;
	std::cout << "\tEstimate Codon Specific Parameters? " << (estimateCodonSpecificParameter ? "TRUE" : "FALSE") << std::endl;
	std::cout << "\tEstimate Hyper Parameters? " << (estimateHyperParameter ? "TRUE" : "FALSE") << std::endl;
	std::cout << "\tEstimate SynthesisRate Parameters? " << (estimateSynthesisRate ? "TRUE" : "FALSE") << std::endl;

	// set the last iteration to the max iterations, this way if the MCMC doesn't exit based on Geweke score, it will use the max iteration for posterior means
	model.setLastIteration(samples);
	std::cout << "\tStarting MCMC with " << maximumIterations << " iterations\n";
	//std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
	for(unsigned iteration = 1u; iteration <= maximumIterations; iteration++)
	{
		if (writeRestartFile)
		{
			if ((iteration) % fileWriteInterval  == 0u)
			{
				std::cout <<"Writing restart file!\n";
				if (multipleFiles)
				{
					std::ostringstream oss;
					oss << (iteration) / thining << file;
					std::string tmp = oss.str();
					model.writeRestartFile(tmp);
				}
				else
				{
					model.writeRestartFile(file);
				}
			}
		}
		if( (iteration) % 100u == 0u)
		{
			std::cout << "Status at iteration " << (iteration) << std::endl;
			std::cout << "\t current logLikelihood: " << likelihoodTrace[(iteration/thining) - 1] << std::endl;
			model.printHyperParameters();
			for(unsigned i = 0u; i < model.getNumMixtureElements(); i++)
			{
				std::cout << "\t current Mixture element probability for element " << i << ": " << model.getCategoryProbability(i) << std::endl;
			}
		}
		if(estimateCodonSpecificParameter) 
		{
			model.proposeCodonSpecificParameter();
			acceptRejectCodonSpecificParameter(genome, model, iteration);
			if( ( (iteration) % adaptiveWidth) == 0u)
			{
				model.adaptCodonSpecificParameterProposalWidth(adaptiveWidth);
			}
		}
		// update hyper parameter
		if(estimateHyperParameter)
		{
			model.updateGibbsSampledHyperParameters(genome);
			model.proposeHyperParameters();
			acceptRejectHyperParameter(genome, model, iteration);
			if( ( (iteration) % adaptiveWidth) == 0u)
			{
				model.adaptHyperParameterProposalWidths(adaptiveWidth);
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
			if( ( (iteration) % adaptiveWidth) == 0u)
			{
				model.adaptSynthesisRateProposalWidth(adaptiveWidth);
			}
		}


		/*if( ( (iteration) % (50*adaptiveWidth)) == 0u)
		{
			double gewekeScore = calculateGewekeScore(iteration/thining);
			std::cout << "##################################################" << "\n";
			std::cout << "Geweke Score after " << iteration << " iterations: " << gewekeScore << "\n";
			std::cout << "##################################################" << "\n";

			if(std::abs(gewekeScore) < 1.96)
			{
				std::cout << "Stopping run based on convergence after " << iteration << " iterations\n" << std::endl;
				// Comment out this break to keep the run from stopping on convergence
				model.setLastIteration(iteration/thining);
				//break;
			}
		}*/
	} // end MCMC loop

	std::cout << "leaving MCMC loop" << std::endl;
	//std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
	//std::chrono::duration<int> time_span = std::chrono::duration_cast<std::chrono::duration<int>> (t2 - t1);
	//std::cout << "The MCMC took " << (unsigned) time_span.count() / 3600 << " hours, " << (unsigned) (time_span.count() / 60) % 60 << 
		//" minutes, and " << (unsigned) time_span.count() % 60 << " seconds." << std::endl;
	//NOTE: The following files used to be written here:
	//selectionParamTrace_#.csv
	//phiTrace_nmix_#.csv
	//mutationParamTrace_#.csv
	//liklihoodTrace.csv
	//sphiTrace.csv
	//expressionLevelTrace.csv

}

double MCMCAlgorithm::calculateGewekeScore(unsigned current_iteration)
{
	double posteriorMean1 = 0.0;
	double posteriorMean2 = 0.0;
	double posteriorVariance1 = 0.0;
	double posteriorVariance2 = 0.0;

	unsigned end1 = std::round( (current_iteration - lastConvergenceTest) * 0.1) + lastConvergenceTest;
	unsigned start2 = std::round(current_iteration - (current_iteration * 0.5));

	double numSamples1 = (double) (end1 - lastConvergenceTest);
	double numSamples2 = (double) std::round(current_iteration * 0.5);

	// calculate mean and and variance of first part of likelihood trace
	for(unsigned i = lastConvergenceTest; i < end1; i++)
	{
		posteriorMean1 += likelihoodTrace[i];
	}
	posteriorMean1 = posteriorMean1 / numSamples1;
	for(unsigned i = lastConvergenceTest; i < end1; i++)
	{
		posteriorVariance1 += (likelihoodTrace[i] - posteriorMean1) * (likelihoodTrace[i] - posteriorMean1);
	}
	posteriorVariance1 = posteriorVariance1 / numSamples1;
	// calculate mean and and variance of last part of likelihood trace
	for(unsigned i = start2; i < current_iteration; i++)
	{
		posteriorMean2 += likelihoodTrace[i];
	}
	posteriorMean2 = posteriorMean2 / numSamples2;
	for(unsigned i = start2; i < current_iteration; i++)
	{
		posteriorVariance2 += (likelihoodTrace[i] - posteriorMean2) * (likelihoodTrace[i] - posteriorMean2);
	}
	posteriorVariance2 = posteriorVariance2 / numSamples2;

	lastConvergenceTest = current_iteration;
	// Geweke score
	return (posteriorMean1 - posteriorMean2) / std::sqrt( ( posteriorVariance1 / numSamples1 ) + ( posteriorVariance2 / numSamples2 ) );
}


std::vector<std::vector<double>> MCMCAlgorithm::solveToeplitzMatrix(int lr, std::vector<double> r, std::vector<double> g)
{
// TODO switch f from 2d index to 1d index
	//choleskiMatrix[i * numVariates + k]
	//      solves Toeplitz matrix equation toep(r)f=g(1+.)
	//		by Levinson's algorithm
	//      a is a workspace of size lr, the number of equations
	std::vector<double> f(lr*lr, 0.0);
	std::vector<double> var(lr, 0.0);
	std::vector<std::vector<double>> returnVec(2);

	double* a = new double[lr]();

	unsigned l1,l2,k;
	double v, d, q, hold;
	v = r[0];
	d = r[1];
	a[0] = 1.0;
	f[0] = g[1]/v;
	q = f[0]*r[1];
	var[0] = (1 - f[0]*f[0])*r[0];

	if (lr == 1) return returnVec;
	for(unsigned l = 1; l < lr; l++)
	{
	  a[l] = -d/v;
	  if (l > 1)
	  {
		l1 = (l - 2)/2;
		l2 = l1 + 1;
		for(unsigned j = 1; j < l2; j++)
		{
		  hold = a[j];
		  k = l - j + 1;
		  a[j] = a[j] + a[l]*a[k];
		  a[k] = a[k] + a[l]*hold;
		}
		if (2*l1 != l - 2) a[l2+1] = a[l2+1]*(1.0 + a[l]);
	  }
	  v = v + a[l]*d;
	  f[l*lr + l] = (g[l+1] - q)/v;
	  for(unsigned j = 0; j < (l-1); j++)
	  {
		  f[l*lr + j] = f[(l-1)*lr +j] + f[l*lr + l]*a[l-j+1];
	  }
	//  estimate the innovations variance
	  var[l] = var[l-1] * (1 - f[l*lr + l]*f[l*lr + l]);
	  if (l == lr) return returnVec;
	  d = 0.0;
	  q = 0.0;
	  for(unsigned i = 0; i < l; i++)
	  {
		k = l-i+2;
		d = d + a[i]*r[k];
		q = q + f[l*lr + i]*r[k];
	  }
	}

	returnVec[0] = f;
	returnVec[1] = var;
	return returnVec;
}


std::vector<double> MCMCAlgorithm::acf(std::vector<double>& x, int nrows, int ncols, int lagmax, bool correlation, bool demean)
{
	if(demean)
	{
		double sum = 0.0;
		for(unsigned i = 0u; i < x.size(); i++) sum += x[i];
		double mean = sum / (double)x.size();
		for(unsigned i = 0u; i < x.size(); i++) x[i] = mean;
	}

	std::vector<double> acf(lagmax, 1.0);
    int d1 = lagmax + 1, d2 = ncols*d1;

    for(int u = 0; u < ncols; u++)
    {
		for(int v = 0; v < ncols; v++)
		{
			for(int lag = 0; lag <= lagmax; lag++)
			{
				double sum = 0.0; int nu = 0;
				for(int i = 0; i < nrows-lag; i++)
				{
					nu++;
					sum += x[i + lag + nrows*u] * x[i + nrows*v];
				}
				acf[lag + d1*u + d2*v] = sum/(nu + lag);
			}
		}
    }
    if(correlation) {
		if(nrows == 1) {
			for(int u = 0; u < ncols; u++)
			{
				acf[0 + d1*u + d2*u] = 1.0;
			}
		} else {
			double *se = new double[ncols]();
			for(int u = 0; u < ncols; u++)
			{
				se[u] = sqrt(acf[0 + d1*u + d2*u]);
			}
			for(int u = 0; u < ncols; u++)
			{
				for(int v = 0; v < ncols; v++)
				{
					for(int lag = 0; lag <= lagmax; lag++) // ensure correlations remain in  [-1,1] :
					{
						double a = acf[lag + d1*u + d2*v] / (se[u]*se[v]);
						acf[lag + d1*u + d2*v] = (a > 1.) ? 1. : ((a < -1.) ? -1. : a);
					}
				}
			}
		}
    }
    return acf;
}


double MCMCAlgorithm::getLogLikelihoodPosteriorMean(unsigned _samples)
{
	double posteriorMean = 0.0;
	unsigned traceLength = likelihoodTrace.size();


	if(_samples > traceLength)
	{
		std::cerr << "Warning in MCMCAlgorithm::getLogLikelihoodPosteriorMean throws: Number of anticipated samples (" <<
			_samples << ") is greater than the length of the available trace (" << traceLength << ")." << "Whole trace is used for posterior estimate! \n";
		_samples = traceLength;
	}
	unsigned start = traceLength - _samples;
	for(unsigned i = start; i < traceLength; i++)
	{
		posteriorMean += likelihoodTrace[i];
	}

	return posteriorMean / (double)_samples;
}

void MCMCAlgorithm::setRestartFileSettings(std::string filename, unsigned interval, bool multiple)
{
	file = filename;
	fileWriteInterval = interval * thining;
	multipleFiles = multiple;
	writeRestartFile = true;
}




unsigned MCMCAlgorithm::getSamples()
{
    return samples;
}


unsigned MCMCAlgorithm::getThining()
{
    return thining;
}


unsigned MCMCAlgorithm::getAdaptiveWidth()
{
    return adaptiveWidth;
}

void MCMCAlgorithm::setSamples(unsigned _samples)
{
    samples = _samples;
}


void MCMCAlgorithm::setThining(unsigned _thining)
{
    thining = _thining;
}


void MCMCAlgorithm::setAdaptiveWidth(unsigned _adaptiveWidth)
{
    adaptiveWidth = _adaptiveWidth;
}

void MCMCAlgorithm::setLogLikelihoodTrace(std::vector<double> _likelihoodTrace)
{
    likelihoodTrace = _likelihoodTrace;
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
		.method("setEstimateMixtureAssignment", &MCMCAlgorithm::setEstimateMixtureAssignment)
        .method("getSamples", &MCMCAlgorithm::getSamples)
        .method("getThining", &MCMCAlgorithm::getThining)
        .method("getAdaptiveWidth", &MCMCAlgorithm::getAdaptiveWidth)
        .method("setSamples", &MCMCAlgorithm::setSamples)
        .method("setThining", &MCMCAlgorithm::setThining)
        .method("setAdaptiveWidth", &MCMCAlgorithm::setAdaptiveWidth)
        .method("setLogLikelihoodTrace", &MCMCAlgorithm::setLogLikelihoodTrace)
		;

	Rcpp::function("TestACF", &MCMCAlgorithm::acf); //TEST THAT ONLY!
	Rcpp::function("solveToeplitzMatrix", &MCMCAlgorithm::solveToeplitzMatrix); //TEST THAT ONLY!
}
#endif
