#include "include/RFP/RFPModel.h"
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

RFPModel::RFPModel() : Model()
{
	//ctor
}


RFPModel::~RFPModel()
{
	//dtor
	//TODO: call Parent's deconstructor
}


double RFPModel::calculateLogLikelihoodPerCodonPerGene(double currAlpha, double currLambdaPrime, unsigned currRFPObserved, 
		unsigned currNumCodonsInMRNA, double phiValue)
{
	double logLikelihood = ((std::lgamma((currNumCodonsInMRNA * currAlpha) + currRFPObserved)) - (std::lgamma(currNumCodonsInMRNA * currAlpha)))
		+ (currRFPObserved * (std::log(phiValue) - std::log(currLambdaPrime + phiValue))) + ((currNumCodonsInMRNA * currAlpha) * (std::log(currLambdaPrime) -
					std::log(currLambdaPrime + phiValue)));

	return logLikelihood;
}


void RFPModel::setParameter(RFPParameter &_parameter)
{
	parameter = &_parameter;
}


void RFPModel::calculateLogLikelihoodRatioPerGene(Gene& gene, int geneIndex, unsigned k, double* logProbabilityRatio)
{
	//TODO: Many functions used here are unimplimented as of now
	double logLikelihood = 0.0;
	double logLikelihood_proposed = 0.0;

	unsigned alphaCategory = parameter->getMutationCategory(k);
	unsigned lambdaPrimeCategory = parameter->getSelectionCategory(k);
	unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(k);


	double phiValue = parameter->getSynthesisRate(geneIndex, synthesisRateCategory, false);
	double phiValue_proposed = parameter->getSynthesisRate(geneIndex, synthesisRateCategory, true);


	for (unsigned index = 0; index < getGroupListSize(); index++) //number of codons, without the stop codons
	{
		std::string codon = getGrouping(index);

		double currAlpha = getParameterForCategory(alphaCategory, RFPParameter::alp, codon, false);
		double currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, RFPParameter::lmPri, codon, false);
		unsigned currRFPObserved = gene.geneData.getRFPObserved(index);

		unsigned currNumCodonsInMRNA = gene.geneData.getCodonCountForCodon(index);
		if (currNumCodonsInMRNA == 0) continue;

		logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, currRFPObserved, currNumCodonsInMRNA, phiValue);
		logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, currRFPObserved, currNumCodonsInMRNA, phiValue_proposed);
	}

	double sPhi = parameter->getSphi(false);
	double logPhiProbability = std::log(Parameter::densityLogNorm(phiValue, (-(sPhi * sPhi) / 2), sPhi));
	double logPhiProbability_proposed = std::log(Parameter::densityLogNorm(phiValue_proposed, (-(sPhi * sPhi) / 2), sPhi));
	double currentLogLikelihood = (logLikelihood + logPhiProbability);
	double proposedLogLikelihood = (logLikelihood_proposed + logPhiProbability_proposed);

	logProbabilityRatio[0] = (proposedLogLikelihood - currentLogLikelihood) - (std::log(phiValue) - std::log(phiValue_proposed));
	logProbabilityRatio[1] = currentLogLikelihood - std::log(phiValue_proposed);
	logProbabilityRatio[2] = proposedLogLikelihood - std::log(phiValue);
}


void RFPModel::calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome, double& logAcceptanceRatioForAllMixtures)
{
	double logLikelihood = 0.0;
	double logLikelihood_proposed = 0.0;
	Gene *gene;
	unsigned index = SequenceSummary::CodonToIndex(grouping);

	Gene *gene;
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		gene = &genome.getGene(i);
		// which mixture element does this gene belong to
		unsigned mixtureElement = parameter->getMixtureAssignment(i);
		// how is the mixture element defined. Which categories make it up
		unsigned alphaCategory = parameter->getMutationCategory(mixtureElement);
		unsigned lambdaPrimeCategory = parameter->getSelectionCategory(mixtureElement);
		unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(mixtureElement);
		// get non codon specific values, calculate likelihood conditional on these
		double phiValue = parameter->getSynthesisRate(i, synthesisRateCategory, false);
		unsigned currRFPObserved = gene->geneData.getRFPObserved(index);
		unsigned currNumCodonsInMRNA = gene->geneData.getCodonCountForCodon(index);
		if (currNumCodonsInMRNA == 0) continue;


		double currAlpha = getParameterForCategory(alphaCategory, RFPParameter::alp, grouping, false);
		double currLambdaPrime = getParameterForCategory(lambdaPrimeCategory, RFPParameter::lmPri, grouping, false);

		double propAlpha = getParameterForCategory(alphaCategory, RFPParameter::alp, grouping, true);
		double propLambdaPrime = getParameterForCategory(lambdaPrimeCategory, RFPParameter::lmPri, grouping, true);

	/*	std::cout <<"currAlpha: " << currAlpha <<"\n";
		std::cout <<"currLambdaPrime: " << currLambdaPrime <<"\n";
		std::cout <<"propAlpha: " << propAlpha <<"\n";
		std::cout <<"propLambdaPrime: " << propLambdaPrime <<"\n";
		std::cout <<"phiValue: " << phiValue <<"\n";
		std::cout <<"currRFPObserved: " << currRFPObserved <<"\n";
		std::cout <<"currNumCodonsInMRNA: " << currNumCodonsInMRNA <<"\n";
		*/logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, currRFPObserved, currNumCodonsInMRNA, phiValue);
		logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(propAlpha, propLambdaPrime, currRFPObserved, currNumCodonsInMRNA, phiValue);
//		std::cout << "logLikelihood: " << logLikelihood <<"\n";
//		std::cout << "logLikelihood_proposed: " << logLikelihood_proposed <<"\n";
	}
	logAcceptanceRatioForAllMixtures = logLikelihood_proposed - logLikelihood;
}


void RFPModel::simulateGenome(Genome &genome)
{
	for (unsigned geneIndex = 0; geneIndex < genome.getGenomeSize(); geneIndex++)
	{
		unsigned mixtureElement = getMixtureAssignment(geneIndex);
		Gene gene = genome.getGene(geneIndex);
		Gene tmpGene = gene;
		for (unsigned codonIndex = 0; codonIndex < 64; codonIndex++)
		{
			std::string codon = SequenceSummary::codonArray[codonIndex];
			unsigned alphaCat = parameter -> getMutationCategory(mixtureElement);
			unsigned lambdaPrimeCat = parameter -> getSelectionCategory(mixtureElement);

			double alpha = getParameterForCategory(alphaCat, RFPParameter::alp, codon, false);
			double lambdaPrime = getParameterForCategory(lambdaPrimeCat, RFPParameter::lmPri, codon, false);

			double alphaPrime = alpha * gene.geneData.getCodonCount(codon);

			#ifndef STANDALONE
				RNGScope scope;
				NumericVector xx(1);
				xx = rgamma(1, alphaPrime, lambdaPrime);
				xx = rpois(1, xx[0]);
				gene.geneData.RFPObserved[codonIndex] = xx[0];
			#else

				std::gamma_distribution<double> GDistribution(alphaPrime, lambdaPrime);
				double tmp = GDistribution(Parameter::generator);
				std::poisson_distribution<unsigned> PDistribution(tmp);
				unsigned simulatedValue = PDistribution(Parameter::generator);
				tmpGene.geneData.setRFPObserved(codonIndex, simulatedValue);
			#endif
		}
		genome.addGene(tmpGene, true);
	}
}
