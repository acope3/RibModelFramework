#include "include/FONSE/FONSEModel.h"
#include <vector>
#include <math.h>
#include <cfloat>
#include <iostream>
#include <array>
#include <limits.h>

double FONSEModel::calculateLogLikelihoodPerPositionPerGene(unsigned position, unsigned numCodons, double phiValue, double codonProb[])
{
	return 0.0;
}

FONSEModel::FONSEModel() : Model()
{
	parameter = nullptr;
}

FONSEModel::~FONSEModel()
{
	//dtor
}

void FONSEModel::calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k, double* logProbabilityRatio)
{
	double likelihood = 0.0;
	double likelihood_proposed = 0.0;
	std::string curAA;
	std::vector <unsigned> positions;
	double mutation[5];
	double selection[5];

	SequenceSummary seqsum = gene.getSequenceSummary();

	// get correct index for everything
	unsigned mutationCategory = parameter->getMutationCategory(k);
	unsigned selectionCategory = parameter->getSelectionCategory(k);
	unsigned expressionCategory = parameter->getSynthesisRateCategory(k);

	double phiValue = parameter->getSynthesisRate(geneIndex, expressionCategory, false);
	double phiValue_proposed = parameter->getSynthesisRate(geneIndex, expressionCategory, true);

	for (unsigned i = 0; i < getGroupListSize(); i++)
	{
		curAA = getGrouping(i);

		parameter->getParameterForCategory(mutationCategory, FONSEParameter::dM, curAA, false, mutation);
		parameter->getParameterForCategory(selectionCategory, FONSEParameter::dOmega, curAA, false, selection);

		likelihood += calculateLogLikelihoodRatioPerAA(gene, curAA, mutation, selection, phiValue);
		likelihood_proposed += calculateLogLikelihoodRatioPerAA(gene, curAA, mutation, selection, phiValue_proposed);
	}
	
	//std::cout << logLikelihood << " " << logLikelihood_proposed << std::endl;

	double sPhi = parameter->getSphi(false);
	double logPhiProbability = std::log(FONSEParameter::densityLogNorm(phiValue, (-(sPhi * sPhi) / 2), sPhi));
	double logPhiProbability_proposed = std::log(Parameter::densityLogNorm(phiValue_proposed, (-(sPhi * sPhi) / 2), sPhi));
	double currentLogLikelihood = (likelihood + logPhiProbability);
	double proposedLogLikelihood = (likelihood_proposed + logPhiProbability_proposed);

	logProbabilityRatio[0] = (proposedLogLikelihood - currentLogLikelihood) - (std::log(phiValue) - std::log(phiValue_proposed));
	logProbabilityRatio[1] = currentLogLikelihood - std::log(phiValue_proposed);
	logProbabilityRatio[2] = proposedLogLikelihood - std::log(phiValue);
}

double FONSEModel::calculateLogLikelihoodRatioPerAA(Gene& gene, std::string grouping, double *mutation, double *selection, double phiValue)
{
	int numCodons = SequenceSummary::GetNumCodonsForAA(grouping);
	double logLikelihood = 0.0;
	unsigned codonRange[2];
	std::vector <unsigned> positions;
	double *codonProb;

	SequenceSummary::AAToCodonRange(grouping, false, codonRange);

	unsigned maxIndexVal = 0u;
	for (int i = 1; i < (numCodons - 1); i++)
	{
		if (selection[maxIndexVal] < selection[i])
		{
			maxIndexVal = i;
		}
	}

	for (unsigned i = codonRange[0]; i < codonRange[1]; i++) {
		positions = gene.geneData.getCodonPositions(i);
		for (unsigned j = 0; j < positions.size(); j++) {
			codonProb = new double[numCodons]();
			calculateCodonProbabilityVector(numCodons, positions[j], maxIndexVal, mutation, selection, phiValue, codonProb);

			for (int k = 0; k < numCodons; k++) {
				if (isnan(codonProb[k])) {
					std::cout << "nan in per aa" << std::endl;
				}
				logLikelihood += std::log(codonProb[k]);
			}
			delete[] codonProb;
		} 
		positions.clear();
	}

	return logLikelihood;
}

void FONSEModel::calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome, double& logAcceptanceRatioForAllMixtures)
{
	int numGenes = genome.getGenomeSize();
	int numCodons = SequenceSummary::GetNumCodonsForAA(grouping);
	double likelihood = 0.0;
	double likelihood_proposed = 0.0;

	double mutation[5];
	double selection[5];
	double mutation_proposed[5];
	double selection_proposed[5];

	std::string curAA;

	Gene *gene;
	SequenceSummary *seqsum;

#ifndef __APPLE__
	#pragma omp parallel for private(mutation, selection, mutation_proposed, selection_proposed) reduction(+:likelihood,likelihood_proposed)
#endif
	for (int i = 0; i < numGenes; i++)
	{
		gene = &genome.getGene(i);
		seqsum = &gene->getSequenceSummary();
		if (seqsum->getAACountForAA(grouping) == 0) continue;

		// which mixture element does this gene belong to
		unsigned mixtureElement = parameter->getMixtureAssignment(i);
		// how is the mixture element defined. Which categories make it up
		unsigned mutationCategory = parameter->getMutationCategory(mixtureElement);
		unsigned selectionCategory = parameter->getSelectionCategory(mixtureElement);
		unsigned expressionCategory = parameter->getSynthesisRateCategory(mixtureElement);
		// get phi value, calculate likelihood conditional on phi
		double phiValue = parameter->getSynthesisRate(i, expressionCategory, false);


		// get current mutation and selection parameter
		parameter->getParameterForCategory(mutationCategory, FONSEParameter::dM, grouping, false, mutation);
		parameter->getParameterForCategory(selectionCategory, FONSEParameter::dOmega, grouping, false, selection);

		// get proposed mutation and selection parameter
		parameter->getParameterForCategory(mutationCategory, FONSEParameter::dM, grouping, true, mutation_proposed);
		parameter->getParameterForCategory(selectionCategory, FONSEParameter::dOmega, grouping, true, selection_proposed);

		likelihood += calculateLogLikelihoodRatioPerAA(*gene, grouping, mutation, selection, phiValue);
		likelihood_proposed += calculateLogLikelihoodRatioPerAA(*gene, grouping, mutation_proposed, selection_proposed, phiValue);

	}
	logAcceptanceRatioForAllMixtures = likelihood_proposed - likelihood;
}

void FONSEModel::setParameter(FONSEParameter &_parameter)
{
	parameter = &_parameter;
}

void FONSEModel::calculateCodonProbabilityVector(unsigned numCodons, unsigned position, unsigned maxIndexValue, double *mutation, double *selection, double phi, double codonProb[])
{
	double denominator;

	if (selection[maxIndexValue] > 0.0) {
		denominator = 0.0;
		for (unsigned i = 0u; i < (numCodons - 1); i++) {
			codonProb[i] = std::exp(((mutation[i] - mutation[maxIndexValue])) + (phi * position * (selection[i] - selection[maxIndexValue])));
			denominator += codonProb[i];
		}
		codonProb[numCodons - 1] = std::exp((-1.0 * mutation[maxIndexValue]) - (phi * position * selection[maxIndexValue]));
		denominator += codonProb[numCodons - 1];
	}
	else {
		denominator = 1.0;
		for (unsigned i = 0u; i < (numCodons - 1); i++) {
			codonProb[i] = std::exp((mutation[i]) + (phi * position * selection[i]));
			denominator += codonProb[i];
		}
		codonProb[numCodons - 1] = 1.0;
	}

	for (unsigned i = 0; i < numCodons; i++) {
		codonProb[i] /= denominator;
		if (isnan(codonProb[i])) {
			std::cout << "nan" << std::endl;
		}
	}
}

/*std::vector<double> FONSEModel::CalculateProbabilitiesForCodons(std::vector<double> mutation, std::vector<double> selection, double phi)
{
unsigned numCodons = mutation.size() + 1;
double* _mutation = &mutation[0];
double* _selection = &selection[0];
double* codonProb = new double[numCodons]();
calculateCodonProbability(numCodons, _mutation, _selection, phi, codonProb);
std::vector<double> returnVector(codonProb, codonProb + numCodons);
return returnVector;
}*/
