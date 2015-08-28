#include "include/FONSE/FONSEModel.h"
#include <vector>
#include <math.h>
#include <cfloat>
#include <iostream>
#include <array>
#include <limits.h>

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

#ifndef __APPLE__
#pragma omp parallel for private(mutation, selection, positions, curAA) reduction(+:likelihood,likelihood_proposed)
#endif
	for (int i = 0; i < getGroupListSize(); i++)
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
	int numCodons = parameter -> codonTable -> getNumCodons(grouping);
	double logLikelihood = 0.0;

	std::vector <unsigned> positions;
	double codonProb[6];

	std::array <unsigned, 8> codonRange = parameter -> codonTable -> AAToCodonRange(grouping); //checked

	unsigned maxIndexVal = 0u;
	for (int i = 1; i < (numCodons - 1); i++)
	{
		if (selection[maxIndexVal] < selection[i])
		{
			maxIndexVal = i;
		}
	}

	for (unsigned i = 0; i < 8; i++) {
		if (codonRange[i] == 100) break;
		positions = gene.geneData.getCodonPositions(codonRange[i]);
		for (unsigned j = 0; j < positions.size(); j++) {
			calculateCodonProbabilityVector(numCodons, positions[j], maxIndexVal, mutation, selection, phiValue, codonProb);
			for (int k = 0; k < numCodons; k++) {
				if (codonProb[k] == 0) continue;
				logLikelihood += std::log(codonProb[k]);
			}
		} 
		positions.clear();
	}

	return logLikelihood;
}

void FONSEModel::calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome, double& logAcceptanceRatioForAllMixtures)
{
	int numGenes = genome.getGenomeSize();
//	int numCodons = SequenceSummary::GetNumCodonsForAA(grouping);
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
	#pragma omp parallel for private(mutation, selection, mutation_proposed, selection_proposed, curAA, gene, seqsum) reduction(+:likelihood,likelihood_proposed)
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

void FONSEModel::calculateLogLikelihoodRatioForHyperParameters(Genome &genome, unsigned iteration, std::vector <double> & logProbabilityRatio)
{
	double currentSphi = getSphi(false);
	double currentMPhi = -(currentSphi * currentSphi) / 2;
	double lpr = 0.0; // this variable is only needed because OpenMP doesn't allow variables in reduction clause to be reference
	double proposedSphi = getSphi(true);
	double proposedMPhi = -(proposedSphi * proposedSphi) / 2;

	logProbabilityRatio.resize(1);

#ifndef __APPLE__
#pragma omp parallel for reduction(+:lpr)
#endif
	for (int i = 0u; i < genome.getGenomeSize(); i++)
	{
		unsigned mixture = getMixtureAssignment(i);
		mixture = getSynthesisRateCategory(mixture);
		double phi = getSynthesisRate(i, mixture, false);
		lpr += std::log(Parameter::densityLogNorm(phi, proposedMPhi, proposedSphi)) - std::log(Parameter::densityLogNorm(phi, currentMPhi, currentSphi));
	}

	lpr -= (std::log(currentSphi) - std::log(proposedSphi));
	logProbabilityRatio[0] = lpr;
}

void FONSEModel::adaptHyperParameterProposalWidths(unsigned adaptiveWidth)
{
	adaptSphiProposalWidth(adaptiveWidth);
}

void FONSEModel::updateHyperParameter(unsigned hp)
{
	switch (hp) {
	case 0:
		updateSphi();
		break;
	}
}

void FONSEModel::updateHyperParameterTraces(unsigned sample)
{
	updateSphiTrace(sample);
}

void FONSEModel::printHyperParameters()
{
	std::cout << "\t current Sphi estimate: " << getSphi() << std::endl;
	std::cout << "\t current Sphi proposal width: " << getCurrentSphiProposalWidth() << std::endl;
}

void FONSEModel::setParameter(FONSEParameter &_parameter)
{
	parameter = &_parameter;
}

void FONSEModel::calculateCodonProbabilityVector(unsigned numCodons, unsigned position, unsigned maxIndexValue, double *mutation, double *selection, double phi, double codonProb[])
{
	double denominator;

	/* c_i = exp[\Delta M - (\phi * \beta(i) * \Delta \omega)],                 *
	 * where \beta(i) = a_1 + (i * a_2)                                         *
	 *                                                                          *
	 * Right now a_1 and a_2 are set to 4.0. However, we are planning on making *
	 * them hyperparameters in the future, since they are constant for the      *
	 * entire genome.                                                           */

	if (selection[maxIndexValue] > 0.0) {
		denominator = 0.0;
		for (unsigned i = 0u; i < (numCodons - 1); i++) {
			codonProb[i] = std::exp(((mutation[i] - mutation[maxIndexValue])) + (phi * (4.0 + (4.0 * position)) * (selection[i] - selection[maxIndexValue])));
			denominator += codonProb[i];
		}
		codonProb[numCodons - 1] = std::exp((-1.0 * mutation[maxIndexValue]) - (phi * (4.0 + (4.0 * position)) * selection[maxIndexValue]));
		denominator += codonProb[numCodons - 1];
	}
	else {
		denominator = 1.0;
		for (unsigned i = 0u; i < (numCodons - 1); i++) {
			codonProb[i] = std::exp((mutation[i]) + (phi * (4.0 + (4.0 * position)) * selection[i]));
			denominator += codonProb[i];
		}
		codonProb[numCodons - 1] = 1.0;
	}

	for (unsigned i = 0; i < numCodons; i++) {
		codonProb[i] /= denominator;
	}
}
