#include "include/FONSE/FONSEModel.h"
#include <vector>
#include <math.h>
#include <cfloat>
#include <iostream>
#include <array>

FONSEModel::FONSEModel() : Model()
{
	parameter = nullptr;
}

FONSEModel::~FONSEModel()
{
	//dtor
}

void FONSEModel::calculateLogLikelihoodRatioPerGene(Gene& gene, int geneIndex, unsigned k, double* logProbabilityRatio)
{
	double logLikelihood = 0.0;
	double logLikelihood_proposed = 0.0;

	SequenceSummary seqsum = gene.getSequenceSummary();

	// get correct index for everything
	unsigned mutationCategory = parameter->getMutationCategory(k);
	unsigned selectionCategory = parameter->getSelectionCategory(k);
	unsigned expressionCategory = parameter->getSynthesisRateCategory(k);

	double phiValue = parameter->getSynthesisRate(geneIndex, expressionCategory, false);
	double phiValue_proposed = parameter->getSynthesisRate(geneIndex, expressionCategory, true);

	std::vector <double> *mutation;
	std::vector <double> *selection;
	mutation = parameter->getParameterForCategory(mutationCategory, FONSEParameter::dM, false);
	selection = parameter->getParameterForCategory(selectionCategory, FONSEParameter::dOmega, false);
	for (unsigned i = 0; i < 65; i++)
	{
		std::vector <unsigned> positions = seqsum.getCodonPositions(i);
		if (positions.size() == 0) continue;

		for (unsigned j = 0; j < positions.size(); j++) {
			logLikelihood += calculateLogLikelihoodPerPositionPerGene(positions[j], i, mutation, selection, phiValue);
			logLikelihood_proposed += calculateLogLikelihoodPerPositionPerGene(positions[j], i, mutation, selection, phiValue_proposed);
		}
	}

	double sPhi = parameter->getSphi(false);
	double logPhiProbability = std::log(FONSEParameter::densityLogNorm(phiValue, (-(sPhi * sPhi) / 2), sPhi));
	double logPhiProbability_proposed = std::log(Parameter::densityLogNorm(phiValue_proposed, (-(sPhi * sPhi) / 2), sPhi));
	double currentLogLikelihood = (logLikelihood + logPhiProbability);
	double proposedLogLikelihood = (logLikelihood_proposed + logPhiProbability_proposed);

	logProbabilityRatio[0] = (proposedLogLikelihood - currentLogLikelihood) - (std::log(phiValue) - std::log(phiValue_proposed));
	logProbabilityRatio[1] = currentLogLikelihood - std::log(phiValue_proposed);
	logProbabilityRatio[2] = proposedLogLikelihood - std::log(phiValue);
}

void FONSEModel::calculateCodonProbabilityVector(unsigned numCodons, double mutation[], double selection[], double phi, double codonProb[])
{
	// calculate numerator and denominator for codon probabilities
	unsigned minIndexVal = 0u;
	double denominator;
	for (unsigned i = 1u; i < (numCodons - 1); i++)
	{
		if (selection[minIndexVal] > selection[i])
		{
			minIndexVal = i;
		}
	}

	// if the min(selection) is less than zero than we have to adjust the reference codon.
	// if the reference codon is the min value (0) than we do not have to adjust the reference codon.
	// This is necessary to deal with very large phi values (> 10^4) and avoid  producing Inf which then
	// causes the denominator to be Inf (Inf / Inf = NaN).
	if (selection[minIndexVal] < 0.0)
	{
		denominator = 0.0;
		for (unsigned i = 0; i < (numCodons - 1); i++)
		{
			codonProb[i] = std::exp(-(mutation[i] - mutation[minIndexVal]) - ((selection[i] - selection[minIndexVal]) * phi));
			//codonProb[i] = std::exp( -mutation[i] - (selection[i] * phi) );
			denominator += codonProb[i];
		}
		// alphabetically last codon is reference codon!
		codonProb[numCodons - 1] = std::exp(mutation[minIndexVal] + selection[minIndexVal] * phi);
		denominator += codonProb[numCodons - 1];
	}
	else{
		denominator = 1.0;
		for (unsigned i = 0; i < (numCodons - 1); i++)
		{
			codonProb[i] = std::exp(-mutation[i] - (selection[i] * phi));
			denominator += codonProb[i];
		}
		// alphabetically last codon is reference codon!
		codonProb[numCodons - 1] = 1.0;
	}
	// normalize codon probabilities
	for (unsigned i = 0; i < numCodons; i++)
	{
		codonProb[i] = codonProb[i] / denominator;
	}
}

double FONSEModel::calculateLogLikelihoodPerPositionPerGene(unsigned position, unsigned codonIndex, std::vector <double> *mutation, std::vector <double> *selection, double phiValue)
{
	double numerator = 0.0;
	double denominator = 0.0;
	std::vector <unsigned> synonymous;

	synonymous = SequenceSummary::getSynonymousCodonIndicesByCodonIndex(codonIndex);

	double mut = mutation->at(codonIndex);
	double sel = selection->at(codonIndex);
	double a1 = 4.0;
	double a2 = 4.0;
	double q = 1;
	double Ne = 1;

	numerator = std::exp(std::log(mut) + (sel * (a1 - a2) * (-1.0 * q * Ne * phiValue)) + (sel * a2 * (-1.0 * q * Ne * phiValue) * position));

	for (unsigned i = 0u; i < synonymous.size(); i++) {
		denominator += std::exp(std::log(mutation->at(synonymous[i])) + (selection->at(synonymous[i]) * (a1 - a2) * (-1.0 * q * Ne * phiValue)) + (selection->at(synonymous[i]) * a2 * (-1.0 * q * Ne * phiValue)));
	}
	return std::log(numerator / denominator);
}

void FONSEModel::calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome, double& logAcceptanceRatioForAllMixtures)
{
	int numGenes = genome.getGenomeSize();
	int numCodons = SequenceSummary::GetNumCodonsForAA(grouping);
	double likelihood = 0.0;
	double likelihood_proposed = 0.0;

	std::vector <double> *mutation;
	std::vector <double> *selection;
	std::vector <double> *mutation_proposed;
	std::vector <double> *selection_proposed;
	std::vector <unsigned> synonymous;
	std::vector <unsigned> positions;
#ifndef __APPLE__
	//#pragma omp parallel for private(mutation, selection, mutation_proposed, selection_proposed, codonCount) reduction(+:likelihood,likelihood_proposed)
#endif
	for (int i = 0; i < numGenes; i++)
	{
		Gene gene = genome.getGene(i);
		SequenceSummary seqsum = gene.getSequenceSummary();
		if (seqsum.getAAcountForAA(grouping) == 0) continue;

		// which mixture element does this gene belong to
		unsigned mixtureElement = parameter->getMixtureAssignment(i);
		// how is the mixture element defined. Which categories make it up
		unsigned mutationCategory = parameter->getMutationCategory(mixtureElement);
		unsigned selectionCategory = parameter->getSelectionCategory(mixtureElement);
		unsigned expressionCategory = parameter->getSynthesisRateCategory(mixtureElement);
		// get phi value, calculate likelihood conditional on phi
		double phiValue = parameter->getSynthesisRate(i, expressionCategory, false);

		// get current mutation and selection parameter
		mutation = parameter->getParameterForCategory(mutationCategory, FONSEParameter::dM, false);
		selection = parameter->getParameterForCategory(selectionCategory, FONSEParameter::dOmega, false);

		// get proposed mutation and selection parameter
		mutation_proposed = parameter->getParameterForCategory(mutationCategory, FONSEParameter::dM, true);
		selection_proposed = parameter->getParameterForCategory(selectionCategory, FONSEParameter::dOmega, true);

		for (unsigned i = 0; i < 65; i++)
		{
			std::vector <unsigned> positions = seqsum.getCodonPositions(i);
			if (positions.size() == 0) continue;

			for (unsigned j = 0; j < positions.size(); j++) {
				likelihood += calculateLogLikelihoodPerPositionPerGene(positions[j], i, mutation, selection, phiValue);
				likelihood_proposed += calculateLogLikelihoodPerPositionPerGene(positions[j], i, mutation_proposed, selection_proposed, phiValue);
			}

		}
	}
	logAcceptanceRatioForAllMixtures = likelihood_proposed - likelihood;
}

void FONSEModel::obtainCodonCount(SequenceSummary& seqsum, std::string curAA, int codonCount[])
{
	unsigned codonRange[2];
	SequenceSummary::AAToCodonRange(curAA, false, codonRange);
	// get codon counts for AA
	unsigned j = 0u;
	for (unsigned i = codonRange[0]; i < codonRange[1]; i++, j++)
	{
		codonCount[j] = seqsum.getCodonCountForCodon(i);
	}
}

void FONSEModel::setParameter(FONSEParameter &_parameter)
{
	parameter = &_parameter;
}

std::vector<double> FONSEModel::CalculateProbabilitiesForCodons(std::vector<double> mutation, std::vector<double> selection, double phi)
{
	unsigned numCodons = mutation.size() + 1;
	double* _mutation = &mutation[0];
	double* _selection = &selection[0];
	double* codonProb = new double[numCodons]();
	calculateCodonProbabilityVector(numCodons, _mutation, _selection, phi, codonProb);
	std::vector<double> returnVector(codonProb, codonProb + numCodons);
	return returnVector;
}
