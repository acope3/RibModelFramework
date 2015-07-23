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
	std::string curAA;
	unsigned codonRange[2];
	std::vector <unsigned> positions;
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
	for (unsigned i = 0; i < getGroupListSize(); i++)
	{
		curAA = getGrouping(i);

		SequenceSummary::AAToCodonRange(curAA, false, codonRange);

		for (unsigned j = codonRange[0]; j < codonRange[1]; j++) {
			positions = seqsum.getCodonPositions(j);
			if (positions.size() == 0) continue;
			for (unsigned k = 0; k < positions.size(); k++) {
				logLikelihood += calculateLogLikelihoodPerPositionPerGene(positions[k], curAA, j, mutation, selection, phiValue);
				logLikelihood_proposed += calculateLogLikelihoodPerPositionPerGene(positions[k], curAA, j, mutation, selection, phiValue_proposed);
			}
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

double FONSEModel::calculateCodonProbability(unsigned position, double mutation[], double selection[], double phi)
{
	return 0.0;
}

double FONSEModel::calculateLogLikelihoodPerPositionPerGene(unsigned position, std::string curAA, unsigned codonIndex, std::vector <double> *mutation, std::vector <double> *selection, double phiValue)
{
	std::vector <unsigned> synonymous;
	unsigned paramRange[2];
	unsigned codonRange[2];
	unsigned paramIndex;
	double mut = 0.0;
	double sel = 0.0;

	double numerator = 0.0;
	double denominator = 1.0;

	// if we are the last codon alphabetically, then we are the reference codon
	SequenceSummary::AAToCodonRange(curAA, true, paramRange);
	SequenceSummary::AAToCodonRange(curAA, false, codonRange);

	paramIndex = paramRange[0] + (codonIndex - codonRange[0]);

	// if we are the reference codon, the numerator is simply 1
	if (codonIndex = paramRange[1] - 1) {
		numerator = 1.0;
	}
	else {
		numerator = (-1.0 * mutation->at(paramIndex)) - (phiValue * position * selection->at(paramIndex));
	}

	for (unsigned i = paramRange[0]; i < paramRange[1]; i++) {
		denominator += (-1.0 * mutation->at(i)) - (phiValue * position * selection->at(i));
	}

	/* Below is the version of the formulation from the REU paper. However, the old FONSE used a much simpler formulation
	numerator = std::exp(std::log(mut) + (sel * (a1 - a2) * (-1.0 * q * Ne * phiValue)) + (sel * a2 * (-1.0 * q * Ne * phiValue) * position));

	for (unsigned i = aaRange[0]; i < aaRange[1]; i++) {
	denominator += std::exp(std::log(mutation->at(i)) + (selection->at(i) * (a1 - a2) * (-1.0 * q * Ne * phiValue)) + (selection->at(i) * a2 * (-1.0 * q * Ne * phiValue)));
	}

	*/


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
	unsigned paramRange[2];
	unsigned aaRange[2];
	std::string curAA;

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

		SequenceSummary::AAToCodonRange(grouping, false, aaRange);

		for (unsigned j = aaRange[0]; j < aaRange[1]; j++) {
			positions = seqsum.getCodonPositions(j);
			if (positions.size() == 0) continue;
			for (unsigned k = 0; k < positions.size(); k++) {
				likelihood += calculateLogLikelihoodPerPositionPerGene(positions[k], grouping, j, mutation, selection, phiValue);
				likelihood_proposed += calculateLogLikelihoodPerPositionPerGene(positions[k], grouping, j, mutation_proposed, selection_proposed, phiValue);
			}
			positions.clear();
		}
	}
	logAcceptanceRatioForAllMixtures = likelihood_proposed - likelihood;
}

void FONSEModel::setParameter(FONSEParameter &_parameter)
{
	parameter = &_parameter;
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
