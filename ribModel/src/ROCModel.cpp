#include "include/ROC/ROCModel.h"
#include <vector>
#include <math.h>
#include <cfloat>
#include <iostream>
#include <array>

ROCModel::ROCModel() : Model()
{
	parameter = nullptr;
}

ROCModel::~ROCModel()
{
	//dtor
}

void ROCModel::calculateLogLikelihoodRatioPerGene(Gene& gene, int geneIndex, unsigned k, double* logProbabilityRatio)
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

	double mutation[5];
	double selection[5];
	int codonCount[6];
#ifndef __APPLE__
#pragma omp parallel for private(mutation, selection, codonCount) reduction(+:logLikelihood,logLikelihood_proposed)
#endif
	for(unsigned i = 0; i < getGroupListSize(); i++)
	{
		std::string curAA = getGrouping(i);

		// skip amino acids which do not occur in current gene. Avoid useless calculations and multiplying by 0
		if(seqsum.getAACountForAA(curAA) == 0) continue;

		// get codon count (total count not parameter->count)
		int numCodons = seqsum.GetNumCodonsForAA(curAA);
		// get mutation and selection parameter->for gene
		//double* mutation = new double[numCodons - 1]();
		parameter->getParameterForCategory(mutationCategory, ROCParameter::dM, curAA, false, mutation);
		//double* selection = new double[numCodons - 1]();
		parameter->getParameterForCategory(selectionCategory, ROCParameter::dEta, curAA, false, selection);

		// prepare array for codon counts for AA
		//int* codonCount = new int[numCodons]();
		obtainCodonCount(seqsum, curAA, codonCount);

		logLikelihood += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation, selection, phiValue);
		logLikelihood_proposed += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation, selection, phiValue_proposed);

		//delete [] &mutation;
		//delete [] &selection;
		//delete [] &codonCount;
	}

	double sPhi = parameter->getSphi(false);
	double logPhiProbability = std::log(ROCParameter::densityLogNorm(phiValue, (-(sPhi * sPhi) / 2), sPhi));
	double logPhiProbability_proposed = std::log(Parameter::densityLogNorm(phiValue_proposed, (-(sPhi * sPhi) / 2), sPhi));
	double currentLogLikelihood = (logLikelihood + logPhiProbability);
	double proposedLogLikelihood = (logLikelihood_proposed + logPhiProbability_proposed);

	logProbabilityRatio[0] = (proposedLogLikelihood - currentLogLikelihood) - (std::log(phiValue) - std::log(phiValue_proposed));
	logProbabilityRatio[1] = currentLogLikelihood - std::log(phiValue_proposed);
	logProbabilityRatio[2] = proposedLogLikelihood - std::log(phiValue);
}

void ROCModel::calculateCodonProbabilityVector(unsigned numCodons, double mutation[], double selection[], double phi, double codonProb[])
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
	if(selection[minIndexVal] < 0.0)
	{
		denominator = 0.0;
		for(unsigned i = 0; i < (numCodons - 1); i++)
		{
			codonProb[i] = std::exp( -(mutation[i] - mutation[minIndexVal]) - ((selection[i] - selection[minIndexVal]) * phi) );
			//codonProb[i] = std::exp( -mutation[i] - (selection[i] * phi) );
			denominator += codonProb[i];
		}
		// alphabetically last codon is reference codon!
		codonProb[numCodons - 1] = std::exp(mutation[minIndexVal] + selection[minIndexVal] * phi);
		denominator += codonProb[numCodons - 1];
	}else{
		denominator = 1.0;
		for(unsigned i = 0; i < (numCodons - 1); i++)
		{
			codonProb[i] = std::exp( -mutation[i] - (selection[i] * phi) );
			denominator += codonProb[i];
		}
		// alphabetically last codon is reference codon!
		codonProb[numCodons - 1] = 1.0;
	}
	// normalize codon probabilities
	for(unsigned i = 0; i < numCodons; i++)
	{
		codonProb[i] = codonProb[i] / denominator;
	}
}

double ROCModel::calculateLogLikelihoodPerAAPerGene(unsigned numCodons, int codonCount[], double mutation[], double selection[], double phiValue)
{
	double logLikelihood = 0.0;
	// calculate codon probabilities
	double* codonProbabilities = new double[numCodons]();
	calculateCodonProbabilityVector(numCodons, mutation, selection, phiValue, codonProbabilities);

	// calculate likelihood for current AA for this combination of selection and mutation category
	for(unsigned i = 0; i < numCodons; i++)
	{
		if (codonCount[i] == 0) continue;
		logLikelihood += std::log(codonProbabilities[i]) * codonCount[i];
	}
	//std::cout <<"deleting codonProbabilities\n";
	delete [] codonProbabilities;
	//std::cout <<"DONEdeleting codonProbabilities\n";
	return logLikelihood;
}

void ROCModel::calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome, double& logAcceptanceRatioForAllMixtures)
{
	int numGenes = genome.getGenomeSize();
	int numCodons = SequenceSummary::GetNumCodonsForAA(grouping);
	double likelihood = 0.0;
	double likelihood_proposed = 0.0;

	double mutation[5];
	double selection[5];
	double mutation_proposed[5];
	double selection_proposed[5];
	int codonCount[6];
#ifndef __APPLE__
#pragma omp parallel for private(mutation, selection, mutation_proposed, selection_proposed, codonCount) reduction(+:likelihood,likelihood_proposed)
#endif
	for(int i = 0; i < numGenes; i++)
	{
		Gene gene = genome.getGene(i);
		SequenceSummary seqsum = gene.getSequenceSummary();
		if(seqsum.getAACountForAA(grouping) == 0) continue;

		// which mixture element does this gene belong to
		unsigned mixtureElement = parameter->getMixtureAssignment(i);
		// how is the mixture element defined. Which categories make it up
		unsigned mutationCategory = parameter->getMutationCategory(mixtureElement);
		unsigned selectionCategory = parameter->getSelectionCategory(mixtureElement);
		unsigned expressionCategory = parameter->getSynthesisRateCategory(mixtureElement);
		// get phi value, calculate likelihood conditional on phi
		double phiValue = parameter->getSynthesisRate(i, expressionCategory, false);

		// get current mutation and selection parameter
		//double* mutation = new double[numCodons - 1]();
		parameter->getParameterForCategory(mutationCategory, ROCParameter::dM, grouping, false, mutation);
		//double* selection = new double[numCodons - 1]();
		parameter->getParameterForCategory(selectionCategory, ROCParameter::dEta, grouping, false, selection);

		// get proposed mutation and selection parameter
		//double* mutation_proposed = new double[numCodons - 1]();
		parameter->getParameterForCategory(mutationCategory, ROCParameter::dM, grouping, true, mutation_proposed);
		//double* selection_proposed = new double[numCodons - 1]();
		parameter->getParameterForCategory(selectionCategory, ROCParameter::dEta, grouping, true, selection_proposed);

		//int* codonCount = new int[numCodons]();
		obtainCodonCount(seqsum, grouping, codonCount);
		likelihood += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation, selection, phiValue);
		likelihood_proposed += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation_proposed, selection_proposed, phiValue);
		//delete [] &codonCount;
		//delete [] &mutation;
		//delete [] &selection;
		//delete [] &mutation_proposed;
		//delete [] &selection_proposed;
	}
	logAcceptanceRatioForAllMixtures = likelihood_proposed - likelihood;
}

void ROCModel::obtainCodonCount(SequenceSummary& seqsum, std::string curAA, int codonCount[])
{
	unsigned codonRange[2];
	SequenceSummary::AAToCodonRange(curAA, false, codonRange);
	// get codon counts for AA
	unsigned j = 0u;
	for(unsigned i = codonRange[0]; i < codonRange[1]; i++, j++)
	{
		codonCount[j] = seqsum.getCodonCountForCodon(i);
	}
}


void ROCModel::setParameter(ROCParameter &_parameter)
{
	parameter = &_parameter;
}

std::vector<double> ROCModel::CalculateProbabilitiesForCodons(std::vector<double> mutation, std::vector<double> selection, double phi)
{
	unsigned numCodons = mutation.size() + 1;
	double* _mutation = &mutation[0];
	double* _selection = &selection[0];
	double* codonProb = new double[numCodons]();
	calculateCodonProbabilityVector(numCodons, _mutation, _selection, phi, codonProb);
	std::vector<double> returnVector(codonProb, codonProb + numCodons);
	return returnVector;
}



void ROCModel::simulateGenome(Genome &genome)
{


     int k;
     int aaCount;

     unsigned codonIndex;
     unsigned aaRange[2];
     std::string curAA;


	//std::srand(std::time(0));
	//simulatedGenes.resize(genes.size());
	std::string tmpDesc = "Simulated Gene";


	for (unsigned geneIndex = 0; geneIndex < genome.getGenomeSize(); geneIndex++) //loop over all genes in the genome
	{
		Gene gene = genome.getGene(geneIndex);
		SequenceSummary seqSum = gene.geneData;
		std::string tmpSeq = "ATG"; //Always will have the start amino acid


		unsigned mixtureElement = getMixtureAssignment(geneIndex);
		unsigned mutationCategory = getMutationCategory(mixtureElement);
		unsigned selectionCategory = getSelectionCategory(mixtureElement);
		unsigned synthesisRateCategory = getSynthesisRateCategory(mixtureElement);
		double phi = getSynthesisRate(geneIndex, synthesisRateCategory, false);

		std::string geneSeq = gene.getSequence();
		for (unsigned position = 1; position < (geneSeq.size() / 3); position++)
	 	{
	 		std::string codon = geneSeq.substr((position * 3), 3);
			std::string aa = SequenceSummary::CodonToAA(codon);
			unsigned numCodons = SequenceSummary::GetNumCodonsForAA(aa);

			double* codonProb = new double[numCodons](); //size the arrays to the proper size based on # of codons.
			double* mutation = new double[numCodons - 1]();
			double* selection = new double[numCodons - 1]();

			getParameterForCategory(mutationCategory, ROCParameter::dM, curAA, false, mutation);
			getParameterForCategory(selectionCategory, ROCParameter::dEta, curAA, false, selection);
			calculateCodonProbabilityVector(numCodons, mutation, selection, phi, codonProb);

			codonIndex = Parameter::randMultinom(codonProb, numCodons);
			SequenceSummary::AAToCodonRange(curAA, false, aaRange); //need the first spot in the array where the codons for curAA are
			codon = seqSum.IndexToCodon(aaRange[0] + codonIndex);//get the correct codon based off codonIndex
			tmpSeq += codon;
	 	}


/*

		std::ostringstream strstream;
		strstream << mixtureElement;
		std::string tmpID = gene.getId() + "_MixtureElement" + strstream.str();
		for (j = 0; j < 22; j++) //loop over each amino acid, naa[]
		{
			aaCount = seqSum.getAAcountForAA(j);
			curAA = seqSum.AminoAcidArray[j];
			if (curAA == "X") continue;
			numCodons = seqSum.GetNumCodonsForAA(curAA);
			if (curAA == "M") aaCount -= 1;
			double* codonProb = new double[numCodons](); //size the arrays to the proper size based on # of codons.
			double* mutation = new double[numCodons - 1]();
			double* selection = new double[numCodons - 1]();

			//get the probability vector for each amino acid
			if (curAA == "M" || curAA == "W")
			{
				codonProb[0] = 1;
			}
			else
			{
				model.getParameterForCategory(mutationCategory, ROCParameter::dM, curAA, false, mutation);
				model.getParameterForCategory(selectionCategory, ROCParameter::dEta, curAA, false, selection);
				model.calculateCodonProbabilityVector(numCodons, mutation, selection, phi, codonProb);
			}
			for (k = 0; k < aaCount; k++)
			{
				codonIndex = ROCParameter::randMultinom(codonProb, numCodons);
				seqSum.AAToCodonRange(curAA, false, aaRange); //need the first spot in the array where the codons for curAA are
				codon = seqSum.IndexToCodon(aaRange[0] + codonIndex);//get the correct codon based off codonIndex
				tmpSeq += codon;
			}
		}

		codon =	seqSum.IndexToCodon((rand() % 3) + 61); //randomly choose a stop codon, from range 61-63
		tmpSeq += codon;
		Gene tmpGene(tmpSeq, tmpID, tmpDesc);
		simulatedGenes[i] = tmpGene;
	}
*/}
}
