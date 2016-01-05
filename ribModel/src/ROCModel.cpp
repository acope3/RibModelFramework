#include "include/ROC/ROCModel.h"

ROCModel::ROCModel(bool _withPhi) : Model()
{
	parameter = 0;
	withPhi = _withPhi;
}

ROCModel::~ROCModel()
{
	//dtor
}

void ROCModel::calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k, double* logProbabilityRatio)
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
	for(int i = 0; i < getGroupListSize(); i++)
	{
		std::string curAA = getGrouping(i);

		// skip amino acids which do not occur in current gene. Avoid useless calculations and multiplying by 0
		if(seqsum.getAACountForAA(curAA) == 0) continue;

		// get number of codons for AA (total number not parameter->count)
		unsigned numCodons = seqsum.GetNumCodonsForAA(curAA);
		// get mutation and selection parameter->for gene
		parameter->getParameterForCategory(mutationCategory, ROCParameter::dM, curAA, false, mutation);
		parameter->getParameterForCategory(selectionCategory, ROCParameter::dEta, curAA, false, selection);
		// get codon occurence in sequence
		obtainCodonCount(seqsum, curAA, codonCount);

		logLikelihood += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation, selection, phiValue);
		logLikelihood_proposed += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation, selection, phiValue_proposed);
	}
	unsigned mixture = getMixtureAssignment(geneIndex);
	mixture = getSynthesisRateCategory(mixture);
	double sPhi = parameter->getSphi(mixture, false);
	double mPhi = (-(sPhi * sPhi) / 2);
	double logPhiProbability = Parameter::densityLogNorm(phiValue, mPhi, sPhi, true);
	double logPhiProbability_proposed = Parameter::densityLogNorm(phiValue_proposed, mPhi, sPhi, true);

	// TODO: make this work for more than one phi value, or for genes that don't have phi values
	if (withPhi) {
		for (unsigned i = 0; i < parameter->getNumObservedPhiSets(); i++) {
			double obsPhi = gene.getObservedSynthesisRate(i);
			if (obsPhi > -1.0) {
				//logPhiProbability += Parameter::densityLogNorm(obsPhi + getAphi(i), std::log(phiValue), getSepsilon(i), true);
				//logPhiProbability_proposed += Parameter::densityLogNorm(obsPhi + getAphi(i), std::log(phiValue_proposed), getSepsilon(i), true);

				logPhiProbability += Parameter::densityLogNorm(obsPhi, std::log(phiValue) + getAphi(i), getSepsilon(i), true);
				logPhiProbability_proposed += Parameter::densityLogNorm(obsPhi, std::log(phiValue_proposed) + getAphi(i), getSepsilon(i), true);

			}
		}
	}

	double currentLogLikelihood = (logLikelihood + logPhiProbability);
	double proposedLogLikelihood = (logLikelihood_proposed + logPhiProbability_proposed);

	logProbabilityRatio[0] = (proposedLogLikelihood - currentLogLikelihood) - (std::log(phiValue) - std::log(phiValue_proposed));
	logProbabilityRatio[1] = currentLogLikelihood - std::log(phiValue_proposed);
	logProbabilityRatio[2] = proposedLogLikelihood - std::log(phiValue);
	logProbabilityRatio[3] = currentLogLikelihood;
	logProbabilityRatio[4] = proposedLogLikelihood;

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
	// This is necessary to deal with very large phi values (> 10^4) and avoid producing Inf which then
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
	double codonProbabilities[6];
	calculateCodonProbabilityVector(numCodons, mutation, selection, phiValue, codonProbabilities);

	// calculate likelihood for current AA for this combination of selection and mutation category
	for(unsigned i = 0; i < numCodons; i++)
	{
		if (codonCount[i] == 0) continue;
		logLikelihood += std::log(codonProbabilities[i]) * codonCount[i];
	}
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
	Gene *gene;
	SequenceSummary *seqsum;
#ifndef __APPLE__
//#pragma omp parallel for private(mutation, selection, mutation_proposed, selection_proposed, codonCount, gene, seqsum) reduction(+:likelihood,likelihood_proposed)
#endif
	for(int i = 0; i < numGenes; i++)
	{
		gene = &genome.getGene(i);
		seqsum = &gene->getSequenceSummary();
		if(seqsum->getAACountForAA(grouping) == 0) continue;

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
		obtainCodonCount(*seqsum, grouping, codonCount);
		likelihood += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation, selection, phiValue);
		likelihood_proposed += calculateLogLikelihoodPerAAPerGene(numCodons, codonCount, mutation_proposed, selection_proposed, phiValue);
	}

	likelihood_proposed = likelihood_proposed + calculateMutationPrior(grouping, true);
	likelihood = likelihood + calculateMutationPrior(grouping, false);

	logAcceptanceRatioForAllMixtures = (likelihood_proposed - likelihood);
}

double ROCModel::calculateAllPriors()
{
	double prior = 0.0;
	unsigned size = getGroupListSize();

	for(unsigned i = 0; i < size; i++)
	{
		std::string grouping = getGrouping(i);
		prior += calculateMutationPrior(grouping, false);
	}

	// add more priors if necessary.

	return prior;
}

double ROCModel::calculateMutationPrior(std::string grouping, bool proposed)
{
	unsigned numCodons = SequenceSummary::GetNumCodonsForAA(grouping);
	double mutation[5];

	double priorValue = 0.0;

	unsigned numMutCat = parameter->getNumMutationCategories();
	double mutation_prior_sd = parameter->getMutationPriorStandardDeviation();
	for(unsigned i = 0u; i < numMutCat; i++)
	{
		parameter->getParameterForCategory(i, ROCParameter::dM, grouping, proposed, mutation);
		for(unsigned k = 0u; k < numCodons; k++)
		{
			priorValue += Parameter::densityNorm(mutation[k], 0.0, mutation_prior_sd, true);
		}
	}
	return priorValue;
}

void ROCModel::obtainCodonCount(SequenceSummary& seqsum, std::string curAA, int codonCount[])
{
	std::array <unsigned, 2> codonRange = SequenceSummary::AAToCodonRange(curAA);
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

void ROCModel::printHyperParameters()
{
	for(unsigned i = 0u; i < getNumSynthesisRateCategories(); i++)
	{
#ifndef STANDALONE
		Rprintf("Current Sphi estimate for selection category %d: %f\n", i, getSphi(i, false));
#else
		std::cout << "Current Sphi estimate for selection category " << i << ": " << getSphi(i, false) << std::endl;
#endif
	}
#ifndef STANDALONE
	Rprintf("\t current Sphi proposal width: %f\n", getCurrentSphiProposalWidth());
#else
	std::cout << "\t current Sphi proposal width: " << getCurrentSphiProposalWidth() << std::endl;
#endif
	if(withPhi)
	{
#ifndef STANDALONE
	Rprintf("\t current Aphi estimates:");
#else
	std::cout << "\t current Aphi estimates:";
#endif
		for (unsigned i = 0; i < getNumPhiGroupings(); i++)
		{
#ifndef STANDALONE
			Rprintf(" %f", getAphi(i, false));
#else
			std::cout << " " << getAphi(i, false);
#endif
		}
#ifndef STANDALONE
		Rprintf("\n\t current Aphi proposal widths:");
#else
		std::cout << "\n\t current Aphi proposal widths:";
#endif
		for (unsigned i = 0; i < getNumPhiGroupings(); i++)
		{
#ifndef STANDALONE
			Rprintf(" %f", getCurrentAphiProposalWidth(i));
#else
			std::cout << " " << getCurrentAphiProposalWidth(i);
#endif
		}
#ifndef STANDALONE
		Rprintf("\n\t current Sepsilon estimates:");
#else
		std::cout << "\n\t current Sepsilon estimates:";
#endif
		for (unsigned i = 0; i < getNumPhiGroupings(); i++)
		{
#ifndef STANDALONE
			Rprintf(" %f", getSepsilon(i));
#else
			std::cout << " " << getSepsilon(i);
#endif
		}
#ifndef STANDALONE
		Rprintf("\n");
#else
		std::cout << std::endl;
#endif
	}
}


void ROCModel::simulateGenome(Genome &genome)
{
	unsigned codonIndex;

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
			std::string aa = SequenceSummary::codonToAA(codon);

			if (aa == "X") continue;

			unsigned numCodons = SequenceSummary::GetNumCodonsForAA(aa);

			double* codonProb = new double[numCodons](); //size the arrays to the proper size based on # of codons.
			double* mutation = new double[numCodons - 1]();
			double* selection = new double[numCodons - 1]();


			if (aa == "M" || aa == "W")
			{
				codonProb[0] = 1;
			}
			else
			{
				getParameterForCategory(mutationCategory, ROCParameter::dM, aa, false, mutation);
				getParameterForCategory(selectionCategory, ROCParameter::dEta, aa, false, selection);
				calculateCodonProbabilityVector(numCodons, mutation, selection, phi, codonProb);
			}


			codonIndex = Parameter::randMultinom(codonProb, numCodons);
			std::array <unsigned, 2> aaRange = SequenceSummary::AAToCodonRange(aa); //need the first spot in the array where the codons for curAA are
			codon = seqSum.indexToCodon(aaRange[0] + codonIndex);//get the correct codon based off codonIndex
			tmpSeq += codon;
	 	}
		std::string codon =	seqSum.indexToCodon((unsigned)Parameter::randUnif(61.0, 63.0)); //randomly choose a stop codon, from range 61-63
		tmpSeq += codon;
		Gene simulatedGene(tmpSeq, tmpDesc, gene.getId());
		genome.addGene(simulatedGene, true);
	}
}

void ROCModel::calculateLogLikelihoodRatioForHyperParameters(Genome &genome, unsigned iteration, std::vector <double> &logProbabilityRatio)
{	
	double lpr = 0.0;
	unsigned selectionCategory = getNumSynthesisRateCategories();
	std::vector<double> currentSphi(selectionCategory, 0.0);
	std::vector<double> currentMphi(selectionCategory, 0.0);
	std::vector<double> proposedSphi(selectionCategory, 0.0);
	std::vector<double> proposedMphi(selectionCategory, 0.0);


	for(unsigned i = 0u; i < selectionCategory; i++)
	{
		currentSphi[i] = getSphi(i, false);
		currentMphi[i] = -((currentSphi[i] * currentSphi[i]) / 2);
		proposedSphi[i] = getSphi(i, true);
		proposedMphi[i] = -((proposedSphi[i] * proposedSphi[i]) / 2);
		// take the jacobian into account for the non-linear transformation from logN to N distribution
		lpr -= (std::log(currentSphi[i]) - std::log(proposedSphi[i]));

		// take prior into account
		//TODO(Cedric): make sure you can control that prior from R
		//lpr -= Parameter::densityNorm(currentSphi[i], 1.0, 0.1, true) - Parameter::densityNorm(proposedSphi[i], 1.0, 0.1, true);
	}

	if (withPhi) {
		// one for each Aphi, and one for Sphi
		logProbabilityRatio.resize(getNumPhiGroupings()+1);
	}
	else {
		logProbabilityRatio.resize(1);
	}
#ifndef __APPLE__
//#pragma omp parallel for reduction(+:lpr)
#endif
	for (int i = 0; i < genome.getGenomeSize(); i++)
	{
		unsigned mixture = getMixtureAssignment(i);
		mixture = getSynthesisRateCategory(mixture);
		double phi = getSynthesisRate(i, mixture, false);
		if (!std::isfinite(phi))
		{
#ifndef STANDALONE
			Rf_error("Phi value for gene %d is not finite (%f)!", i, phi);
#else
			std::cerr << "phi " << i << " not finite! " << phi << "\n";
#endif
		}
		lpr += Parameter::densityLogNorm(phi, proposedMphi[mixture], proposedSphi[mixture], true)
				- Parameter::densityLogNorm(phi, currentMphi[mixture], currentSphi[mixture], true);
	}

	// TODO: USE CONSTANTS INSTEAD OF 0
	logProbabilityRatio[0] = lpr;

	if (withPhi)
	{
		for (unsigned i = 0; i < parameter->getNumObservedPhiSets(); i++) {
			lpr = 0.0;
			double Aphi = getAphi(i, false);
			double Aphi_proposed = getAphi(i, true);
			double Sepsilon = getSepsilon(i);
#ifndef __APPLE__
//#pragma omp parallel for reduction(+:lpr)
#endif
			for (int j = 0; j < genome.getGenomeSize(); j++) {
				unsigned mixtureAssignment = getMixtureAssignment(j);
				mixtureAssignment = getSynthesisRateCategory(mixtureAssignment);
				double logphi = std::log(getSynthesisRate(j, mixtureAssignment, false));
				double obsPhi = genome.getGene(j).getObservedSynthesisRate(i);
				if (obsPhi > -1.0) {
					double logobsPhi = std::log(obsPhi);
					double proposed = Parameter::densityNorm(logobsPhi, logphi + Aphi_proposed, Sepsilon, true);
					double current = Parameter::densityNorm(logobsPhi, logphi + Aphi, Sepsilon, true);
					lpr += proposed - current;
				}
			}
			logProbabilityRatio[i+1] = lpr;
		}
	}
}

void ROCModel::updateGibbsSampledHyperParameters(Genome &genome)
{
	// TODO: Fix this for any numbers of phi values
	if (withPhi) {
		double shape = ((double)genome.getGenomeSize() - 1.0) / 2.0;
		for (unsigned i = 0; i < parameter->getNumObservedPhiSets(); i++) {
			double rate = 0.0;
			unsigned mixtureAssignment;
			double aphi = getAphi(i);
			for (unsigned j = 0; j < genome.getGenomeSize(); j++) {
				mixtureAssignment = getMixtureAssignment(i);
				double obsPhi = genome.getGene(j).getObservedSynthesisRate(i);
				if (obsPhi > -1.0) {
					double sum = std::log(obsPhi) - aphi - std::log(getSynthesisRate(j, mixtureAssignment, false));
					rate += sum * sum;
				}
			}
			rate /= 2;
			// TODO(Cedric): rate has to be 1/rate? see relatioship gamma to inv-gamma
			double rand = parameter->randGamma(shape, rate);
			parameter->setSepsilon(i, std::sqrt(1 / rand));
		}
	}
}

void ROCModel::proposeHyperParameters()
{
	parameter->proposeSphi();
	if (withPhi) {
		parameter->proposeAphi();
	}
}

void ROCModel::adaptHyperParameterProposalWidths(unsigned adaptiveWidth)
{
	adaptSphiProposalWidth(adaptiveWidth);
	if (withPhi) {
		adaptAphiProposalWidth(adaptiveWidth);
	}
}

void ROCModel::updateAllHyperParameter()
{
	updateSphi();
	for (unsigned i = 0; i < parameter->getNumObservedPhiSets(); i++) {
		updateAphi(i);
	}
}

void ROCModel::updateHyperParameter(unsigned hp)
{
	if (hp == 0) {
		updateSphi();
	}
	else {
		updateAphi(hp - 1);
	}
}

void ROCModel::updateHyperParameterTraces(unsigned sample)
{
	updateSphiTrace(sample);
	for (unsigned i = 0; i < parameter->getNumObservedPhiSets(); i++) {
		updateAphiTrace(i, sample);
		updateSepsilonTrace(i, sample);
	}
}

void ROCModel::updateTracesWithInitialValues(Genome &genome)
{
	std::vector <std::string> groupList = parameter->getGroupList();

	for (unsigned i = 0; i < genome.getGenomeSize(); i++)
	{
		parameter->updateSynthesisRateTrace(0, i);
		parameter->updateMixtureAssignmentTrace(0, i);
	}

	for (unsigned i = 0; i < groupList.size(); i++)
	{
		std::array <unsigned, 2> codonRange = SequenceSummary::AAToCodonRange(groupList[i], true);
		for (unsigned j = codonRange[0]; j < codonRange[1]; j++)
		{
			std::string codon = SequenceSummary::indexToCodon(j, true);
			parameter->updateCodonSpecificParameterTrace(0, codon);
		}
	}
}
