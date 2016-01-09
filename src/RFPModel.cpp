#include "include/RFP/RFPModel.h"

//R runs only
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif


//--------------------------------------------------//
//----------- Constructors & Destructors ---------- //
//--------------------------------------------------//


RFPModel::RFPModel() : Model()
{
	parameter = 0;
	//ctor
}


RFPModel::~RFPModel()
{
	//dtor
	//TODO: call Parent's deconstructor
}





double RFPModel::calculateLogLikelihoodPerCodonPerGene(double currAlpha, double currLambdaPrime,
													   unsigned currRFPObserved, unsigned currNumCodonsInMRNA, double phiValue)
{
	double logLikelihood = ((std::lgamma((currNumCodonsInMRNA * currAlpha) + currRFPObserved)) - (std::lgamma(currNumCodonsInMRNA * currAlpha)))
						   + (currRFPObserved * (std::log(phiValue) - std::log(currLambdaPrime + phiValue))) + ((currNumCodonsInMRNA * currAlpha) * (std::log(currLambdaPrime) -
																																					 std::log(currLambdaPrime + phiValue)));

	return logLikelihood;
}





//------------------------------------------------//
//---------- Likelihood Ratio Functions ----------//
//------------------------------------------------//


void RFPModel::calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k, double* logProbabilityRatio)
{

	double logLikelihood = 0.0;
	double logLikelihood_proposed = 0.0;

	unsigned alphaCategory = parameter->getMutationCategory(k);
	unsigned lambdaPrimeCategory = parameter->getSelectionCategory(k);
	unsigned synthesisRateCategory = parameter->getSynthesisRateCategory(k);


	double phiValue = parameter->getSynthesisRate(geneIndex, synthesisRateCategory, false);
	double phiValue_proposed = parameter->getSynthesisRate(geneIndex, synthesisRateCategory, true);

#ifndef __APPLE__
#pragma omp parallel for reduction(+:logLikelihood,logLikelihood_proposed)
#endif
	for (int index = 0; index < getGroupListSize(); index++) //number of codons, without the stop codons
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

	double stdDevSynthesisRate = parameter->getStdDevSynthesisRate(false);
	double logPhiProbability = Parameter::densityLogNorm(phiValue, (-(stdDevSynthesisRate * stdDevSynthesisRate) / 2), stdDevSynthesisRate, true);
	double logPhiProbability_proposed = Parameter::densityLogNorm(phiValue_proposed, (-(stdDevSynthesisRate * stdDevSynthesisRate) / 2), stdDevSynthesisRate, true);
	double currentLogLikelihood = (logLikelihood + logPhiProbability);
	double proposedLogLikelihood = (logLikelihood_proposed + logPhiProbability_proposed);

	logProbabilityRatio[0] = (proposedLogLikelihood - currentLogLikelihood) - (std::log(phiValue) - std::log(phiValue_proposed));
	logProbabilityRatio[1] = currentLogLikelihood - std::log(phiValue_proposed);
	logProbabilityRatio[2] = proposedLogLikelihood - std::log(phiValue);
	logProbabilityRatio[3] = currentLogLikelihood;
	logProbabilityRatio[4] = proposedLogLikelihood;
}


void RFPModel::calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome, double& logAcceptanceRatioForAllMixtures)
{
	double logLikelihood = 0.0;
	double logLikelihood_proposed = 0.0;
	Gene *gene;
	unsigned index = SequenceSummary::codonToIndex(grouping);


#ifndef __APPLE__
#pragma omp parallel for private(gene) reduction(+:logLikelihood,logLikelihood_proposed)
#endif
	for (int i = 0u; i < genome.getGenomeSize(); i++)
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


		logLikelihood += calculateLogLikelihoodPerCodonPerGene(currAlpha, currLambdaPrime, currRFPObserved, currNumCodonsInMRNA, phiValue);
		logLikelihood_proposed += calculateLogLikelihoodPerCodonPerGene(propAlpha, propLambdaPrime, currRFPObserved, currNumCodonsInMRNA, phiValue);
	}
	logAcceptanceRatioForAllMixtures = logLikelihood_proposed - logLikelihood;
}


void RFPModel::calculateLogLikelihoodRatioForHyperParameters(Genome &genome, unsigned iteration, std::vector <double> & logProbabilityRatio)
{

	double lpr = 0.0; // this variable is only needed because OpenMP doesn't allow variables in reduction clause to be reference

	unsigned selectionCategory = getNumSynthesisRateCategories();
	std::vector<double> currentStdDevSynthesisRate(selectionCategory, 0.0);
	std::vector<double> currentMphi(selectionCategory, 0.0);
	std::vector<double> proposedStdDevSynthesisRate(selectionCategory, 0.0);
	std::vector<double> proposedMphi(selectionCategory, 0.0);
	for(unsigned i = 0u; i < selectionCategory; i++)
	{
		currentStdDevSynthesisRate[i] = getStdDevSynthesisRate(i, false);
		currentMphi[i] = -((currentStdDevSynthesisRate[i] * currentStdDevSynthesisRate[i]) / 2);
		proposedStdDevSynthesisRate[i] = getStdDevSynthesisRate(i, true);
		proposedMphi[i] = -((proposedStdDevSynthesisRate[i] * proposedStdDevSynthesisRate[i]) / 2);
		// take the jacobian into account for the non-linear transformation from logN to N distribution
		lpr -= (std::log(currentStdDevSynthesisRate[i]) - std::log(proposedStdDevSynthesisRate[i]));
		// take prior into account
		//TODO(Cedric): make sure you can control that prior from R
		lpr -= Parameter::densityNorm(currentStdDevSynthesisRate[i], 1.0, 0.1, true) - Parameter::densityNorm(proposedStdDevSynthesisRate[i], 1.0, 0.1, true);
	}




	logProbabilityRatio.resize(1);
#ifndef __APPLE__
#pragma omp parallel for reduction(+:lpr)
#endif
	for (int i = 0u; i < genome.getGenomeSize(); i++)
	{
		unsigned mixture = getMixtureAssignment(i);
		mixture = getSynthesisRateCategory(mixture);
		double phi = getSynthesisRate(i, mixture, false);
		if (i == 0) {
	//	std::cout <<"proposed: " << Parameter::densityLogNorm(phi, proposedMPhi[mixture], proposedStdDevSynthesisRate[mixture], false) <<"\n";
		//std::cout <<"current: " << Parameter::densityLogNorm(phi, currentMPhi, currentStdDevSynthesisRate, false) <<"\n";
		}
		lpr += Parameter::densityLogNorm(phi, proposedMphi[mixture], proposedStdDevSynthesisRate[mixture], true) -
				Parameter::densityLogNorm(phi, currentMphi[mixture], currentStdDevSynthesisRate[mixture], true);
		//std::cout <<"LPR: " << lpr <<"\n";
	}

	logProbabilityRatio[0] = lpr;
}





//----------------------------------------------------------//
//---------- Initialization and Restart Functions ----------//
//----------------------------------------------------------//


void RFPModel::initTraces(unsigned samples, unsigned num_genes)
{
	parameter->initAllTraces(samples, num_genes);
}


void RFPModel::writeRestartFile(std::string filename)
{
	return parameter->writeEntireRestartFile(filename);
}






//----------------------------------------//
//---------- Category Functions ----------//
//----------------------------------------//


double RFPModel::getCategoryProbability(unsigned i)
{
	return parameter->getCategoryProbability(i);
}


unsigned RFPModel::getMutationCategory(unsigned mixture)
{
	return parameter->getMutationCategory(mixture);
}


unsigned RFPModel::getSelectionCategory(unsigned mixture)
{
	return parameter->getSelectionCategory(mixture);
}


unsigned RFPModel::getSynthesisRateCategory(unsigned mixture)
{
	return parameter->getSynthesisRateCategory(mixture);
}


std::vector<unsigned> RFPModel::getMixtureElementsOfSelectionCategory(unsigned k)
{
	return parameter->getMixtureElementsOfSelectionCategory(k);
}





//------------------------------------------//
//---------- Group List Functions ----------//
//------------------------------------------//


unsigned RFPModel::getGroupListSize()
{
	return parameter->getGroupListSize();
}


std::string RFPModel::getGrouping(unsigned index)
{
	return parameter->getGrouping(index);
}





//---------------------------------------------------//
//---------- stdDevSynthesisRate Functions ----------//
//---------------------------------------------------//


double RFPModel::getStdDevSynthesisRate(unsigned selectionCategory, bool proposed)
{
	return parameter->getStdDevSynthesisRate(selectionCategory, proposed);
}


double RFPModel::getCurrentStdDevSynthesisRateProposalWidth()
{
	return parameter->getCurrentStdDevSynthesisRateProposalWidth();
}


void RFPModel::updateStdDevSynthesisRate()
{
	parameter->updateStdDevSynthesisRate();
}





//----------------------------------------------//
//---------- Synthesis Rate Functions ----------//
//----------------------------------------------//


double RFPModel::getSynthesisRate(unsigned index, unsigned mixture, bool proposed)
{
	return parameter->getSynthesisRate(index, mixture, proposed);
}


void RFPModel::updateSynthesisRate(unsigned i, unsigned k)
{
	parameter->updateSynthesisRate(i, k);
}





//-----------------------------------------//
//---------- Iteration Functions ----------//
//-----------------------------------------//


unsigned RFPModel::getLastIteration()
{
	return parameter->getLastIteration();
}


void RFPModel::setLastIteration(unsigned iteration)
{
	parameter->setLastIteration(iteration);
}





//-------------------------------------//
//---------- Trace Functions ----------//
//-------------------------------------//


void RFPModel::updateStdDevSynthesisRateTrace(unsigned sample)
{
	parameter->updateStdDevSynthesisRateTrace(sample);
}


void RFPModel::updateSynthesisRateTrace(unsigned sample, unsigned i)
{
	parameter->updateSynthesisRateTrace(sample, i);
}


void RFPModel::updateMixtureAssignmentTrace(unsigned sample, unsigned i)
{
	parameter->updateMixtureAssignmentTrace(sample, i);
}


void RFPModel::updateMixtureProbabilitiesTrace(unsigned sample)
{
	parameter->updateMixtureProbabilitiesTrace(sample);
}


void RFPModel::updateCodonSpecificParameterTrace(unsigned sample, std::string codon)
{
	parameter->updateCodonSpecificParameterTrace(sample, codon);
}


void RFPModel::updateHyperParameterTraces(unsigned sample)
{
	updateStdDevSynthesisRateTrace(sample);
}


void RFPModel::updateTracesWithInitialValues(Genome & genome)
{
	return; //TODO: fill in
}





//----------------------------------------------//
//---------- Adaptive Width Functions ----------//
//----------------------------------------------//


void RFPModel::adaptStdDevSynthesisRateProposalWidth(unsigned adaptiveWidth)
{
	parameter->adaptStdDevSynthesisRateProposalWidth(adaptiveWidth);
}


void RFPModel::adaptSynthesisRateProposalWidth(unsigned adaptiveWidth)
{
	parameter->adaptSynthesisRateProposalWidth(adaptiveWidth);
}


void RFPModel::adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth)
{
	parameter->adaptCodonSpecificParameterProposalWidth(adaptiveWidth);
}


void RFPModel::adaptHyperParameterProposalWidths(unsigned adaptiveWidth)
{
	adaptStdDevSynthesisRateProposalWidth(adaptiveWidth);
}





//-------------------------------------//
//---------- Other Functions ----------//
//-------------------------------------//


void RFPModel::proposeCodonSpecificParameter()
{
	parameter->proposeCodonSpecificParameter();
}


void RFPModel::proposeHyperParameters()
{
	parameter->proposeStdDevSynthesisRate();
}


void RFPModel::proposeSynthesisRateLevels()
{
	parameter->proposeSynthesisRateLevels();
}


unsigned RFPModel::getNumPhiGroupings()
{
	return parameter->getNumObservedPhiSets();
}


unsigned RFPModel::getMixtureAssignment(unsigned index)
{
	return parameter->getMixtureAssignment(index);
}


unsigned RFPModel::getNumMixtureElements()
{
	return parameter->getNumMixtureElements();
}


unsigned RFPModel::getNumSynthesisRateCategories()
{
	return parameter->getNumSynthesisRateCategories();
}


void RFPModel::setNumPhiGroupings(unsigned value)
{
	parameter->setNumObservedPhiSets(value);
}


void RFPModel::setMixtureAssignment(unsigned i, unsigned catOfGene)
{
	parameter->setMixtureAssignment(i, catOfGene);
}


void RFPModel::setCategoryProbability(unsigned mixture, double value)
{
	parameter->setCategoryProbability(mixture, value);
}


void RFPModel::updateCodonSpecificParameter(std::string aa)
{
	parameter->updateCodonSpecificParameter(aa);
}


void RFPModel::updateGibbsSampledHyperParameters(Genome &genome)
{
	//TODO: fill in
}


void RFPModel::updateAllHyperParameter()
{
	updateStdDevSynthesisRate();
}


void RFPModel::updateHyperParameter(unsigned hp)
{
	// NOTE: when adding additional hyper parameter, also add to updateAllHyperParameter()
	switch (hp) {
		case 0:
			updateStdDevSynthesisRate();
			break;
	}
}


void RFPModel::simulateGenome(Genome &genome)
{
	for (unsigned geneIndex = 0; geneIndex < genome.getGenomeSize(); geneIndex++)
	{
		unsigned mixtureElement = getMixtureAssignment(geneIndex);
		Gene gene = genome.getGene(geneIndex);
		double phi = parameter -> getSynthesisRate(geneIndex, mixtureElement, false);
		Gene tmpGene = gene;
		for (unsigned codonIndex = 0; codonIndex < 61; codonIndex++)
		{
			std::string codon = SequenceSummary::codonArray[codonIndex];
			unsigned alphaCat = parameter -> getMutationCategory(mixtureElement);
			unsigned lambdaPrimeCat = parameter -> getSelectionCategory(mixtureElement);

			double alpha = getParameterForCategory(alphaCat, RFPParameter::alp, codon, false);
			double lambdaPrime = getParameterForCategory(lambdaPrimeCat, RFPParameter::lmPri, codon, false);

			double alphaPrime = alpha * gene.geneData.getCodonCountForCodon(codon);

#ifndef STANDALONE
				RNGScope scope;
				NumericVector xx(1);
				xx = rgamma(1, alphaPrime, 1.0/lambdaPrime);
				xx = rpois(1, xx[0] * phi);
				tmpGene.geneData.setRFPObserved(codonIndex, xx[0]);
			#else
			std::gamma_distribution<double> GDistribution(alphaPrime,1.0/lambdaPrime);
			double tmp = GDistribution(Parameter::generator);
			std::poisson_distribution<unsigned> PDistribution(phi * tmp);
			unsigned simulatedValue = PDistribution(Parameter::generator);
			tmpGene.geneData.setRFPObserved(codonIndex, simulatedValue);
#endif
		}
		genome.addGene(tmpGene, true);
	}
}


void RFPModel::printHyperParameters()
{
	for(unsigned i = 0u; i < getNumSynthesisRateCategories(); i++)
	{
		std::cout << "stdDevSynthesisRate posterior estimate for selection category " << i << ": " << parameter -> getStdDevSynthesisRate(i) << std::endl;
	}
	std::cout << "\t current stdDevSynthesisRate proposal width: " << getCurrentStdDevSynthesisRateProposalWidth() << std::endl;
}


void RFPModel::setParameter(RFPParameter &_parameter)
{
	parameter = &_parameter;
}


double RFPModel::calculateAllPriors()
{
	return 0.0; //TODO(Cedric): implement me, see ROCModel
}


double RFPModel::getParameterForCategory(unsigned category, unsigned param, std::string codon, bool proposal)
{
	return parameter->getParameterForCategory(category, param, codon, proposal);
}








