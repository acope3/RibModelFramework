#include "include/base/Model.h"

Model::Model()
{
  //ctor
  parameter = 0;
  withPhi = false;
  fix_sEpsilon = false;
  parameter_types = {};
}

// TODO: Rule of Three dictates we may need a copy assignment operator as well (operator=)

Model::~Model()
{
//dtor
}




//Cedric: This functions will replace calculateMutationPrior in ROC/FONSE model and allows us to more generally use priors on codon specific parameters.
//			We have to first change how current and proposed csp values are stored to move the function getParameterForCategory up into the base parameter class.

double Model::calculatePriorForCodonSpecificParam(Parameter *parameter, std::string grouping, unsigned paramType, bool proposed)
{
	unsigned numCodons = SequenceSummary::GetNumCodonsForAA(grouping, true); // TODO(Cedric): rename getNumCodonsForGrouping and have it return 1 if grouping is a codon to make it applicable for RFP
	double parameterValues[5];

	double priorValue = 0.0;

	unsigned numCat = parameter->getNumMutationCategories(); // TODO(Cedric): rename this function getNumCSPone or something like that
	double prior_sd = parameter->getCodonSpecificPriorStdDev(paramType);
	for (unsigned i = 0u; i < numCat; i++)
	{
		//parameter->getParameterForCategory(i, paramType, grouping, proposed, mutation); // TODO(Cedric): we have to change how csp are stored first!
		for (unsigned k = 0u; k < numCodons; k++)
		{
			priorValue += Parameter::densityNorm(parameterValues[k], 0.0, prior_sd, true);
		}
	}
	return priorValue;
}

std::vector<std::string> Model::getParameterTypeList()
{
	return parameter_types;
}


bool Model::isShared(std::string csp_parameters)
{
	return false;
}

void Model::fillMatrices(Genome& genome)
{
  //do nothing
}
void Model::clearMatrices()
{
	//do nothing
}



//TO DO, Alex: As of now, seems like base Model object cannot see the parameter object initialized in the derived model object.
// To make withPhi fittings work for all models, I am unable to get this to work. We have a less than ideal where we have the same code in all of the model classes
//, which call functions initialized in in the base Parameter class


double Model::getNoiseOffset(unsigned index, bool proposed)
{
	return parameter->getNoiseOffset(index, proposed);
}


double Model::getObservedSynthesisNoise(unsigned index)
{
	return parameter->getObservedSynthesisNoise(index);
}


double Model::getCurrentNoiseOffsetProposalWidth(unsigned index)
{
	return parameter->getCurrentNoiseOffsetProposalWidth(index);
}


void Model::updateNoiseOffset(unsigned index)
{
	parameter->updateNoiseOffset(index);
}


void Model::updateNoiseOffsetTrace(unsigned sample)
{
	parameter->updateNoiseOffsetTraces(sample);
}


void Model::updateObservedSynthesisNoiseTrace(unsigned sample)
{
	parameter->updateObservedSynthesisNoiseTraces(sample);
}


void Model::adaptNoiseOffsetProposalWidth(unsigned adaptiveWidth, bool adapt)
{
	parameter->adaptNoiseOffsetProposalWidth(adaptiveWidth, adapt);
}



void Model::updateGibbsSampledHyperParameters(Genome &genome)
{
	parameter->updateGibbsSampledHyperParameters(genome, withPhi, fix_sEpsilon);
}
