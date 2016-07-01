#include "include/base/Trace.h"
#include "include/SequenceSummary.h"


#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif



//--------------------------------------------------//
//----------- Constructors & Destructors -----------//
//--------------------------------------------------//
 
 
Trace::Trace()
{
	categories = 0;
	numCodonSpecificParamTypes = 2;
	codonSpecificParameterTrace.resize(numCodonSpecificParamTypes);
	// TODO: fill this
}


Trace::~Trace()
{
	//dtor
}


Trace::Trace(unsigned _numCodonSpecificParamTypes)
{
	categories = 0;
	numCodonSpecificParamTypes = _numCodonSpecificParamTypes;
	codonSpecificParameterTrace.resize(numCodonSpecificParamTypes);
}


//-----------------------------------------------------//
//----------Private Initialization Functions ----------//
//-----------------------------------------------------//


void Trace::initializeSharedTraces(unsigned samples, unsigned num_genes, unsigned numSelectionCategories,
	unsigned numMixtures, std::vector<mixtureDefinition> &_categories, unsigned maxGrouping)
{
	my_print("maxGrouping: %\n", maxGrouping);

	//numSelectionCategories always == numSynthesisRateCategories, so only one is passed in for convenience
	
	initStdDevSynthesisRateTrace(numSelectionCategories, samples);
	initSynthesisRateAcceptanceRatioTrace(num_genes, numSelectionCategories);
	codonSpecificAcceptanceRatioTrace.resize(maxGrouping);
	initSynthesisRateTrace(samples, num_genes, numSelectionCategories);
	initMixtureAssignmentTrace(samples, num_genes);
	initMixtureProbabilitiesTrace(samples, numMixtures);

	categories = &_categories;
}


void Trace::initStdDevSynthesisRateTrace(unsigned numSelectionCategories, unsigned samples)
{
	stdDevSynthesisRateTrace.resize(numSelectionCategories);
	for (unsigned i = 0u; i < numSelectionCategories; i++)
	{
		std::vector<double> temp(samples, 0.0);
		stdDevSynthesisRateTrace[i] = temp;
	}
}


void Trace::initSynthesisRateAcceptanceRatioTrace(unsigned num_genes, unsigned numSynthesisRateCategories)
{
	synthesisRateAcceptanceRatioTrace.resize(numSynthesisRateCategories);
	for (unsigned category = 0; category < numSynthesisRateCategories; category++)
	{
		synthesisRateAcceptanceRatioTrace[category].resize(num_genes);
	}
	//NOTE: this is not sized to samples because push_back takes care of the initialization
}


void Trace::initSynthesisRateTrace(unsigned samples, unsigned num_genes, unsigned numSynthesisRateCategories)
{
	synthesisRateTrace.resize(numSynthesisRateCategories);
	for (unsigned category = 0; category < numSynthesisRateCategories; category++)
	{
		synthesisRateTrace[category].resize(num_genes);
		for (unsigned i = 0; i < num_genes; i++)
		{
			std::vector<double> tempExpr(samples, 0.0);
			synthesisRateTrace[category][i] = tempExpr;
		}
	}
}


void Trace::initMixtureAssignmentTrace(unsigned samples, unsigned num_genes)
{
	mixtureAssignmentTrace.resize(num_genes);
	for (unsigned i = 0u; i < num_genes; i++)
	{
		mixtureAssignmentTrace[i].resize(samples);
	}
}


void Trace::initMixtureProbabilitiesTrace(unsigned samples, unsigned numMixtures)
{
	mixtureProbabilitiesTrace.resize(numMixtures);
	for (unsigned i = 0u; i < numMixtures; i++)
	{
		mixtureProbabilitiesTrace[i].resize(samples, 0.0);
	}
}


void Trace::initCodonSpecificParameterTrace(unsigned samples, unsigned numCategories, unsigned numParam, unsigned paramType)
{
	std::vector <std::vector <std::vector <double>>> tmp;
	tmp.resize(numCategories);
	for (unsigned category = 0; category < numCategories; category++)
	{
		tmp[category].resize(numParam);
		for (unsigned i = 0; i < numParam; i++)
		{
			std::vector <double> temp(samples, 0.0);
			tmp[category][i] = temp;
		}
	}


	//TODO: R output for error message here
	codonSpecificParameterTrace[paramType] = tmp;
	/*
	switch (paramType) {
	case 0:
		codonSpecificParameterTraceOne = tmp;
		break;
	case 1:
		codonSpecificParameterTraceTwo = tmp;
		break;
	default:
		my_printError("ERROR: Invalid paramType given, codon specific parameter trace not initialized.\n");
		break;
	}
	*/
}





//----------------------------------//
//---------- ROC Specific ----------//
//----------------------------------//

void Trace::initSynthesisOffsetTrace(unsigned samples, unsigned numPhiGroupings)
{
	synthesisOffsetTrace.resize(numPhiGroupings);
	for (unsigned i = 0; i < numPhiGroupings; i++)
	{
		synthesisOffsetTrace[i].resize(samples);
	}

	synthesisOffsetAcceptanceRatioTrace.resize(numPhiGroupings);
}


void Trace::initObservedSynthesisNoiseTrace(unsigned samples, unsigned numPhiGroupings)
{
	observedSynthesisNoiseTrace.resize(numPhiGroupings);
	for (unsigned i = 0; i < numPhiGroupings; i++)
	{
		observedSynthesisNoiseTrace[i].resize(samples);
	}
}





//----------------------------------------------------//
//---------- Model Initialization Functions ----------//
//----------------------------------------------------//

/* Note 1) -- on these Model Initialization Functions' initCodonSpecificParameterTrace functions
 * The last argument specifies the codon specific parameter type (dM and dEta for ROC)
 * You can check Parameter.cpp to check what the values of dM and dEta are.
 * The constants aren't used here because they are not available to the Trace object.
*/

void Trace::initializeRFPTrace(unsigned samples, unsigned num_genes, unsigned numAlphaCategories,
	unsigned numLambdaPrimeCategories, unsigned numParam, unsigned numMixtures,
	std::vector<mixtureDefinition> &_categories, unsigned maxGrouping)
{
	initializeSharedTraces(samples, num_genes, numLambdaPrimeCategories, numMixtures,
		_categories, maxGrouping);

	// See Note 1) above.
	initCodonSpecificParameterTrace(samples, numAlphaCategories,  numParam, 0u); // alp
	initCodonSpecificParameterTrace(samples, numLambdaPrimeCategories, numParam, 1u); // lmPri
}


void Trace::initializeROCTrace(unsigned samples, unsigned num_genes, unsigned numMutationCategories,
	unsigned numSelectionCategories, unsigned numParam, unsigned numMixtures,
	std::vector<mixtureDefinition> &_categories, unsigned maxGrouping, unsigned numObservedPhiSets)
{
	initializeSharedTraces(samples, num_genes, numSelectionCategories, numMixtures, _categories, maxGrouping);

	// See Note 1) above.
	initCodonSpecificParameterTrace(samples, numMutationCategories, numParam, 0u); // dM
	initCodonSpecificParameterTrace(samples, numSelectionCategories, numParam, 1u); // dEta

	initSynthesisOffsetTrace(samples, numObservedPhiSets);
	initObservedSynthesisNoiseTrace(samples, numObservedPhiSets);
}


void Trace::initializeFONSETrace(unsigned samples, unsigned num_genes, unsigned numMutationCategories,
	unsigned numSelectionCategories, unsigned numParam, unsigned numMixtures,
	std::vector<mixtureDefinition> &_categories, unsigned maxGrouping)
{
	initializeSharedTraces(samples, num_genes, numSelectionCategories, numMixtures,
		 _categories, maxGrouping);

	// See Note 1) above.
	initCodonSpecificParameterTrace(samples, numMutationCategories, numParam, 0u); // dM
	initCodonSpecificParameterTrace(samples, numSelectionCategories, numParam, 1u); // dOmega
}


void Trace::initializePANSETrace(unsigned samples, unsigned num_genes, unsigned numAlphaCategories,
	unsigned numLambdaPrimeCategories, unsigned numParam, unsigned numMixtures,
	std::vector<mixtureDefinition> &_categories, unsigned maxGrouping)
{
	initializeSharedTraces(samples, num_genes, numLambdaPrimeCategories, numMixtures,
		_categories, maxGrouping);

	// See Note 1) above.
	initCodonSpecificParameterTrace(samples, numAlphaCategories,  numParam, 0u); //alp
	initCodonSpecificParameterTrace(samples, numLambdaPrimeCategories, numParam, 1u); //lmPri
}


//--------------------------------------//
//---------- Getter Functions ----------//
//--------------------------------------//


std::vector<double> Trace::getStdDevSynthesisRateTrace(unsigned selectionCategory) 
{ 
	return stdDevSynthesisRateTrace[selectionCategory]; 
}


std::vector<double> Trace::getExpectedSynthesisRateTrace()
{
	unsigned numGenes = synthesisRateTrace[0].size(); //number of genes
	unsigned samples = synthesisRateTrace[0][0].size(); //number of samples
	std::vector<double> RV(samples, 0.0);
	for (unsigned sample = 0; sample < samples; sample++)
	{
		for (unsigned geneIndex = 0; geneIndex < numGenes; geneIndex++)
		{
			unsigned mixtureElement = mixtureAssignmentTrace[geneIndex][sample];
			unsigned category = getSynthesisRateCategory(mixtureElement);
			RV[sample] += synthesisRateTrace[category][geneIndex][sample];
		}
		RV[sample] /= numGenes;
	}
	return RV;
}


std::vector<double> Trace::getStdDevSynthesisRateAcceptanceRatioTrace()
{
	return stdDevSynthesisRateAcceptanceRatioTrace;
}


std::vector<std::vector<std::vector<double>>> Trace::getSynthesisRateTrace()
{
	return synthesisRateTrace;
}


std::vector<double> Trace::getSynthesisRateAcceptanceRatioTraceByMixtureElementForGene(unsigned mixtureElement,
	unsigned geneIndex)
{
	unsigned category = getSynthesisRateCategory(mixtureElement);
	return synthesisRateAcceptanceRatioTrace[category][geneIndex];
}


std::vector<std::vector<std::vector<double>>> Trace::getSynthesisRateAcceptanceRatioTrace()
{
	return synthesisRateAcceptanceRatioTrace;
}


std::vector<double> Trace::getCodonSpecificAcceptanceRatioTraceForAA(std::string aa)
{
	aa[0] = (char)std::toupper(aa[0]);
	unsigned aaIndex = SequenceSummary::aaToIndex.find(aa)->second;
	return codonSpecificAcceptanceRatioTrace[aaIndex];
}


std::vector<double> Trace::getSynthesisRateTraceForGene(unsigned geneIndex)
{
	unsigned traceLength = synthesisRateTrace[0][0].size();

	std::vector<double> returnVector(traceLength, 0.0);
	for (unsigned i = 0u; i < traceLength; i++)
	{
		unsigned mixtureElement = mixtureAssignmentTrace[geneIndex][i];
		unsigned category = getSynthesisRateCategory(mixtureElement);
		returnVector[i] = synthesisRateTrace[category][geneIndex][i];
	}
	return returnVector;
}


std::vector<double> Trace::getSynthesisRateTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex)
{
	unsigned category = getSynthesisRateCategory(mixtureElement);
	return synthesisRateTrace[category][geneIndex];
}


std::vector<unsigned> Trace::getMixtureAssignmentTraceForGene(unsigned geneIndex)
{
	return mixtureAssignmentTrace[geneIndex];
}


std::vector<double> Trace::getMixtureProbabilitiesTraceForMixture(unsigned mixtureIndex)
{
	return mixtureProbabilitiesTrace[mixtureIndex];
}


std::vector<std::vector<unsigned>> Trace::getMixtureAssignmentTrace()
{
	return mixtureAssignmentTrace;
}


std::vector<std::vector<double>> Trace::getMixtureProbabilitiesTrace()
{
	return mixtureProbabilitiesTrace;
}


std::vector<std::vector<double>> Trace::getCodonSpecificAcceptanceRatioTrace()
{
	return codonSpecificAcceptanceRatioTrace;
}


unsigned Trace::getSynthesisRateCategory(unsigned mixtureElement)
{
	return categories->at(mixtureElement).delEta;
}


std::vector<std::vector<std::vector<std::vector<double>>>>* Trace::getCodonSpecificParameterTrace() 
{
	return &codonSpecificParameterTrace;
}


//----------------------------------//
//---------- ROC Specific ----------//
//----------------------------------//

std::vector<double> Trace::getCodonSpecificParameterTraceByMixtureElementForCodon(unsigned mixtureElement,
	std::string& codon, unsigned paramType, bool withoutReference)
{
	std::vector <double> rv;
	unsigned codonIndex = SequenceSummary::codonToIndex(codon, withoutReference);
	unsigned category = getCodonSpecificCategory(mixtureElement, paramType);
	rv = codonSpecificParameterTrace[paramType][category][codonIndex];
	/*
	switch (paramType) {
	case 0:
		rv = codonSpecificParameterTraceOne[category][codonIndex];
		break;
	case 1:
		rv = codonSpecificParameterTraceTwo[category][codonIndex];
		break;
	default:
		my_printError("ERROR: Unknown Parameter type in getCodonSpecificParameterTraceByMixtureElementForCodon\n");
		break;
	}
	*/
	return rv;
}


std::vector<double> Trace::getSynthesisOffsetTrace(unsigned index)
{
	return synthesisOffsetTrace[index];
}


std::vector<double> Trace::getSynthesisOffsetAcceptanceRatioTraceForIndex(unsigned index)
{
	return synthesisOffsetAcceptanceRatioTrace[index];
}


std::vector<double> Trace::getObservedSynthesisNoiseTrace(unsigned index)
{
	return observedSynthesisNoiseTrace[index];
}


std::vector<std::vector<std::vector<double>>> Trace::getCodonSpecificParameterTraceByParamType(unsigned paramType)
{
	/*
	std::vector<std::vector<std::vector<double>>> rv;
	switch (paramType) {
	case 0:
		rv = codonSpecificParameterTraceOne;
		break;
	case 1:
		rv = codonSpecificParameterTraceTwo;
		break;
	default:
		my_printError("ERROR: Unknown Parameter type in getCodonSpecificParameterTraceByParamType\n";
		break;
	}
	*/
	return codonSpecificParameterTrace[paramType];
}


std::vector<std::vector<double>> Trace::getSynthesisOffsetAcceptanceRatioTrace()
{
	return synthesisOffsetAcceptanceRatioTrace;
}


unsigned Trace::getCodonSpecificCategory(unsigned mixtureElement, unsigned paramType)
{
	unsigned rv = 0;
	switch (paramType) {
	case 0:
		rv = categories->at(mixtureElement).delM;
		break;
	case 1:
		rv = categories->at(mixtureElement).delEta;
		break;
	default:
		my_printError("ERROR: Unknown parameter type in getCodonSpecificCategory\n");
		break;
	}
	return rv;
}





//--------------------------------------//
//---------- Update Functions ----------//
//--------------------------------------//


void Trace::updateStdDevSynthesisRateTrace(unsigned sample, double stdDevSynthesisRate, unsigned synthesisRateCategory)
{
	stdDevSynthesisRateTrace[synthesisRateCategory][sample] = stdDevSynthesisRate;
}


void Trace::updateStdDevSynthesisRateAcceptanceRatioTrace(double acceptanceLevel)
{
	stdDevSynthesisRateAcceptanceRatioTrace.push_back(acceptanceLevel);
}


void Trace::updateSynthesisRateAcceptanceRatioTrace(unsigned category, unsigned geneIndex, double acceptanceLevel)
{
	synthesisRateAcceptanceRatioTrace[category][geneIndex].push_back(acceptanceLevel);
}


void Trace::updateCodonSpecificAcceptanceRatioTrace(unsigned codonIndex, double acceptanceLevel)
{
	codonSpecificAcceptanceRatioTrace[codonIndex].push_back(acceptanceLevel);
}


void Trace::updateSynthesisRateTrace(unsigned sample, unsigned geneIndex,
	std::vector<std::vector <double>> &currentSynthesisRateLevel)
{
	for (unsigned category = 0; category < synthesisRateTrace.size(); category++)
	{
		synthesisRateTrace[category][geneIndex][sample] = currentSynthesisRateLevel[category][geneIndex];
	}
}


void Trace::updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex, unsigned value)
{
	mixtureAssignmentTrace[geneIndex][sample] = value;
}


void Trace::updateMixtureProbabilitiesTrace(unsigned samples, std::vector<double> &categoryProbabilities)
{
	for (unsigned category = 0; category < mixtureProbabilitiesTrace.size(); category++)
	{
		mixtureProbabilitiesTrace[category][samples] = categoryProbabilities[category];
	}
}


//----------------------------------//
//---------- ROC Specific ----------//
//----------------------------------//


void Trace::updateCodonSpecificParameterTraceForAA(unsigned sample, std::string aa,
	std::vector<std::vector<double>> &curParam, unsigned paramType)
{
	unsigned aaStart;
	unsigned aaEnd;
	SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);
	for (unsigned category = 0; category < codonSpecificParameterTrace[paramType].size(); category++)
	{
		for (unsigned i = aaStart; i < aaEnd; i++)
		{
			codonSpecificParameterTrace[paramType][category][i][sample] = curParam[category][i];
		}
	}
	/*
	switch (paramType) {
	case 0: 
		for (unsigned category = 0; category < codonSpecificParameterTraceOne.size(); category++)
		{
			for (unsigned i = aaStart; i < aaEnd; i++)
			{
				codonSpecificParameterTraceOne[category][i][sample] = curParam[category][i];
			}
		}
		break;
	case 1:
		for (unsigned category = 0; category < codonSpecificParameterTraceTwo.size(); category++)
		{
			for (unsigned i = aaStart; i < aaEnd; i++)
			{
				codonSpecificParameterTraceTwo[category][i][sample] = curParam[category][i];
			}
		}
		break;
	default:
		my_printError("ERROR: Unknown parameter type in updateCodonSpecificParameterTraceForAA\n");
		break;
	}
	*/
}


void Trace::updateSynthesisOffsetTrace(unsigned index, unsigned sample, double value)
{
	synthesisOffsetTrace[index][sample] = value;
}


void Trace::updateSynthesisOffsetAcceptanceRatioTrace(unsigned index, double value)
{
	synthesisOffsetAcceptanceRatioTrace[index].push_back(value);
}


void Trace::updateObservedSynthesisNoiseTrace(unsigned index, unsigned sample, double value)
{
	observedSynthesisNoiseTrace[index][sample] = value;
}


//-------------------------------------//
//---------- RFP Specific -------------//
//-------------------------------------//


void Trace::updateCodonSpecificParameterTraceForCodon(unsigned sample, std::string codon,
	std::vector<std::vector<double>> &curParam, unsigned paramType)
{
	unsigned i = SequenceSummary::codonToIndex(codon);
	for (unsigned category = 0; category < codonSpecificParameterTrace[paramType].size(); category++)
	{
		codonSpecificParameterTrace[paramType][category][i][sample] = curParam[category][i];
	}

	/*
	switch (paramType)
	{
		case 0:
			for (unsigned category = 0; category < codonSpecificParameterTraceOne.size(); category++)
			{
				codonSpecificParameterTraceOne[category][i][sample] = curParam[category][i];
			}
			break;
		case 1:
			for (unsigned category = 0; category < codonSpecificParameterTraceTwo.size(); category++)
			{
				codonSpecificParameterTraceTwo[category][i][sample] = curParam[category][i];
			}
			break;
		default:
			my_printError("ERROR: Unknown parameter type in updateCodonSpecificParameterTraceForCodon\n");
			break;
	}
	*/

}


// -----------------------------------------------------------------------------------------------------//
// ---------------------------------------- R SECTION --------------------------------------------------//
// -----------------------------------------------------------------------------------------------------//

#ifndef STANDALONE

//--------------------------------------//
//---------- Getter Functions ----------//
//--------------------------------------//
std::vector<double> Trace::getSynthesisRateAcceptanceRatioTraceByMixtureElementForGeneR(unsigned mixtureElement,
	unsigned geneIndex)
{
	std::vector<double> RV;
	bool checkGene = checkIndex(geneIndex, 1, synthesisRateAcceptanceRatioTrace.size());
	bool checkMixtureElement = checkIndex(mixtureElement, 1, mixtureProbabilitiesTrace.size());
	if (checkGene && checkMixtureElement)
	{
		RV = getSynthesisRateAcceptanceRatioTraceByMixtureElementForGene(mixtureElement - 1, geneIndex - 1);
	}
	return RV;
}


std::vector<double> Trace::getSynthesisRateTraceForGeneR(unsigned geneIndex)
{
	std::vector<double> RV;
	bool checkGene = checkIndex(geneIndex, 1, synthesisRateTrace[0].size());
	if (checkGene)
	{
		RV = getSynthesisRateTraceForGene(geneIndex - 1);
	}
	return RV;
}


std::vector<double> Trace::getSynthesisRateTraceByMixtureElementForGeneR(unsigned mixtureElement, unsigned geneIndex)
{
	std::vector<double> RV;
	bool checkMixtureElement = checkIndex(mixtureElement, 1, mixtureProbabilitiesTrace.size());
	bool checkGene = checkIndex(geneIndex, 1, synthesisRateTrace[0].size());
	if (checkMixtureElement && checkGene)
	{
		RV = getSynthesisRateTraceByMixtureElementForGene(mixtureElement - 1, geneIndex - 1);
	}
	return RV;
}


std::vector<unsigned> Trace::getMixtureAssignmentTraceForGeneR(unsigned geneIndex)
{
	std::vector <unsigned> RV;
	bool checkGene = checkIndex(geneIndex, 1, mixtureAssignmentTrace[0].size());
	if (checkGene)
	{
		RV = getMixtureAssignmentTraceForGene(geneIndex - 1);
	}
	return RV;
}


std::vector<double> Trace::getMixtureProbabilitiesTraceForMixtureR(unsigned mixtureIndex)
{
	std::vector<double> RV;
	bool check = checkIndex(mixtureIndex, 1, mixtureProbabilitiesTrace.size());
	if (check)
	{
		RV = getMixtureProbabilitiesTraceForMixture(mixtureIndex - 1);
	}
	return RV;
}


std::vector<std::vector<double>> Trace::getStdDevSynthesisRateTraces()
{
	return stdDevSynthesisRateTrace;
}


unsigned Trace::getNumberOfMixtures()
{
	return mixtureProbabilitiesTrace.size();
}

//--------------------------------------//
//---------- Setter Functions ----------//
//--------------------------------------//
void Trace::setStdDevSynthesisRateTraces(std::vector<std::vector<double>> _stdDevSynthesisRateTrace)
{
	stdDevSynthesisRateTrace = _stdDevSynthesisRateTrace;
}


void Trace::setStdDevSynthesisRateAcceptanceRatioTrace(std::vector<double> _stdDevSynthesisRateAcceptanceRatioTrace)
{
	stdDevSynthesisRateAcceptanceRatioTrace = _stdDevSynthesisRateAcceptanceRatioTrace;
}


void Trace::setSynthesisRateTrace(std::vector<std::vector<std::vector<double>>> _synthesisRateTrace)
{
	synthesisRateTrace = _synthesisRateTrace;
}


void Trace::setSynthesisRateAcceptanceRatioTrace(std::vector<std::vector<std::vector<double>>>_synthesisRateAcceptanceRatioTrace)
{
	synthesisRateAcceptanceRatioTrace = _synthesisRateAcceptanceRatioTrace;
}


void Trace::setMixtureAssignmentTrace(std::vector<std::vector<unsigned>> _mixtureAssignmentTrace)
{
	mixtureAssignmentTrace = _mixtureAssignmentTrace;
}


void Trace::setMixtureProbabilitiesTrace(std::vector<std::vector<double>> _mixtureProbabilitiesTrace)
{
	mixtureProbabilitiesTrace = _mixtureProbabilitiesTrace;
}


void Trace::setCodonSpecificAcceptanceRatioTrace(std::vector<std::vector<double>> _cspAcceptanceRatioTrace)
{
	codonSpecificAcceptanceRatioTrace = _cspAcceptanceRatioTrace;
}


void Trace::setCategories(std::vector<mixtureDefinition> &_categories)
{
	categories = &_categories;
}

//----------------------------------//
//---------- ROC Specific ----------//
//----------------------------------//
std::vector<double> Trace::getCodonSpecificParameterTraceByMixtureElementForCodonR(unsigned mixtureElement,
	std::string& codon, unsigned paramType, bool withoutReference)
{
	std::vector<double> RV;
	bool checkMixtureElement = checkIndex(mixtureElement, 1, getNumberOfMixtures());
	if (checkMixtureElement)
	{
		RV = getCodonSpecificParameterTraceByMixtureElementForCodon(mixtureElement - 1, codon, paramType, withoutReference);
	}
	return RV;
}


std::vector<std::vector<double>> Trace::getSynthesisOffsetTraceR()
{
	return synthesisOffsetTrace;
}


std::vector<std::vector<double>> Trace::getObservedSynthesisNoiseTraceR()
{
	return observedSynthesisNoiseTrace;
}


void Trace::setSynthesisOffsetTrace(std::vector<std::vector <double> > _NoiseOffsetTrace)
{
	synthesisOffsetTrace = _NoiseOffsetTrace;
}


void Trace::setSynthesisOffsetAcceptanceRatioTrace(std::vector<std::vector <double> > _NoiseOffsetAcceptanceRatioTrace)
{
	synthesisOffsetAcceptanceRatioTrace = _NoiseOffsetAcceptanceRatioTrace;
}


void Trace::setObservedSynthesisNoiseTrace(std::vector<std::vector <double> > _ObservedSynthesisNoiseTrace)
{
	observedSynthesisNoiseTrace = _ObservedSynthesisNoiseTrace;
}


void Trace::setCodonSpecificParameterTrace(std::vector<std::vector<std::vector<double>>> _parameterTrace, unsigned paramType)
{
	codonSpecificParameterTrace[paramType] = _parameterTrace;
	/*
	switch (paramType) {
	case 0:
		codonSpecificParameterTraceOne = _parameterTrace;
		break;
	case 1:
		codonSpecificParameterTraceTwo = _parameterTrace;
		break;
	default:
		my_printError("ERROR: Unknown parameter type in setCodonSpecificParameterTrace.\n");
		break;
	}
	*/
}


bool Trace::checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound)
{
	bool check = false;
	if (lowerbound <= index && index <= upperbound)
	{
		check = true;
	}
	else
	{
		my_printError("ERROR: Index: % is out of bounds. Index must be between % & %\n", index, lowerbound, upperbound);
	}

	return check;
}
#endif
