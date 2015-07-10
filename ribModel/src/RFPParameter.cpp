#include "include/RFP/RFPParameter.h"
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

//Constructors & Destructors:
RFPParameter::RFPParameter() : Parameter()
{
}


RFPParameter::RFPParameter(std::string filename) : Parameter(64)
{
	initFromRestartFile(filename);
	numParam = 61;
}


RFPParameter::RFPParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, std::vector<std::vector<unsigned>> thetaKMatrix, 
		bool splitSer, std::string _mutationSelectionState) : Parameter(64)
{
	initParameterSet(sphi, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
	initRFPParameterSet();
}


RFPParameter::~RFPParameter()
{
	//dtor 
	//TODO: Need to call Parameter's deconstructor
}

RFPParameter& RFPParameter::operator=(const RFPParameter& rhs)
{
	if (this == &rhs) return *this; // handle self assignment

	Parameter::operator=(rhs);

	currentAlphaParameter = rhs.currentAlphaParameter;
	proposedAlphaParameter = rhs.proposedAlphaParameter;
	currentLambdaPrimeParameter = rhs.proposedLambdaPrimeParameter;
	proposedLambdaPrimeParameter = rhs.proposedLambdaPrimeParameter;

	lambdaValues = rhs.lambdaValues;
	numAcceptForAlphaAndLambdaPrime = rhs.numAcceptForAlphaAndLambdaPrime;

	bias_csp = rhs.bias_csp;
	std_csp = rhs.std_csp;

	return *this;
}


//Initialization functions:
void RFPParameter::initRFPParameterSet()
{

	unsigned alphaCategories = getNumMutationCategories();
	unsigned lambdaPrimeCategories = getNumSelectionCategories();

	currentAlphaParameter.resize(alphaCategories);
	proposedAlphaParameter.resize(alphaCategories);
	currentLambdaPrimeParameter.resize(lambdaPrimeCategories);
	proposedLambdaPrimeParameter.resize(lambdaPrimeCategories);
	lambdaValues.resize(lambdaPrimeCategories);
	numParam = 61;

	for (unsigned i = 0; i < alphaCategories; i++)
	{
		std::vector <double> tmp(numParam,1.0);
		currentAlphaParameter[i] = tmp;
		proposedAlphaParameter[i] = tmp;
	}
	for (unsigned i = 0; i < lambdaPrimeCategories; i++)
	{
		std::vector <double> tmp(numParam,1.0);
		currentLambdaPrimeParameter[i] = tmp;
		proposedLambdaPrimeParameter[i] = tmp;
		lambdaValues[i] = tmp; //Maybe we don't initialize this one? or we do it differently?
	}

	numAcceptForAlphaAndLambdaPrime.resize(maxGrouping, 0u);
	bias_csp = 0;
	std_csp.resize(numParam, 0.1);

	groupList = {"GCA", "GCC", "GCG", "GCT", "TGC", "TGT", "GAC", "GAT", "GAA", "GAG",
		"TTC", "TTT", "GGA", "GGC", "GGG", "GGT", "CAC", "CAT", "ATA", "ATC",
		"ATT", "AAA", "AAG", "CTA", "CTC", "CTG", "CTT", "TTA", "TTG", "ATG",
		"AAC", "AAT", "CCA", "CCC", "CCG", "CCT", "CAA", "CAG", "AGA", "AGG",
		"CGA", "CGC", "CGG", "CGT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC",
		"ACG", "ACT", "GTA", "GTC", "GTG", "GTT", "TGG", "TAC", "TAT", "AGC",
		"AGT"};
}


void RFPParameter::initAlpha(double alphaValue, unsigned mixtureElement, std::string codon)
{
	unsigned category = getMutationCategory(mixtureElement);
	unsigned index = SequenceSummary::CodonToIndex(codon);
	currentAlphaParameter[category][index] = alphaValue;
}


void RFPParameter::initLambdaPrime(double lambdaPrimeValue, unsigned mixtureElement, std::string codon)
{
	unsigned category = getMutationCategory(mixtureElement);
	unsigned index = SequenceSummary::CodonToIndex(codon);
	currentLambdaPrimeParameter[category][index] = lambdaPrimeValue;
}


void RFPParameter::initMutationSelectionCategories(std::vector<std::string> files, unsigned numCategories, unsigned paramType)
{
	//TODO: Not sure we need this function here and in ROCParameter. This could possibly move up to Parameter. Some changes may
	//need to be made to the function in order for it to function. Mainly, the paramType needs to be better defined.

	//Not needed yet because we have no "true" values to give it.
}


//Restart file functions:
void RFPParameter::writeEntireRestartFile(std::string filename)
{
	//TODO: not a priority right now
}


void RFPParameter::writeRFPRestartFile(std::string filename)
{
	//TODO: not a priority right now
}


void RFPParameter::initFromRestartFile(std::string filename)
{
	//TODO: not a priority right now
}


void RFPParameter::initRFPValuesFromFile(std::string filename)
{
	//TODO: not a priority right now
}


//Codon Specific Parameter functions:
void RFPParameter::updateCodonSpecificParameter(std::string grouping)
{
	unsigned i = SequenceSummary::CodonToIndex(grouping);
	numAcceptForAlphaAndLambdaPrime[i]++;

	for(unsigned k = 0u; k < numMutationCategories; k++)
	{
		currentAlphaParameter[k][i] = proposedAlphaParameter[k][i];
	}
	for(unsigned k = 0u; k < numSelectionCategories; k++)
	{
		currentLambdaPrimeParameter[k][i] = proposedLambdaPrimeParameter[k][i];
	}
}


double RFPParameter::getCurrentCodonSpecificProposalWidth(unsigned index)
{
	return std_csp[index];
}


//TODO: Are we wanting to use a Covaraince Matrix structure?
void RFPParameter::proposeCodonSpecificParameter()
{
	unsigned numAlpha = currentAlphaParameter[0].size();
	unsigned numLambdaPrime = currentLambdaPrimeParameter[0].size();

	for (unsigned i = 0; i < numMutationCategories; i++)
	{
		for (unsigned j = 0; j < numAlpha; j++)
		{
			proposedAlphaParameter[i][j] = std::exp( randNorm( std::log(currentAlphaParameter[i][j]) , std_csp[j]) );
		}
	}

	for (unsigned i = 0; i < numSelectionCategories; i++)
	{
		for (unsigned j = 0; j < numLambdaPrime; j++)
		{
			proposedLambdaPrimeParameter[i][j] = std::exp( randNorm( std::log(currentLambdaPrimeParameter[i][j]) , std_csp[j]) );
		}
	}
}


//functions to manage proposal widths:
//TODO: Does not use a Covaraince Matrix if we need them.
void RFPParameter::adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth)
{

	std::cout << "acceptance ratio for codon:\n";
	for (unsigned i = 0; i < groupList.size(); i++)
	{
		std::cout << groupList[i] << "\t";

		unsigned codonIndex = SequenceSummary::CodonToIndex(groupList[i]);
		double acceptanceLevel = (double)numAcceptForAlphaAndLambdaPrime[codonIndex] / (double)adaptationWidth;
		std::cout << acceptanceLevel << "\n";
		traces.updateCspAcceptanceRatioTrace(codonIndex, acceptanceLevel);
		if (acceptanceLevel < 0.2)
		{
			std_csp[i] *= 0.8;
		}
		if (acceptanceLevel > 0.3)
		{
			std_csp[i] *= 1.2;
		}
		numAcceptForAlphaAndLambdaPrime[codonIndex] = 0u;
	}
	std::cout << "\n";
}


//TODO: The only thing stopping this from moving up to Parameter is the trace
//stuff. Is there a way around this?
void RFPParameter::adaptSphiProposalWidth(unsigned adaptationWidth)
{
	double acceptanceLevel = (double)numAcceptForSphi / (double)adaptationWidth;
	traces.updateSphiAcceptanceRatioTrace(acceptanceLevel);
	if(acceptanceLevel < 0.2)
	{
		std_sphi *= 0.8;
	}
	if(acceptanceLevel > 0.3)
	{
		std_sphi *= 1.2;
	}
	numAcceptForSphi = 0u;
}


//TODO: The only thing stopping this from moving up to Parameter is the trace
//stuff. Is there a way around this?
void RFPParameter::adaptSynthesisRateProposalWidth(unsigned adaptationWidth)
{
	for(unsigned cat = 0u; cat < numSelectionCategories; cat ++)
	{
		unsigned numGenes = numAcceptForSynthesisRate[cat].size();
		for(unsigned i = 0; i < numGenes; i++)
		{
			double acceptanceLevel = (double)numAcceptForSynthesisRate[cat][i] / (double)adaptationWidth;
			traces.updateSynthesisRateAcceptanceRatioTrace(cat, i, acceptanceLevel);
			if(acceptanceLevel < 0.2)
			{
				std_phi[cat][i] *= 0.8;
			}
			if(acceptanceLevel > 0.3)
			{
				std_phi[cat][i] *= 1.2;
			}
			numAcceptForSynthesisRate[cat][i] = 0u;
		}
	}
}


//Posterior Mean Functions:
//TODO: Traces prevent this from being in the parent class
double RFPParameter::getSphiPosteriorMean(unsigned samples)
{
	double posteriorMean = 0.0;
	std::vector<double> sPhiTrace = traces.getSPhiTrace();
	unsigned traceLength = sPhiTrace.size();

	if(samples > traceLength)
	{
		std::cerr << "Warning in Parameter::getSphiPosteriorMean throws: Number of anticipated samples (" <<
			samples << ") is greater than the length of the available trace (" << traceLength << ")." << "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}
	unsigned start = traceLength - samples;
	for(unsigned i = start; i < traceLength; i++)
	{
		posteriorMean += sPhiTrace[i];
	}
	return posteriorMean / (double)samples;
}


//TODO: Traces prevent this from being in the parent class
double RFPParameter::getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement)
{
	unsigned expressionCategory = getSynthesisRateCategory(mixtureElement);
	double posteriorMean = 0.0;
	std::vector<double> synthesisRateTrace = traces.getSynthesisRateTraceByMixtureElementForGene(mixtureElement, geneIndex);
	unsigned traceLength = synthesisRateTrace.size();


	if(samples > traceLength)
	{
		std::cerr << "Warning in Parameter::getSynthesisRatePosteriorMean throws: Number of anticipated samples (" <<
			samples << ") is greater than the length of the available trace (" << traceLength << ")." << "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}
	unsigned start = traceLength - samples;
	unsigned category = 0u;
	unsigned usedSamples = 0u;
	std::vector<unsigned> mixtureAssignmentTrace = traces.getMixtureAssignmentTraceForGene(geneIndex);
	for(unsigned i = start; i < traceLength; i++)
	{
		category = mixtureAssignmentTrace[i];
		category = getSynthesisRateCategory(category);
		if(category == expressionCategory)
		{
			posteriorMean += synthesisRateTrace[i];
			usedSamples++;
		}
	}
	// Can return NaN if gene was never in category! But that is Ok.
	return posteriorMean / (double)usedSamples;
}


double RFPParameter::getAlphaPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon)
{
	double posteriorMean = 0.0;
	std::vector<double> alphaParameterTrace = traces.getAlphaParameterTraceByMixtureElementForCodon(mixtureElement, codon);
	unsigned traceLength = alphaParameterTrace.size();

	if(samples > traceLength)
	{
		std::cerr << "Warning in RFPParameter::getAlphaPosteriorMean throws: Number of anticipated samples (" <<
			samples << ") is greater than the length of the available trace (" << traceLength << ")." << "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}
	unsigned start = traceLength - samples;
	for(unsigned i = start; i < traceLength; i++)
	{
		posteriorMean += alphaParameterTrace[i];
	}
	return posteriorMean / (double)samples;
}


double RFPParameter::getLambdaPrimePosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon)
{
	double posteriorMean = 0.0;
	std::vector<double> lambdaPrimeParameterTrace = traces.getLambdaPrimeParameterTraceByMixtureElementForCodon(mixtureElement, codon);
	unsigned traceLength = lambdaPrimeParameterTrace.size();

	if(samples > traceLength)
	{
		std::cerr << "Warning in RFPParameter::getLambdaPrimePosteriorMean throws: Number of anticipated samples (" <<
			samples << ") is greater than the length of the available trace (" << traceLength << ")." << "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}
	unsigned start = traceLength - samples;
	for(unsigned i = start; i < traceLength; i++)
	{
		posteriorMean += lambdaPrimeParameterTrace[i];
	}
	return posteriorMean / (double)samples;
}



//Variance Functions:
//TODO: Traces prevent this from being in the parent class
double RFPParameter::getSphiVariance(unsigned samples, bool unbiased)
{
	std::vector<double> sPhiTrace = traces.getSPhiTrace();
	unsigned traceLength = sPhiTrace.size();
	if(samples > traceLength)
	{
		std::cerr << "Warning in Parameter::getSphiVariance throws: Number of anticipated samples (" <<
			samples << ") is greater than the length of the available trace (" << traceLength << ")." << "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}
	double posteriorMean = getSphiPosteriorMean(samples);

	double posteriorVariance = 0.0;

	unsigned start = traceLength - samples;
	for(unsigned i = start; i < traceLength; i++)
	{
		double difference = sPhiTrace[i] - posteriorMean;
		posteriorVariance += difference * difference;
	}
	double normalizationTerm = unbiased ? (1/((double)samples-1.0)) : (1/(double)samples);
	return normalizationTerm * posteriorVariance;
}


//TODO: Traces prevent this from being in the parent class
double RFPParameter::getSynthesisRateVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement, bool unbiased)
{
	std::vector<double> synthesisRateTrace = traces.getSynthesisRateTraceByMixtureElementForGene(mixtureElement, geneIndex);
	unsigned traceLength = synthesisRateTrace.size();
	if(samples > traceLength)
	{
		std::cerr << "Warning in Parameter::getSynthesisRateVariance throws: Number of anticipated samples (" <<
			samples << ") is greater than the length of the available trace (" << traceLength << ")." << "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}

	double posteriorMean = getSynthesisRatePosteriorMean(samples, geneIndex, mixtureElement);

	double posteriorVariance = 0.0;
	if(!std::isnan(posteriorMean))
	{
		unsigned start = traceLength - samples;
		unsigned category = 0u;
		double difference = 0.0;
		std::vector<unsigned> mixtureAssignmentTrace = traces.getMixtureAssignmentTraceForGene(geneIndex);
		for(unsigned i = start; i < traceLength; i++)
		{
			category = mixtureAssignmentTrace[i];
			category = getSynthesisRateCategory(category);
			difference = synthesisRateTrace[i] - posteriorMean;
			posteriorVariance += difference * difference;
		}
	}
	double normalizationTerm = unbiased ? (1/((double)samples-1.0)) : (1/(double)samples);
	return normalizationTerm * posteriorVariance;
}


double RFPParameter::getAlphaVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased)
{
	std::vector<double> alphaParameterTrace = traces.getAlphaParameterTraceByMixtureElementForCodon(mixtureElement, codon);
	unsigned traceLength = alphaParameterTrace.size();
	if(samples > traceLength)
	{
		std::cerr << "Warning in RFPParameter::getAlphaVariance throws: Number of anticipated samples (" <<
			samples << ") is greater than the length of the available trace (" << traceLength << ")." << "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}

	double posteriorMean = getAlphaPosteriorMean(mixtureElement, samples, codon);

	double posteriorVariance = 0.0;

	unsigned start = traceLength - samples;
	double difference = 0.0;
	for(unsigned i = start; i < traceLength; i++)
	{
		difference = alphaParameterTrace[i] - posteriorMean;
		posteriorVariance += difference * difference;
	}
	double normalizationTerm = unbiased ? (1/((double)samples-1.0)) : (1/(double)samples);
	return normalizationTerm * posteriorVariance;
}


double RFPParameter::getLambdaPrimeVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased)
{
	std::vector<double> lambdaPrimeParameterTrace = traces.getLambdaPrimeParameterTraceByMixtureElementForCodon(mixtureElement, codon);
	unsigned traceLength = lambdaPrimeParameterTrace.size();
	if(samples > traceLength)
	{
		std::cerr << "Warning in RFPParameter::getSelectionVariance throws: Number of anticipated samples (" <<
			samples << ") is greater than the length of the available trace (" << traceLength << ")." << "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}

	double posteriorMean = getLambdaPrimePosteriorMean(mixtureElement, samples, codon);

	double posteriorVariance = 0.0;

	unsigned start = traceLength - samples;
	double difference = 0.0;
	for(unsigned i = start; i < traceLength; i++)
	{
		difference = lambdaPrimeParameterTrace[i] - posteriorMean;
		posteriorVariance += difference * difference;
	}
	double normalizationTerm = unbiased ? (1/((double)samples-1.0)) : (1/(double)samples);
	return normalizationTerm * posteriorVariance;
}


//Other functions:
double RFPParameter::getParameterForCategory(unsigned category, unsigned paramType, std::string codon, bool proposal)
{
	double rv;
	unsigned codonIndex = SequenceSummary::CodonToIndex(codon);
	if (paramType == RFPParameter::alp)
	{
		rv = (proposal ? proposedAlphaParameter[category][codonIndex] : currentAlphaParameter[category][codonIndex]);
	}
	else if (paramType == RFPParameter::lmPri)
	{
		rv = (proposal ? proposedLambdaPrimeParameter[category][codonIndex] : currentLambdaPrimeParameter[category][codonIndex]);
	}
	else
	{
		std::cerr << "Warning in RFPParameter::getParameterForCategory: Unkown parameter type: " << paramType << "\n";
		std::cerr << "\tReturning alpha parameter! \n";
		rv = (proposal ? proposedAlphaParameter[category][codonIndex] : currentAlphaParameter[category][codonIndex]);
	}

	return rv;
}


//TODO: Use of Trace prevents this from being in the base class
std::vector <double> RFPParameter::getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex)
{
	std::vector<unsigned> mixtureAssignmentTrace = traces.getMixtureAssignmentTraceForGene(geneIndex);
	std::vector<double> probabilities(numMixtures, 0.0);
	unsigned traceLength = mixtureAssignmentTrace.size();

	if (samples > traceLength)
	{
		std::cerr << "Warning in RFPParameter::getMixtureAssignmentPosteriorMean throws: Number of anticipated samples (" <<
			samples << ") is greater than the length of the available trace (" << traceLength << ")." << "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}

	unsigned start = traceLength - samples;
	for (unsigned i = start; i < traceLength; i++)
	{
		unsigned value = mixtureAssignmentTrace[i];
		probabilities[value]++;
	}

	for (unsigned i = 0; i < numMixtures; i++)
	{
		probabilities[i] /= (double)samples;
	}
	return probabilities;
}


//Statics:
const unsigned RFPParameter::alp = 0u;
const unsigned RFPParameter::lmPri = 1u;

//R Wrapper functions:
void RFPParameter::initAlphaR(double alphaValue, unsigned mixtureElement, std::string codon)
{
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	if (check)
	{
		mixtureElement--;
		codon[0] = std::toupper(codon[0]);
		codon[1] = std::toupper(codon[1]);
		codon[2] = std::toupper(codon[2]);

		initAlpha(alphaValue, mixtureElement, codon);
	}
}


void RFPParameter::initLambdaPrimeR(double lambdaPrimeValue, unsigned mixtureElement, std::string codon)
{
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	if (check)
	{
		mixtureElement--;
		codon[0] = std::toupper(codon[0]);
		codon[1] = std::toupper(codon[1]);
		codon[2] = std::toupper(codon[2]);

		initLambdaPrime(lambdaPrimeValue, mixtureElement, codon);
	}
}


double RFPParameter::getParameterForCategoryR(unsigned mixtureElement, unsigned paramType, std::string codon, bool proposal)
{
	double rv = 0.0;
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	if (check)
	{
		mixtureElement--;
		unsigned category = 0;
		codon[0] = std::toupper(codon[0]);
		codon[1] = std::toupper(codon[1]);
		codon[2] = std::toupper(codon[2]);
		if (paramType == RFPParameter::alp)
		{
			//THIS NEEDS TO CHANGE!!!!
			category = getMutationCategory(mixtureElement); //really alpha here
		}
		else if (paramType == RFPParameter::lmPri)
		{
			category = getSelectionCategory(mixtureElement);
		}
		rv = getParameterForCategory(category, paramType, codon, proposal);
	}
	return rv;
}


#ifndef STANDALONE
RFPParameter::RFPParameter(double sphi, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix, bool splitSer) : Parameter(64)
{
  unsigned _numMixtures = _matrix.size() / 2;
  std::vector<std::vector<unsigned>> thetaKMatrix;
  thetaKMatrix.resize(_numMixtures);

  unsigned index = 0;
  for (unsigned i = 0; i < _numMixtures; i++)
  {
    for (unsigned j = 0; j < 2; j++, index++)
    {
      thetaKMatrix[i].push_back(_matrix[index]);
    }
  }
  initParameterSet(sphi, _matrix.size() / 2, geneAssignment, thetaKMatrix, splitSer);
  initRFPParameterSet();

}

RFPParameter::RFPParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, bool splitSer, std::string _mutationSelectionState) :
Parameter(64)
{
  std::vector<std::vector<unsigned>> thetaKMatrix;
  initParameterSet(sphi, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
  initRFPParameterSet();
}
#endif
