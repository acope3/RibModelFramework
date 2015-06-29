#include "RFPParameter.h"
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

//Constructors & Destructors:
RFPParameter::RFPParameter(std::string filename) : Parameter()
{
	initFromRestartFile(filename);
}


RFPParameter::RFPParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, std::vector<std::vector<unsigned>> thetaKMatrix, 
		bool splitSer = true, std::string _mutationSelectionState = "allUnique") : Parameter()
{
	initParameterSet(sphi, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
	initRFPParameterSet();
}


RFPParameter::~RFPParameter()
{
	//dtor 
	//TODO: Need to call Parameter's deconstructor
}


RFPParameter::RFPParameter(const RFPParameter& other) : Parameter(other)
{
	//TODO: impliment with the correct variables
}


RFPParameter::RFPParameter& operator=(const RFPParameter& rhs)
{
	if (this == &rhs) return *this; // handle self assignment

	Parameter::operator=(rhs);

	//TODO: finish implimentation
	return *this;
}


//Initialization functions:
void RFPParameter::initRFPParameterSet()
{
	//TODO: set up values
}


void RFPParameter::initAlpha(std::vector<double> alphaValues, unsigned mixtureElement, char aa)
{
	//TODO: set up values
}


void RFPParameter::initLambdaPrime(std::vector<double> lambdaPrimeValues, unsigned mixtureElement, char aa)
{
	//TODO: set up values
}

void ROCParameter::initMutationSelectionCategories(std::vector<std::string> files, unsigned numCategories, unsigned paramType)
{
	//TODO: Not sure we need this function here and in ROCParameter. This could possibly move up to Parameter. Some changes may
	//need to be made to the function in order for it to function. Mainly, the paramType needs to be better defined.
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
	numAcceptForAlphaAndLambdaPrime++;
	unsigned i = SequenceSummary::CodonToIndex(grouping);

	for(unsigned k = 0u; k < numMutationCategories; k++)
	{
		currentAlphaParameter[k][i] = proposedAlphaParameter[k][i];
	}
	for(unsigned k = 0u; k < numSelectionCategories; k++)
	{
		currentLambdaPrimeParameter[k][i] = proposedLambdaPrimeParameter[k][i];
	}
	currentiidSum = proposediidSum;
}


double RFPParameter::getPreviousCodonSpecificProposalWidth(unsigned index)
{
	return prev_std_csp[index];
}


double RFPParameter::getCurrentCodonSpecificProposalWidth(unsigned index)
{
	return curr_std_csp[index];
}


void RFPParameter::proposeCodonSpecificParameter()
{
	//TODO: Needs to be implimented
}


//functions to manage proposal widths:
void RFPParameter::adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth)
{
	//TODO: Needs to be implimented
}

//TODO: The only thing stopping this from moving up to Parameter is the trace
//stuff. Is there a way around this?
void RFPParameter::adaptSphiProposalWidth(unsigned adaptationWidth)
{
	double acceptanceLevel = (double)numAcceptForSphi / (double)adaptationWidth;
	traces.updateSphiAcceptanceRatioTrace(acceptanceLevel);
	if(acceptanceLevel < 0.2)
	{
		prev_std_sphi = std_sphi;
		std_sphi = std::max(0.01, std_sphi * 0.8);
	}
	if(acceptanceLevel > 0.3)
	{
		prev_std_sphi = std_sphi;
		std_sphi = std::min(100.0, std_sphi * 1.2);
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
				prev_std_phi[cat][i] = std_phi[cat][i];
				std_phi[cat][i] = std::max(0.01, std_phi[cat][i] * 0.8);
			}
			if(acceptanceLevel > 0.3)
			{
				prev_std_phi[cat][i] = std_phi[cat][i];
				std_phi[cat][i] = std::min(10.0, std_phi[cat][i] * 1.2);
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
	std::vector<double> lambdaPrimeParameterTrace = traces.getlambdaPrimeParameterTraceByMixtureElementForCodon(mixtureElement, codon);
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
void RFPParameter::getParameterForCategory(unsigned category, unsigned parameter, char aa, bool proposal, double* returnValue)
{
}


std::vector <double> RFPParameter::getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex)
{
}

