#include "include/ROC/ROCParameter.h"

#include <math.h>
#include <ctime>
#include <iostream>
#include <set>
#include <fstream>
#include <sstream>
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

ROCParameter::ROCParameter(std::string filename) : Parameter()
{
	initFromRestartFile(filename);
}

#ifndef STANDALONE
ROCParameter::ROCParameter(double sphi, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix, bool splitSer) : Parameter()
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
	initROCParameterSet();

}

ROCParameter::ROCParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, bool splitSer, std::string _mutationSelectionState) :
	Parameter()
{
	std::vector<std::vector<unsigned>> thetaKMatrix;
	initParameterSet(sphi, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
	initROCParameterSet();
}
#endif



ROCParameter::ROCParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, std::vector<std::vector<unsigned>> thetaKMatrix, 
		bool splitSer, std::string _mutationSelectionState) : Parameter()
{
	initParameterSet(sphi, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
	initROCParameterSet();
}


ROCParameter::~ROCParameter()
{
	//dtor
}

ROCParameter::ROCParameter(const ROCParameter& other) : Parameter(other)
{

	// proposal bias and std for codon specific parameter
	bias_csp = other.bias_csp;
	std_csp = other.std_csp;

	currentMutationParameter = other.currentMutationParameter;
	proposedMutationParameter = other.proposedMutationParameter;

	currentSelectionParameter = other.currentSelectionParameter;
	proposedSelectionParameter = other.proposedSelectionParameter;

	phiEpsilon = other.phiEpsilon;
	phiEpsilon_proposed = other.phiEpsilon_proposed;
}
ROCParameter& ROCParameter::operator=(const ROCParameter& rhs)
{
	if (this == &rhs) return *this; // handle self assignment

	Parameter::operator=(rhs);

	// proposal bias and std for codon specific parameter
	bias_csp = rhs.bias_csp;
	std_csp = rhs.std_csp;

	currentMutationParameter = rhs.currentMutationParameter;
	proposedMutationParameter = rhs.proposedMutationParameter;

	currentSelectionParameter = rhs.currentSelectionParameter;
	proposedSelectionParameter = rhs.proposedSelectionParameter;

	phiEpsilon = rhs.phiEpsilon;
	phiEpsilon_proposed = rhs.phiEpsilon_proposed;

	return *this;
}



void ROCParameter::initROCParameterSet()
{
	// proposal bias and std for codon specific parameter
	bias_csp = 0;
	std_csp.resize(numParam, 0.1); 
	numAcceptForMutationAndSelection.resize(22, 0u);

	phiEpsilon = 0.1;
	phiEpsilon_proposed = 0.1;

	//may need getter fcts
	currentMutationParameter.resize(numMutationCategories);
	proposedMutationParameter.resize(numMutationCategories);

	for (unsigned i = 0u; i < numMutationCategories; i++)
	{
		std::vector<double> tmp(numParam, 0.0);
		currentMutationParameter[i] = tmp;
		proposedMutationParameter[i] = tmp;
	}


	currentSelectionParameter.resize(numSelectionCategories);
	proposedSelectionParameter.resize(numSelectionCategories);


	for (unsigned i = 0u; i < numSelectionCategories; i++)
	{
		std::vector<double> tmp(numParam, 0.0);
		proposedSelectionParameter[i] = tmp;
		currentSelectionParameter[i] = tmp;
	}
}


std::vector<std::vector<double>> ROCParameter::calculateSelectionCoefficients(unsigned sample, unsigned mixture)
{
	unsigned numGenes = mixtureAssignment.size();
	std::vector<std::vector<double>> selectionCoefficients;
	selectionCoefficients.resize(numGenes);
	for (unsigned i = 0; i < numGenes; i++)
	{
		for (unsigned j = 0; j < 22; j++)
		{
			unsigned aaRange[2];
			SequenceSummary::AAindexToCodonRange(j, true, aaRange);
			std::vector <double> tmp;
			double minValue = 0.0;
			for (unsigned k = aaRange[0]; k < aaRange[1]; k++)
			{
				std::string codon = SequenceSummary::codonArrayParameter[k];
				tmp.push_back(getSelectionPosteriorMean(sample, mixture, codon));
				if (tmp[k] < minValue)
				{
					minValue = tmp[k];
				}
			}
			tmp.push_back(0.0);
			double phi = getSynthesisRatePosteriorMean(sample, i, mixture);
			for (unsigned k = 0; k < tmp.size(); k++)
			{
				tmp[k] -= minValue;
				selectionCoefficients[i].push_back(phi * tmp[k]);
			}
		}
	}
	return selectionCoefficients;
}


void ROCParameter::writeEntireRestartFile(std::string filename)
{
	writeBasicRestartFile(filename);
	writeROCRestartFile(filename);
}


void ROCParameter::writeROCRestartFile(std::string filename)
{
	std::ofstream out;
	out.open(filename.c_str(), std::ofstream::app);
	if (out.fail())
	{
		std::cerr <<"Could not open RestartFile.txt to append\n";
		std::exit(1);
	}

	std::ostringstream oss;
	unsigned j;
	oss <<">std_csp:\n";
	for (unsigned i = 0; i < std_csp.size(); i++)
	{
		oss << std_csp[i];
		if ((i + 1) % 10 == 0) oss <<"\n";
		else oss <<" ";
	}
	oss <<">currentMutationParameter:\n";
	for (unsigned i = 0; i < currentMutationParameter.size(); i++)
	{
		oss <<"***\n";
		for (j = 0; j < currentMutationParameter[i].size(); j++)
		{
			oss << currentMutationParameter[i][j];
			if ((j + 1) % 10 == 0) oss <<"\n";
			else oss <<" ";
		}
		if (j % 10 != 0) oss <<"\n";
	}

	oss <<">currentSelectionParameter:\n";
	for (unsigned i = 0; i < currentSelectionParameter.size(); i++)
	{
		oss <<"***\n";
		for (j = 0; j < currentSelectionParameter[i].size(); j++)
		{
			oss << currentSelectionParameter[i][j];
			if ((j + 1) % 10 == 0) oss <<"\n";
			else oss <<" ";
		}
		if (j % 10 != 0) oss <<"\n";
	}
	std::string output = oss.str();
	out << output;
	out.close();	
}


void ROCParameter::initFromRestartFile(std::string filename)
{
	initBaseValuesFromFile(filename);
	initROCValuesFromFile(filename);
}

void ROCParameter::initROCValuesFromFile(std::string filename)
{
	std::ifstream input;
	input.open(filename.c_str());
	if (input.fail())
	{
		std::cerr <<"Could not open RestartFile.txt to initialzie ROC values\n";
		std::exit(1);
	}
	std::string tmp, variableName;
	unsigned cat = 0;
	while (getline(input, tmp))
	{
		int flag;
		if (tmp[0] == '>') flag = 1;
		else if (input.eof() || tmp == "\n") flag = 2;
		else if (tmp[0] == '#') flag = 3;
		else flag = 4;

		if (flag == 1)
		{
			cat = 0;
			variableName = tmp.substr(1,tmp.size()-2);
		}
		else if (flag == 2)
		{
			std::cout <<"here\n";
		}
		else if (flag == 3) //user comment, continue
		{
			continue;
		}
		else 
		{
			std::istringstream iss;      
			if (variableName == "currentMutationParameter")
			{
				if (tmp == "***")
				{
					currentMutationParameter.resize(currentMutationParameter.size() + 1);
					cat++;
				}
				else if (tmp == "\n") continue;
				else
				{
					double val;
					iss.str(tmp);
					while (iss >> val)
					{
						currentMutationParameter[cat - 1].push_back(val);
					}
				}
			}
			else if (variableName == "currentSelectionParameter")
      {
        if (tmp == "***")
        {
          currentSelectionParameter.resize(currentSelectionParameter.size() + 1);
          cat++;
        }
				else if (tmp == "\n") continue;
        else
        {
          double val;
          iss.str(tmp);
          while (iss >> val)
          {
            currentSelectionParameter[cat - 1].push_back(val);
          }
        }
      }
			else if (variableName == "std_csp")
			{
				double val;
				iss.str(tmp);
				while (iss >> val)
				{
					std_csp.push_back(val);
				}
			}
		}
	}
	input.close();

	//init other values
	phiEpsilon = 0.1;
	phiEpsilon_proposed = 0.1;
	bias_csp = 0;
	numAcceptForMutationAndSelection.resize(22, 0u);
	proposedMutationParameter.resize(numMutationCategories);
	proposedSelectionParameter.resize(numSelectionCategories);
	for (unsigned i = 0; i < numMutationCategories; i++)
	{
		proposedMutationParameter[i] = currentMutationParameter[i];
	}
	for (unsigned i = 0; i < numSelectionCategories; i++)
	{
		proposedSelectionParameter[i] = currentSelectionParameter[i];
	}
}
#ifndef STANDALONE
	SEXP ROCParameter::calculateSelectionCoefficientsR(unsigned sample, unsigned mixture)
	{
		NumericMatrix RSelectionCoefficents(mixtureAssignment.size(), 62); //62 due to stop codons
		std::vector<std::vector<double>> selectionCoefficients;
		bool checkMixture = checkIndex(mixture, 1, numMixtures);
		if (checkMixture)
		{
			selectionCoefficients = calculateSelectionCoefficients(sample, mixture - 1);
			unsigned index = 0;
			for (unsigned i = 0; i < selectionCoefficients.size(); i++)
			{
				for(unsigned j = 0; j < selectionCoefficients[i].size(); j++, index++)
				{
					RSelectionCoefficents[index] = selectionCoefficients[i][j];
				}
			}
		}
		return RSelectionCoefficents;
	}
#endif

//This function seems to be used only for the purpose of R
	void ROCParameter::initSelection(std::vector<double> selectionValues, unsigned mixtureElement, char aa)
	{
		//TODO: seperate out the R wrapper functionality and make the wrapper
		//currentSelectionParameter
		bool check = checkIndex(mixtureElement, 1, numMixtures);
		if (check)
		{
			mixtureElement--;
			unsigned aaRange[2];
			int category = getSelectionCategory(mixtureElement);

			aa = std::toupper(aa);
			SequenceSummary::AAToCodonRange(aa, true, aaRange);
			for (unsigned i = aaRange[0], j = 0; i < aaRange[1]; i++, j++)
			{
				currentSelectionParameter[category][i] = selectionValues[j];
			}
		}
	}

//This function seems to be used only for the purpose of R
	void ROCParameter::initMutation(std::vector<double> mutationValues, unsigned mixtureElement, char aa)
	{
		//TODO: seperate out the R wrapper functionality and make the wrapper
		//currentMutationParameter
		bool check = checkIndex(mixtureElement, 1, numMixtures);
		if (check)
		{
			mixtureElement--;

			unsigned aaRange[2];
			unsigned category = getMutationCategory(mixtureElement);

			aa = std::toupper(aa);
			SequenceSummary::AAToCodonRange(aa, true, aaRange);
			for (unsigned i = aaRange[0], j = 0; i < aaRange[1]; i++, j++)
			{
				currentMutationParameter[category][i] = mutationValues[j];
			}
		}
	}


	void ROCParameter::initMutationSelectionCategories(std::vector<std::string> files, unsigned numCategories, unsigned paramType)
	{
		//should noise be a variable?
		unsigned i, j;
		std::size_t pos, pos2;
		std::ifstream currentFile;
		std::string tmpString;
		std::string type;
		if (paramType == ROCParameter::dM) type = "mu";
		else type = "eta";


		for (i = 0; i < numCategories; i++)
		{
			std::vector<double> temp(numParam, 0.0);

			//open the file, make sure it opens
			currentFile.open(files[i].c_str());
			if (currentFile.fail())
			{
				std::cerr <<"Error opening file " << i <<" in the file vector.\n";
				std::exit(1);
			}
			currentFile >> tmpString; //trash the first line, no info given.


			j = 0;
			while (currentFile >> tmpString)
			{
				pos = tmpString.find(",");
				pos2 = tmpString.find(",", pos + 1);
				if (pos != std::string::npos && pos2 != std::string::npos)
				{
					std::string val = tmpString.substr(pos + 1, pos2 - (pos + 1));
					if (tmpString.find(type) != std::string::npos) //mu or eta was found, depending on category
					{
						temp[j] = std::atof(val.c_str()); //std::stod(val);
						//					temp[j] = randNorm(temp[j], 0.3);
						j++;
						if (j == numParam) break;
					}
				}
			}
			unsigned altered = 0u;
			for (j = 0; j < categories.size(); j++)
			{
				if (paramType == ROCParameter::dM && categories[j].delM == i)
				{
					currentMutationParameter[j] = temp;
					proposedMutationParameter[j] = temp;
					altered++;
				}
				else if (paramType == ROCParameter::dEta && categories[j].delEta == i)
				{
					currentSelectionParameter[j] = temp;
					proposedSelectionParameter[j] = temp;
					altered++;
				}
				if (altered == numCategories) break; //to not access indicies out of bounds.
			}
			currentFile.close();
		}
	}


	void ROCParameter::getParameterForCategory(unsigned category, unsigned paramType, char aa, bool proposal, double* returnSet)
	{
		std::vector<double> *tempSet;
		if(paramType == ROCParameter::dM)
		{
			tempSet = (proposal ? &proposedMutationParameter[category] : &currentMutationParameter[category]);
		}
		else if(paramType == ROCParameter::dEta)
		{
			tempSet = (proposal ? &proposedSelectionParameter[category] : &currentSelectionParameter[category]);
		}
		else
		{
			std::cerr << "Warning in ROCParameter::getParameterForCategory: Unkown parameter type: " << paramType << "\n";
			std::cerr << "\tReturning mutation parameter! \n";
			tempSet = (proposal ? &proposedMutationParameter[category] : &currentMutationParameter[category]);
		}
		unsigned aaRange[2];
		SequenceSummary::AAToCodonRange(aa, true, aaRange);

		unsigned j = 0u;
		for(unsigned i = aaRange[0]; i < aaRange[1]; i++, j++)
		{
			returnSet[j] = tempSet -> at(i);
		}
	}


	double ROCParameter::getCurrentCodonSpecificProposalWidth(unsigned aa)
	{
		unsigned codonRange[2];
		SequenceSummary::AAindexToCodonRange(aa, true, codonRange);
		return std_csp[codonRange[0]];
	}

	void ROCParameter::adaptSphiProposalWidth(unsigned adaptationWidth)
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

	void ROCParameter::adaptSynthesisRateProposalWidth(unsigned adaptationWidth)
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


	double ROCParameter::getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement)
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


	double ROCParameter::getSphiPosteriorMean(unsigned samples)
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



	std::vector <double> ROCParameter::getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex)
	{
		std::vector<unsigned> mixtureAssignmentTrace = traces.getMixtureAssignmentTraceForGene(geneIndex);
		std::vector<double> probabilities(numMixtures, 0.0);
		unsigned traceLength = mixtureAssignmentTrace.size();

		if(samples > traceLength)
		{
			std::cerr << "Warning in Parameter::getMixtureAssignmentPosteriorMean throws: Number of anticipated samples (" <<
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


	double ROCParameter::getSphiVariance(unsigned samples, bool unbiased)
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

	double ROCParameter::getSynthesisRateVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement, bool unbiased)
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


	double ROCParameter::getMutationPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon)
	{
		double posteriorMean = 0.0;
		std::vector<double> mutationParameterTrace = traces.getMutationParameterTraceByMixtureElementForCodon(mixtureElement, codon);
		unsigned traceLength = mutationParameterTrace.size();

		if(samples > traceLength)
		{
			std::cerr << "Warning in ROCParameter::getMutationPosteriorMean throws: Number of anticipated samples (" <<
				samples << ") is greater than the length of the available trace (" << traceLength << ")." << "Whole trace is used for posterior estimate! \n";
			samples = traceLength;
		}
		unsigned start = traceLength - samples;
		for(unsigned i = start; i < traceLength; i++)
		{
			posteriorMean += mutationParameterTrace[i];
		}
		return posteriorMean / (double)samples;
	}


	double ROCParameter::getSelectionPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon)
	{
		double posteriorMean = 0.0;
		std::vector<double> selectionParameterTrace = traces.getSelectionParameterTraceByMixtureElementForCodon(mixtureElement, codon);
		unsigned traceLength = selectionParameterTrace.size();

		if(samples > traceLength)
		{
			std::cerr << "Warning in ROCParameter::getSelectionPosteriorMean throws: Number of anticipated samples (" <<
				samples << ") is greater than the length of the available trace (" << traceLength << ")." << "Whole trace is used for posterior estimate! \n";
			samples = traceLength;
		}
		unsigned start = traceLength - samples;
		for(unsigned i = start; i < traceLength; i++)
		{
			posteriorMean += selectionParameterTrace[i];
		}
		return posteriorMean / (double)samples;
	}


	double ROCParameter::getSelectionVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased)
	{
		std::vector<double> selectionParameterTrace = traces.getSelectionParameterTraceByMixtureElementForCodon(mixtureElement, codon);
		unsigned traceLength = selectionParameterTrace.size();
		if(samples > traceLength)
		{
			std::cerr << "Warning in ROCParameter::getSelectionVariance throws: Number of anticipated samples (" <<
				samples << ") is greater than the length of the available trace (" << traceLength << ")." << "Whole trace is used for posterior estimate! \n";
			samples = traceLength;
		}

		double posteriorMean = getSelectionPosteriorMean(mixtureElement, samples, codon);

		double posteriorVariance = 0.0;

		unsigned start = traceLength - samples;
		double difference = 0.0;
		for(unsigned i = start; i < traceLength; i++)
		{
			difference = selectionParameterTrace[i] - posteriorMean;
			posteriorVariance += difference * difference;
		}
		double normalizationTerm = unbiased ? (1/((double)samples-1.0)) : (1/(double)samples);
		return normalizationTerm * posteriorVariance;
	}


	double ROCParameter::getMutationVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased)
	{
		std::vector<double> mutationParameterTrace = traces.getMutationParameterTraceByMixtureElementForCodon(mixtureElement, codon);
		unsigned traceLength = mutationParameterTrace.size();
		if(samples > traceLength)
		{
			std::cerr << "Warning in ROCParameter::getMutationVariance throws: Number of anticipated samples (" <<
				samples << ") is greater than the length of the available trace (" << traceLength << ")." << "Whole trace is used for posterior estimate! \n";
			samples = traceLength;
		}

		double posteriorMean = getMutationPosteriorMean(mixtureElement, samples, codon);

		double posteriorVariance = 0.0;

		unsigned start = traceLength - samples;
		double difference = 0.0;
		for(unsigned i = start; i < traceLength; i++)
		{
			difference = mutationParameterTrace[i] - posteriorMean;
			posteriorVariance += difference * difference;
		}
		double normalizationTerm = unbiased ? (1/((double)samples-1.0)) : (1/(double)samples);
		return normalizationTerm * posteriorVariance;
	}

	std::string ROCParameter::getGrouping(unsigned index)
	{
		return groupList[index];
	}

	// Cedric: I decided to use a normal distribution to propose Sphi and phi instead of a lognormal because:
	// 1. It is a symmetric distribution and you therefore do not have to account for the unsymmetry in jump probabilities
	// 2. The one log and exp operation that have to be performed per parameter are cheaper than the operations necessary to draw from a lognormal
	// 3. phi has to be on a non log scale for the likelihood evaluation thus it does not help to keep phi on th elog scale all the time
	// 4. the adjusment of the likelihood by the jacobian that arises from this transformation is cheap and by grouping everything in one class it takes place more or less at the same place
	void ROCParameter::proposeCodonSpecificParameter()
	{
		/*
			 for(unsigned i = 0; i < numMutationCategories; i++)
			 {
			 proposedMutationParameter[i] = propose(currentMutationParameter[i], ROCParameter::randNorm, bias_csp, std_csp);
			 }
			 for(unsigned i = 0; i < numSelectionCategories; i++)
			 {
			 proposedSelectionParameter[i] = propose(currentSelectionParameter[i], ROCParameter::randNorm, bias_csp, std_csp);
			 }
		 */
		for (unsigned k = 0; k < 22; k++)
		{
			//need to skip W, M, X
			if (k == 21 || k == 18 || k == 10) continue;
			std::vector<double> iidProposed;
			unsigned aaRange[2];
			SequenceSummary::AAindexToCodonRange(k, true, aaRange);
			unsigned numCodons = aaRange[1] - aaRange[0];
			for(unsigned i = 0u; i < numCodons*(numMutationCategories + numSelectionCategories); i++)
			{
				iidProposed.push_back( randNorm(0.0, 1.0) );
			}

			std::vector<double> covaryingNums;
			covaryingNums = covarianceMatrix[k].transformIidNumersIntoCovaryingNumbers(iidProposed);
			for(unsigned i = 0; i < numMutationCategories; i++)
			{
				for(unsigned j = i * numCodons, l = aaRange[0]; j < (i * numCodons) + numCodons; j++, l++)
				{
					proposedMutationParameter[i][l] = currentMutationParameter[i][l] + covaryingNums[j];
				}
			}
			for(unsigned i = 0; i < numSelectionCategories; i++)
			{
				for(unsigned j = i * numCodons, l = aaRange[0]; j < (i * numCodons) + numCodons; j++, l++)
				{
					proposedSelectionParameter[i][l] = currentSelectionParameter[i][l] + covaryingNums[(numMutationCategories * numCodons) + j];
				}
			}
		}


	}


	std::vector<double> ROCParameter::propose(std::vector<double> currentParam, double (*proposal)(double a, double b), double A, std::vector<double> B)
	{
		unsigned numParam = currentParam.size();
		std::vector<double> proposedParam(numParam, 0.0);
		for(unsigned i = 0; i < numParam; i++)
		{
			proposedParam[i] = (*proposal)(A + currentParam[i], B[i]);
		}
		return proposedParam;
	}

	void ROCParameter::adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth)
	{
		unsigned numCSPsets = numAcceptForMutationAndSelection.size();
		std::cout << "acceptance ratio for amino acid:\n\t";
		for(unsigned i = 0; i < numCSPsets; i++)
		{
			if(i == 21 || i == 10 || i == 18) continue;
			std::cout << SequenceSummary::AminoAcidArray[i] << "\t";

			double acceptanceLevel = (double)numAcceptForMutationAndSelection[i] / (double)adaptationWidth;
			std::cout << acceptanceLevel << "\n";
			traces.updateCspAcceptanceRatioTrace(i, acceptanceLevel);
			unsigned codonRange[2];
			SequenceSummary::AAindexToCodonRange(i, true, codonRange);
			for(unsigned k = codonRange[0]; k < codonRange[1]; k++)
			{
				if(acceptanceLevel < 0.2)
				{
					covarianceMatrix[i] *= 0.8;
					covarianceMatrix[i].choleskiDecomposition();
					std_csp[k] *= 0.8;
				}
				if(acceptanceLevel > 0.3)
				{
					covarianceMatrix[i] *= 1.2;
					covarianceMatrix[i].choleskiDecomposition();
					std_csp[k] *= 1.2;
				}
			}
			numAcceptForMutationAndSelection[i] = 0u;
		}
		std::cout << "\n";
	}


	void ROCParameter::updateCodonSpecificParameter(std::string grouping)
	{
		char aa = grouping[0];
		unsigned aaRange[2];
		SequenceSummary::AAToCodonRange(aa, true, aaRange);
		unsigned aaIndex = SequenceSummary::aaToIndex.find(aa)->second;
		numAcceptForMutationAndSelection[aaIndex]++;

		for(unsigned k = 0u; k < numMutationCategories; k++)
		{
			for(unsigned i = aaRange[0]; i < aaRange[1]; i++)
			{
				currentMutationParameter[k][i] = proposedMutationParameter[k][i];
			}
		}
		for(unsigned k = 0u; k < numSelectionCategories; k++)
		{
			for(unsigned i = aaRange[0]; i < aaRange[1]; i++)
			{
				currentSelectionParameter[k][i] = proposedSelectionParameter[k][i];
			}
		}
	}



	void ROCParameter::initMutationSelectionCategoriesR(std::vector<std::string> files, unsigned numCategories, std::string paramType)
	{
		unsigned value;
		bool check = true;
		if (paramType == "Mutation")
		{
			value = ROCParameter::dM;
		}
		else if (paramType == "Selection")
		{
			value = ROCParameter::dEta;
		}
		else
		{
			std::cerr <<"Bad paramType given. Expected \"Mutation\" or \"Selection\".\nFunction not being executed!\n";
			check = false;
		}
		if (files.size() != numCategories) //we have different sizes and need to stop
		{
			std::cerr <<"The number of files given and the number of categories given differ. Function will not be executed!\n";
			check = false;
		}

		if (check)
		{
			initMutationSelectionCategories(files, numCategories, value);
		}
	}

