// FONSEParameter.cpp

#include "include/FONSE/FONSEParameter.h"

#include <cmath>
#include <ctime>
#include <iostream>
#include <set>
#include <fstream>
#include <sstream>

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

const unsigned FONSEParameter::dM = 0;
const unsigned FONSEParameter::dOmega = 1;

FONSEParameter::FONSEParameter() : Parameter()
{
}

FONSEParameter::FONSEParameter(std::string filename) : Parameter()
{
	initFromRestartFile(filename);
}

#ifndef STANDALONE
FONSEParameter::FONSEParameter(double sphi, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix, bool splitSer) 
	: Parameter()
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
	initFONSEParameterSet();

}

FONSEParameter::FONSEParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, bool splitSer, std::string _mutationSelectionState)
	: Parameter()
{
	std::vector<std::vector<unsigned>> thetaKMatrix;
	initParameterSet(sphi, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
	initFONSEParameterSet();
}
#endif

FONSEParameter::FONSEParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
	std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer, std::string _mutationSelectionState) :
	Parameter()
{
	initParameterSet(sphi, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
	initFONSEParameterSet();
}

FONSEParameter::~FONSEParameter()
{
	// destructor
}

FONSEParameter::FONSEParameter(const FONSEParameter &other) : Parameter(other)
{
	bias_csp = other.bias_csp;
	std_csp = other.std_csp;

	currentMutationParameter = other.currentMutationParameter;
	proposedMutationParameter = other.proposedMutationParameter;

	currentSelectionParameter = other.currentSelectionParameter;
	proposedSelectionParameter = other.proposedSelectionParameter;

	phiEpsilon = other.phiEpsilon;
	phiEpsilon_proposed = other.phiEpsilon_proposed;

	numAcceptForMutationAndSelection = other.numAcceptForMutationAndSelection;

}

FONSEParameter& FONSEParameter::operator=(const FONSEParameter& rhs)
{
	if (this == &rhs)
		return *this; // handle self assignment

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

	covarianceMatrix = rhs.covarianceMatrix;

	numAcceptForMutationAndSelection = rhs.numAcceptForMutationAndSelection;

	return *this;
}

void FONSEParameter::initFONSEParameterSet()
{
	groupList = { "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R", "S", "T", "V", "Y", "Z" };
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

	CodonTable *codonTable = CodonTable::getInstance();
	for (unsigned i = 0; i < maxGrouping; i++)
	{
		std::string aa = CodonTable::aminoAcidArray[i];
		unsigned numCodons = codonTable -> getNumCodonsForAA(aa, true);
		CovarianceMatrix m((numMutationCategories + numSelectionCategories) * numCodons);
		m.choleskiDecomposition();
		covarianceMatrix.push_back(m);
	}

}

std::vector <std::vector <double> > FONSEParameter::calculateSelectionCoefficients(unsigned sample, unsigned mixture)
{
	unsigned numGenes = (unsigned)mixtureAssignment.size();
	std::vector<std::vector<double>> selectionCoefficients;
	selectionCoefficients.resize(numGenes);
	for (unsigned i = 0; i < numGenes; i++)
	{
		for (unsigned j = 0; j < getGroupListSize(); j++)
		{

			std::string aa = getGrouping(j);
			CodonTable *codonTable = CodonTable::getInstance();
			std::vector <unsigned> codonRange = codonTable -> AAToCodonRange(aa, true); //checked
			std::vector<double> tmp;
			double minValue = 0.0;
			for (unsigned k = 0; k < codonRange.size(); k++)
			{
				std::string codon = codonTable -> getForParamVectorCodon(codonRange[k]);
				tmp.push_back(getSelectionPosteriorMean(sample, mixture, codon));
				if (tmp[codonRange[k]] < minValue)
				{
					minValue = tmp[codonRange[k]];
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

void FONSEParameter::writeEntireRestartFile(std::string filename)
{
	writeBasicRestartFile(filename);
	writeFONSERestartFile(filename);
}

void FONSEParameter::writeFONSERestartFile(std::string filename)
{
	std::ofstream out;
	out.open(filename.c_str(), std::ofstream::app);
	if (out.fail())
	{
		std::cerr << "Could not open RestartFile.txt to append\n";
		std::exit(1);
	}

	std::ostringstream oss;
	unsigned j;
	oss << ">std_csp:\n";
	for (unsigned i = 0; i < std_csp.size(); i++)
	{
		oss << std_csp[i];
		if ((i + 1) % 10 == 0)
			oss << "\n";
		else
			oss << " ";
	}
	oss << ">currentMutationParameter:\n";
	for (unsigned i = 0; i < currentMutationParameter.size(); i++)
	{
		oss << "***\n";
		for (j = 0; j < currentMutationParameter[i].size(); j++)
		{
			oss << currentMutationParameter[i][j];
			if ((j + 1) % 10 == 0)
				oss << "\n";
			else
				oss << " ";
		}
		if (j % 10 != 0)
			oss << "\n";
	}

	oss << ">currentSelectionParameter:\n";
	for (unsigned i = 0; i < currentSelectionParameter.size(); i++)
	{
		oss << "***\n";
		for (j = 0; j < currentSelectionParameter[i].size(); j++)
		{
			oss << currentSelectionParameter[i][j];
			if ((j + 1) % 10 == 0)
				oss << "\n";
			else
				oss << " ";
		}
		if (j % 10 != 0)
			oss << "\n";
	}
	std::string output = oss.str();
	out << output;
	out.close();
}

void FONSEParameter::initFromRestartFile(std::string filename)
{
	initBaseValuesFromFile(filename);
	initFONSEValuesFromFile(filename);
}

void FONSEParameter::initFONSEValuesFromFile(std::string filename)
{
	std::ifstream input;
	input.open(filename.c_str());
	if (input.fail())
	{
		std::cerr << "Could not open RestartFile.txt to initialize FONSE values\n";
		std::exit(1);
	}
	std::string tmp, variableName;
	unsigned cat = 0;
	while (getline(input, tmp))
	{
		int flag;
		if (tmp[0] == '>')
			flag = 1;
		else if (input.eof() || tmp == "\n")
			flag = 2;
		else if (tmp[0] == '#')
			flag = 3;
		else
			flag = 4;

		if (flag == 1)
		{
			cat = 0;
			variableName = tmp.substr(1, tmp.size() - 2);
		}
		else if (flag == 2)
		{
			std::cout << "here\n";
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
				else if (tmp == "\n")
					continue;
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
				else if (tmp == "\n")
					continue;
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

	groupList = { "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R", "S", "T", "V", "Y", "Z" };
	//groupList = { "C", "D", "E", "F", "H", "K", "M", "N", "Q", "W", "Y" };
}

#ifndef STANDALONE
SEXP FONSEParameter::calculateSelectionCoefficientsR(unsigned sample, unsigned mixture)
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
			for (unsigned j = 0; j < selectionCoefficients[i].size(); j++, index++)
			{
				RSelectionCoefficents[index] = selectionCoefficients[i][j];
			}
		}
	}
	return RSelectionCoefficents;
}

void FONSEParameter::initCovarianceMatrix(SEXP _matrix, std::string aa)
{
	std::vector<double> tmp;
	NumericMatrix matrix(_matrix);

	for (unsigned i = 0u; i < aa.length(); i++)	aa[i] = (char)std::toupper(aa[i]);

	CodonTable *codonTable = CodonTable::getInstance();
	unsigned aaIndex = codonTable -> AAToAAIndex(aa);
	unsigned numRows = matrix.nrow();
	std::vector<double> covMatrix(numRows * numRows);

	//NumericMatrix stores the matrix by column, not by row. The loop
	//below transposes the matrix when it stores it.
	unsigned index = 0;
	for (unsigned i = 0; i < numRows; i++)
	{
		for (unsigned j = i; j < numRows * numRows; j += numRows, index++)
		{
			covMatrix[index] = matrix[j];
		}
	}
	CovarianceMatrix m(covMatrix);
	m.choleskiDecomposition();
	covarianceMatrix[aaIndex] = m;
}
#endif

CovarianceMatrix& FONSEParameter::getCovarianceMatrixForAA(std::string aa)
{
	CodonTable *codonTable = CodonTable::getInstance();
	aa[0] = (char)std::toupper(aa[0]);
	unsigned aaIndex = codonTable -> AAToAAIndex(aa);
	return covarianceMatrix[aaIndex];
}

void FONSEParameter::initSelection(std::vector<double> selectionValues, unsigned mixtureElement, std::string aa)
{
	CodonTable *codonTable = CodonTable::getInstance();
	//TODO: seperate out the R wrapper functionality and make the wrapper
	//currentSelectionParameter
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	if (check)
	{
		mixtureElement--;

		int category = getSelectionCategory(mixtureElement);

		aa[0] = (char)std::toupper(aa[0]);
		std::vector <unsigned> codonRange = codonTable -> AAToCodonRange(aa, true); //checked
		for (unsigned i = 0, j = 0; i < codonRange.size(); i++, j++)
		{
			currentSelectionParameter[category][codonRange[i]] = selectionValues[j];
		}
	}
}

void FONSEParameter::initMutation(std::vector<double> mutationValues, unsigned mixtureElement, std::string aa)
{
	CodonTable *codonTable = CodonTable::getInstance();
	//TODO: seperate out the R wrapper functionality and make the wrapper
	//currentMutationParameter
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	if (check)
	{
		mixtureElement--;
		unsigned category = getMutationCategory(mixtureElement);
		aa[0] = (char)std::toupper(aa[0]);
		std::vector <unsigned> codonRange = codonTable -> AAToCodonRange(aa, true); //checked
		for (unsigned i = 0, j = 0; i < 8; i++, j++)
		{
			currentMutationParameter[category][codonRange[i]] = mutationValues[j];
		}
	}
}

void FONSEParameter::initMutationSelectionCategories(std::vector<std::string> files, unsigned numCategories,
	unsigned paramType)
{
	//should noise be a variable?
	unsigned i, j;
	std::size_t pos, pos2;
	std::ifstream currentFile;
	std::string tmpString;
	std::string type;
	if (paramType == FONSEParameter::dM)
		type = "mu";
	else
		type = "eta";

	for (i = 0; i < numCategories; i++)
	{
		std::vector<double> temp(numParam, 0.0);

		//open the file, make sure it opens
		currentFile.open(files[i].c_str());
		if (currentFile.fail())
		{
			std::cerr << "Error opening file " << i << " in the file vector.\n";
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
					if (j == numParam)
						break;
				}
			}
		}
		unsigned altered = 0u;
		for (j = 0; j < categories.size(); j++)
		{
			if (paramType == FONSEParameter::dM && categories[j].delM == i)
			{
				currentMutationParameter[j] = temp;
				proposedMutationParameter[j] = temp;
				altered++;
			}
			else if (paramType == FONSEParameter::dOmega && categories[j].delEta == i)
			{
				currentSelectionParameter[j] = temp;
				proposedSelectionParameter[j] = temp;
				altered++;
			}
			if (altered == numCategories)
				break; //to not access indicies out of bounds.
		}
		currentFile.close();
	}
}

void FONSEParameter::getParameterForCategory(unsigned category, unsigned paramType, std::string aa, bool proposal,
	double *returnSet)
{
	std::vector<double> *tempSet;
	if (paramType == FONSEParameter::dM)
	{
		tempSet = (proposal ? &proposedMutationParameter[category] : &currentMutationParameter[category]);
	}
	else if (paramType == FONSEParameter::dOmega)
	{
		tempSet = (proposal ? &proposedSelectionParameter[category] : &currentSelectionParameter[category]);
	}
	else
	{
		std::cerr << "Warning in FONSEParameter::getParameterForCategory: Unknown parameter type: " << paramType << "\n";
		std::cerr << "\tReturning mutation parameter! \n";
		tempSet = (proposal ? &proposedMutationParameter[category] : &currentMutationParameter[category]);
	}
	CodonTable *codonTable = CodonTable::getInstance();
	std::vector <unsigned> codonRange = codonTable -> AAToCodonRange(aa, true); //checked

	unsigned j = 0u;
	for (unsigned i = 0; i < codonRange.size(); i++, j++)
	{
		returnSet[j] = tempSet->at(codonRange[i]);
	}
}

double FONSEParameter::getCurrentCodonSpecificProposalWidth(unsigned aa)
{
	CodonTable *codonTable = CodonTable::getInstance();
	std::vector <unsigned> codonRange = codonTable -> AAIndexToCodonRange(aa, true); //checked
	return std_csp[codonRange[0]];
}

void FONSEParameter::setNumPhiGroupings(unsigned _phiGroupings)
{
	phiGroupings = _phiGroupings;
}

void FONSEParameter::adaptSphiProposalWidth(unsigned adaptationWidth)
{
	double acceptanceLevel = (double)numAcceptForSphi / (double)adaptationWidth;
	traces.updateSphiAcceptanceRatioTrace(acceptanceLevel);
	if (acceptanceLevel < 0.2)
	{
		std_sphi *= 0.8;
	}
	if (acceptanceLevel > 0.3)
	{
		std_sphi *= 1.2;
	}
	numAcceptForSphi = 0u;
}

void FONSEParameter::adaptSynthesisRateProposalWidth(unsigned adaptationWidth)
{
	unsigned acceptanceUnder = 0u;
	unsigned acceptanceOver = 0u;
	for (unsigned cat = 0u; cat < numSelectionCategories; cat++)
	{
		unsigned numGenes = numAcceptForSynthesisRate[cat].size();
		for (unsigned i = 0; i < numGenes; i++)
		{
			double acceptanceLevel = (double)numAcceptForSynthesisRate[cat][i] / (double)adaptationWidth;
			traces.updateSynthesisRateAcceptanceRatioTrace(cat, i, acceptanceLevel);
			if (acceptanceLevel < 0.225)
			{
				std_phi[cat][i] *= 0.8;
				if (acceptanceLevel < 0.2) acceptanceUnder++;
			}
			if (acceptanceLevel > 0.275)
			{
				std_phi[cat][i] *= 1.2;
				if (acceptanceLevel > 0.3) acceptanceOver++;
			}
			numAcceptForSynthesisRate[cat][i] = 0u;
		}
	}
	std::cout << "acceptance ratio for synthesis rate:\n";
	std::cout << "\t acceptance ratio to low: " << acceptanceUnder << "\n";
	std::cout << "\t acceptance ratio to high: " << acceptanceOver << "\n";
}

void FONSEParameter::adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth)
{
	CodonTable *codonTable = CodonTable::getInstance();
	unsigned numCSPsets = numAcceptForMutationAndSelection.size();
	std::cout << "acceptance ratio for amino acid:\n\t";
	for (unsigned i = 0; i < numCSPsets; i++)
	{
		if (i == 21 || i == 10 || i == 18)
			continue;
		std::cout << CodonTable::aminoAcidArray[i] << "\t\t";
	}
	std::cout << "\n\t";
	for (unsigned i = 0; i < numCSPsets; i++)
	{
		if (i == 21 || i == 10 || i == 18)
			continue;
		double acceptanceLevel = (double)numAcceptForMutationAndSelection[i] / (double)adaptationWidth;
		std::cout << acceptanceLevel << "\t";
		traces.updateCspAcceptanceRatioTrace(i, acceptanceLevel);
		std::vector <unsigned> codonRange = codonTable -> AAIndexToCodonRange(i, true); //checked
		for (unsigned k = 0; k < 8; k++)
		{
			if (codonRange[k] == 100) break;
			if (acceptanceLevel < 0.2)
			{
				covarianceMatrix[i] *= 0.8;
				covarianceMatrix[i].choleskiDecomposition();
				std_csp[codonRange[k]] *= 0.8;
			}
			if (acceptanceLevel > 0.3)
			{
				covarianceMatrix[i] *= 1.2;
				covarianceMatrix[i].choleskiDecomposition();
				std_csp[codonRange[k]] *= 1.2;
			}
		}
		numAcceptForMutationAndSelection[i] = 0u;
	}
	std::cout << "\n";
}

double FONSEParameter::getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement)
{
	unsigned expressionCategory = getSynthesisRateCategory(mixtureElement);
	double posteriorMean = 0.0;
	std::vector<double> synthesisRateTrace = traces.getSynthesisRateTraceByMixtureElementForGene(mixtureElement,
		geneIndex);
	unsigned traceLength = (unsigned)synthesisRateTrace.size();

	if (samples > traceLength)
	{
		std::cerr << "Warning in Parameter::getSynthesisRatePosteriorMean throws: Number of anticipated samples ("
			<< samples << ") is greater than the length of the available trace (" << traceLength << ")."
			<< "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}
	unsigned start = traceLength - samples;
	unsigned category = 0u;
	unsigned usedSamples = 0u;
	std::vector<unsigned> mixtureAssignmentTrace = traces.getMixtureAssignmentTraceForGene(geneIndex);
	for (unsigned i = start; i < traceLength; i++)
	{
		category = mixtureAssignmentTrace[i];
		category = getSynthesisRateCategory(category);
		if (category == expressionCategory)
		{
			posteriorMean += synthesisRateTrace[i];
			usedSamples++;
		}
	}
	// Can return NaN if gene was never in category! But that is Ok.
	return posteriorMean / (double)usedSamples;
}

double FONSEParameter::getSphiPosteriorMean(unsigned samples)
{
	double posteriorMean = 0.0;
	std::vector<double> sPhiTrace = traces.getSphiTrace();
	unsigned traceLength = sPhiTrace.size();

	if (samples > traceLength)
	{
		std::cerr << "Warning in Parameter::getSphiPosteriorMean throws: Number of anticipated samples (" << samples
			<< ") is greater than the length of the available trace (" << traceLength << ")."
			<< "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}
	unsigned start = traceLength - samples;
	for (unsigned i = start; i < traceLength; i++)
	{
		posteriorMean += sPhiTrace[i];
	}
	return posteriorMean / (double)samples;
}

std::vector<double> FONSEParameter::getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex)
{
	std::vector<unsigned> mixtureAssignmentTrace = traces.getMixtureAssignmentTraceForGene(geneIndex);
	std::vector<double> probabilities(numMixtures, 0.0);
	unsigned traceLength = mixtureAssignmentTrace.size();

	if (samples > traceLength)
	{
		std::cerr
			<< "Warning in FONSEParameter::getMixtureAssignmentPosteriorMean throws: Number of anticipated samples ("
			<< samples << ") is greater than the length of the available trace (" << traceLength << ")."
			<< "Whole trace is used for posterior estimate! \n";
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

double FONSEParameter::getSphiVariance(unsigned samples, bool unbiased)
{
	std::vector<double> sPhiTrace = traces.getSphiTrace();
	unsigned traceLength = sPhiTrace.size();
	if (samples > traceLength)
	{
		std::cerr << "Warning in Parameter::getSphiVariance throws: Number of anticipated samples (" << samples
			<< ") is greater than the length of the available trace (" << traceLength << ")."
			<< "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}
	double posteriorMean = getSphiPosteriorMean(samples);

	double posteriorVariance = 0.0;

	unsigned start = traceLength - samples;
	for (unsigned i = start; i < traceLength; i++)
	{
		double difference = sPhiTrace[i] - posteriorMean;
		posteriorVariance += difference * difference;
	}
	double normalizationTerm = unbiased ? (1 / ((double)samples - 1.0)) : (1 / (double)samples);
	return normalizationTerm * posteriorVariance;
}

double FONSEParameter::getSynthesisRateVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement,
	bool unbiased)
{
	std::vector<double> synthesisRateTrace = traces.getSynthesisRateTraceByMixtureElementForGene(mixtureElement,
		geneIndex);
	unsigned traceLength = synthesisRateTrace.size();
	if (samples > traceLength)
	{
		std::cerr << "Warning in Parameter::getSynthesisRateVariance throws: Number of anticipated samples (" << samples
			<< ") is greater than the length of the available trace (" << traceLength << ")."
			<< "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}

	double posteriorMean = getSynthesisRatePosteriorMean(samples, geneIndex, mixtureElement);

	double posteriorVariance = 0.0;
	if (!std::isnan(posteriorMean))
	{
		unsigned start = traceLength - samples;
		unsigned category = 0u;
		double difference = 0.0;
		std::vector<unsigned> mixtureAssignmentTrace = traces.getMixtureAssignmentTraceForGene(geneIndex);
		for (unsigned i = start; i < traceLength; i++)
		{
			category = mixtureAssignmentTrace[i];
			category = getSynthesisRateCategory(category);
			difference = synthesisRateTrace[i] - posteriorMean;
			posteriorVariance += difference * difference;
		}
	}
	double normalizationTerm = unbiased ? (1 / ((double)samples - 1.0)) : (1 / (double)samples);
	return normalizationTerm * posteriorVariance;
}

double FONSEParameter::getMutationPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon)
{
	double posteriorMean = 0.0;
	std::vector<double> mutationParameterTrace = traces.getMutationParameterTraceByMixtureElementForCodon(
		mixtureElement, codon);
	unsigned traceLength = mutationParameterTrace.size();

	if (samples > traceLength)
	{
		std::cerr << "Warning in FONSEParameter::getMutationPosteriorMean throws: Number of anticipated samples ("
			<< samples << ") is greater than the length of the available trace (" << traceLength << ")."
			<< "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}
	unsigned start = traceLength - samples;
	for (unsigned i = start; i < traceLength; i++)
	{
		posteriorMean += mutationParameterTrace[i];
	}
	return posteriorMean / (double)samples;
}

double FONSEParameter::getSelectionPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon)
{
	double posteriorMean = 0.0;
	std::vector<double> selectionParameterTrace = traces.getSelectionParameterTraceByMixtureElementForCodon(
		mixtureElement, codon);
	unsigned traceLength = selectionParameterTrace.size();

	if (samples > traceLength)
	{
		std::cerr << "Warning in FONSEParameter::getSelectionPosteriorMean throws: Number of anticipated samples ("
			<< samples << ") is greater than the length of the available trace (" << traceLength << ")."
			<< "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}
	unsigned start = traceLength - samples;
	for (unsigned i = start; i < traceLength; i++)
	{
		posteriorMean += selectionParameterTrace[i];
	}
	return posteriorMean / (double)samples;
}

double FONSEParameter::getSelectionVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased)
{
	std::vector<double> selectionParameterTrace = traces.getSelectionParameterTraceByMixtureElementForCodon(
		mixtureElement, codon);
	unsigned traceLength = selectionParameterTrace.size();
	if (samples > traceLength)
	{
		std::cerr << "Warning in FONSEParameter::getSelectionVariance throws: Number of anticipated samples (" << samples
			<< ") is greater than the length of the available trace (" << traceLength << ")."
			<< "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}

	double posteriorMean = getSelectionPosteriorMean(mixtureElement, samples, codon);

	double posteriorVariance = 0.0;

	unsigned start = traceLength - samples;
	double difference = 0.0;
	for (unsigned i = start; i < traceLength; i++)
	{
		difference = selectionParameterTrace[i] - posteriorMean;
		posteriorVariance += difference * difference;
	}
	double normalizationTerm = unbiased ? (1 / ((double)samples - 1.0)) : (1 / (double)samples);
	return normalizationTerm * posteriorVariance;
}

double FONSEParameter::getMutationVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased)
{
	std::vector<double> mutationParameterTrace = traces.getMutationParameterTraceByMixtureElementForCodon(
		mixtureElement, codon);
	unsigned traceLength = mutationParameterTrace.size();
	if (samples > traceLength)
	{
		std::cerr << "Warning in FONSEParameter::getMutationVariance throws: Number of anticipated samples (" << samples
			<< ") is greater than the length of the available trace (" << traceLength << ")."
			<< "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}

	double posteriorMean = getMutationPosteriorMean(mixtureElement, samples, codon);

	double posteriorVariance = 0.0;

	unsigned start = traceLength - samples;
	double difference = 0.0;
	for (unsigned i = start; i < traceLength; i++)
	{
		difference = mutationParameterTrace[i] - posteriorMean;
		posteriorVariance += difference * difference;
	}
	double normalizationTerm = unbiased ? (1 / ((double)samples - 1.0)) : (1 / (double)samples);
	return normalizationTerm * posteriorVariance;
}

void FONSEParameter::proposeCodonSpecificParameter()
{
	for (unsigned k = 0; k < getGroupListSize(); k++)
	{
		std::vector<double> iidProposed;
		std::string aa = getGrouping(k);
		CodonTable *codonTable = CodonTable::getInstance();
		std::vector <unsigned> codonRange = codonTable -> AAToCodonRange(aa, true);//checked
		
		unsigned numCodons = 0; 
		for (unsigned i = 0; i < codonRange.size(); i++)
		{
			numCodons++;
		}
		for (unsigned i = 0u; i < numCodons * (numMutationCategories + numSelectionCategories); i++)
		{
			iidProposed.push_back(randNorm(0.0, 1.0));
		}

		std::vector<double> covaryingNums;
		covaryingNums = covarianceMatrix[codonTable -> AAToAAIndex(aa)].transformIidNumersIntoCovaryingNumbers(
			iidProposed);
		for (unsigned i = 0; i < numMutationCategories; i++)
		{
			for (unsigned j = i * numCodons, l = codonRange[0]; j < (i * numCodons) + numCodons; j++, l++)
			{
				proposedMutationParameter[i][l] = currentMutationParameter[i][l] + covaryingNums[j];
			}
		}
		for (unsigned i = 0; i < numSelectionCategories; i++)
		{
			for (unsigned j = i * numCodons, l = codonRange[0]; j < (i * numCodons) + numCodons; j++, l++)
			{
				proposedSelectionParameter[i][l] = currentSelectionParameter[i][l]
					+ covaryingNums[(numMutationCategories * numCodons) + j];
			}
		}
	}

}

void FONSEParameter::proposeHyperParameters()
{
	Sphi_proposed = std::exp(randNorm(std::log(Sphi), std_sphi));
}

void FONSEParameter::updateCodonSpecificParameter(std::string grouping)
{
	CodonTable *codonTable = CodonTable::getInstance();
	std::vector <unsigned> codonRange = codonTable -> AAToCodonRange(grouping, true); //checked
	unsigned aaIndex = codonTable -> AAToAAIndex(grouping);
	numAcceptForMutationAndSelection[aaIndex]++;

	for (unsigned k = 0u; k < numMutationCategories; k++)
	{
		for (unsigned i = 0; i < 8; i++)
		{
			currentMutationParameter[k][codonRange[i]] = proposedMutationParameter[k][codonRange[i]];
		}
	}
	for (unsigned k = 0u; k < numSelectionCategories; k++)
	{
		for (unsigned i = 0; i < codonRange.size(); i++)
		{
			currentSelectionParameter[k][codonRange[i]] = proposedSelectionParameter[k][codonRange[i]];
		}
	}
}

void FONSEParameter::initMutationSelectionCategoriesR(std::vector<std::string> files, unsigned numCategories,
	std::string paramType)
{
	unsigned value = 0;
	bool check = true;
	if (paramType == "Mutation")
	{
		value = FONSEParameter::dM;
	}
	else if (paramType == "Selection")
	{
		value = FONSEParameter::dOmega;
	}
	else
	{
		std::cerr << "Bad paramType given. Expected \"Mutation\" or \"Selection\".\nFunction not being executed!\n";
		check = false;
	}
	if (files.size() != numCategories) //we have different sizes and need to stop
	{
		std::cerr
			<< "The number of files given and the number of categories given differ. Function will not be executed!\n";
		check = false;
	}

	if (check)
	{
		initMutationSelectionCategories(files, numCategories, value);
	}
}

