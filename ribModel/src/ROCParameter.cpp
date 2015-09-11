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

const unsigned ROCParameter::dM = 0;
const unsigned ROCParameter::dEta = 1;

ROCParameter::ROCParameter() : Parameter()
{
}


ROCParameter::ROCParameter(std::string filename) :
		Parameter()
{
	initFromRestartFile(filename);
}

#ifndef STANDALONE
ROCParameter::ROCParameter(double sphi, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix, bool splitSer) : Parameter()
{
	maxGrouping = 22;
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
	maxGrouping = 22;
	std::vector<std::vector<unsigned>> thetaKMatrix;
	initParameterSet(sphi, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
	initROCParameterSet();
}
#endif

ROCParameter::ROCParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
		std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer, std::string _mutationSelectionState) :
		Parameter()
{
	initParameterSet(sphi, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
	initROCParameterSet();
}

ROCParameter::~ROCParameter()
{
	//dtor
}

ROCParameter& ROCParameter::operator=(const ROCParameter& rhs)
{
	if (this == &rhs)
		return *this; // handle self assignment

	Parameter::operator=(rhs);

	covarianceMatrix = rhs.covarianceMatrix;
	// proposal bias and std for codon specific parameter
	bias_csp = rhs.bias_csp;
	std_csp = rhs.std_csp;

	currentMutationParameter = rhs.currentMutationParameter;
	proposedMutationParameter = rhs.proposedMutationParameter;

	currentSelectionParameter = rhs.currentSelectionParameter;
	proposedSelectionParameter = rhs.proposedSelectionParameter;

	phiEpsilon = rhs.phiEpsilon;
	phiEpsilon_proposed = rhs.phiEpsilon_proposed;
	Aphi = rhs.Aphi;
	Aphi_proposed = rhs.Aphi_proposed;
	std_Aphi = rhs.std_Aphi;
	numAcceptForAphi = rhs.numAcceptForAphi;
	traces = rhs.traces;
	numAcceptForMutationAndSelection = rhs.numAcceptForMutationAndSelection;

	return *this;
}

void ROCParameter::initROCParameterSet()
{
	groupList = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R", "S", "T", "V", "Y", "Z"};
	// proposal bias and std for codon specific parameter
	bias_csp = 0;
	std_csp.resize(numParam, 0.1);
	numAcceptForMutationAndSelection.resize(maxGrouping, 0u);

	phiEpsilon = 0.1;
	phiEpsilon_proposed = 0.1;
	for (unsigned i = 0; i < getNumPhiGroupings(); i++) {
		Aphi[i] = 0.0;
		Aphi_proposed[i] = 0.0;
		std_Aphi[i] = 0.1;
		Sepsilon[i] = 0.0;
		numAcceptForAphi[i] = 0;
	}

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

std::vector<std::vector<double>> ROCParameter::calculateSelectionCoefficients(unsigned sample, unsigned mixture)
{
	unsigned numGenes = (unsigned)mixtureAssignment.size();
	std::vector<std::vector<double>> selectionCoefficients;
	selectionCoefficients.resize(numGenes);
	CodonTable *codonTable = CodonTable::getInstance();
	for (unsigned i = 0; i < numGenes; i++)
	{
		for (unsigned j = 0; j < getGroupListSize(); j++)
		{
			std::string aa = getGrouping(j);
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
		std::cerr << "Could not open RestartFile.txt to append\n";
		std::exit(1);
	}

	std::ostringstream oss;
	unsigned j;
	oss << ">Aphi:\n";
	for (unsigned i = 0; i < Aphi.size(); i++)
	{
		oss << Aphi[i];
		if ((i + 1) % 10 == 0)
			oss << "\n";
		else
			oss << " ";
	}
	oss << ">Sepsilon:\n";
	for (unsigned i = 0; i < Sepsilon.size(); i++)
	{
		oss << Sepsilon[i];
		if ((i + 1) % 10 == 0)
			oss << "\n";
		else
			oss << " ";
	}
	oss << ">std_Aphi:\n";
	for (unsigned i = 0; i < std_Aphi.size(); i++)
	{
		oss << std_Aphi[i];
		if ((i + 1) % 10 == 0)
			oss << "\n";
		else
			oss << " ";
	}
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

	CodonTable *codonTable = CodonTable::getInstance();
	for (unsigned i = 0; i < groupList.size(); i++)
	{
		std::string aa = groupList[i];
		oss <<">covarianceMatrix:\n" << aa <<"\n";
		CovarianceMatrix m = covarianceMatrix[codonTable -> AAToAAIndex(aa)];
		std::vector<double>* tmp = m.getCovMatrix();
		int size = m.getNumVariates();
		for(unsigned k = 0; k < size * size; k++)
		{
			if (k % size == 0 && k != 0) { oss <<"\n"; }
			oss << tmp->at(k) << "\t";
		}
		oss <<"\n***\n";
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
	CodonTable *codonTable = CodonTable::getInstance();
	std::ifstream input;
	covarianceMatrix.resize(maxGrouping);
	std::vector <double> mat;
	input.open(filename.c_str());
	if (input.fail())
	{
		std::cerr << "Could not open RestartFile.txt to initialzie ROC values\n";
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
			mat.clear();
			cat = 0;
			variableName = tmp.substr(1,tmp.size()-2);
			if (variableName == "covarianceMatrix")
			{
				getline(input,tmp);
				//char aa = tmp[0];
				cat = codonTable -> AAToAAIndex(tmp); // ????
				std::cout << cat <<"\n";
			}
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
				while (iss >> val) {
					std_csp.push_back(val);
				}
			}
			else if (variableName == "Aphi")
			{
				double val;
				iss.str(tmp);
				while (iss >> val) {
					Aphi.push_back(val);
				}
			}
			else if (variableName == "Sepsilon")
			{
				double val;
				iss.str(tmp);
				while (iss >> val) {
					Sepsilon.push_back(val);
				}
			}
			else if (variableName == "std_Aphi")
			{
				double val;
				iss.str(tmp);
				while (iss >> val) {
					std_Aphi.push_back(val);
				}
			}
			else if (variableName == "covarianceMatrix")
			{
				if (tmp == "***") //end of matrix
				{
					CovarianceMatrix CM(mat);
					covarianceMatrix[cat] = CM;
				}
				double val;
				iss.str(tmp);
				while (iss >> val)
				{
					mat.push_back(val);
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


#ifndef STANDALONE
using namespace Rcpp;
void ROCParameter::initCovarianceMatrix(SEXP _matrix, std::string aa)
{
	std::vector<double> tmp;
	NumericMatrix matrix(_matrix);

	for(unsigned i = 0u; i < aa.length(); i++)	aa[i] = (char)std::toupper(aa[i]);

	CodonTable *codonTable = CodonTable::getInstance();
	unsigned aaIndex = codonTable -> AAToAAIndex(aa);
	unsigned numRows = matrix.nrow();
	std::vector<double> covMatrix(numRows * numRows);

	//NumericMatrix stores the matrix by column, not by row. The loop
	//below transposes the matrix when it stores it.
	unsigned index = 0;
	for (unsigned i = 0; i < numRows; i++)
	{
		for(unsigned j = i; j < numRows * numRows; j += numRows, index++)
		{
			covMatrix[index] = matrix[j];
		}
	}
	CovarianceMatrix m(covMatrix);
	m.choleskiDecomposition();
	covarianceMatrix[aaIndex] = m;
}
#endif

CovarianceMatrix& ROCParameter::getCovarianceMatrixForAA(std::string aa)
{
	CodonTable *codonTable = CodonTable::getInstance();
	aa[0] = (char) std::toupper(aa[0]);
	unsigned aaIndex = codonTable -> AAToAAIndex(aa);;
	return covarianceMatrix[aaIndex];
}

//This function seems to be used only for the purpose of R
void ROCParameter::initSelection(std::vector<double> selectionValues, unsigned mixtureElement, std::string aa)
{
	//TODO: seperate out the R wrapper functionality and make the wrapper
	//currentSelectionParameter
	CodonTable *codonTable = CodonTable::getInstance();
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	if (check)
	{
		mixtureElement--;
		int category = getSelectionCategory(mixtureElement);

		aa[0] = (char) std::toupper(aa[0]);
		std::vector <unsigned> codonRange = codonTable -> AAToCodonRange(aa, true); //checked
		for (unsigned i = 0, j = 0; i < codonRange.size(); i++, j++)
		{
			currentSelectionParameter[category][codonRange[i]] = selectionValues[j];
		}
	}
}

//This function seems to be used only for the purpose of R
void ROCParameter::initMutation(std::vector<double> mutationValues, unsigned mixtureElement, std::string aa)
{
	CodonTable *codonTable = CodonTable::getInstance();
	//TODO: seperate out the R wrapper functionality and make the wrapper
	//currentMutationParameter
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	if (check)
	{
		mixtureElement--;

		unsigned category = getMutationCategory(mixtureElement);
		aa[0] = (char) std::toupper(aa[0]);
		std::vector <unsigned> codonRange = codonTable -> AAToCodonRange(aa, true); //checked
		for (unsigned i = 0, j = 0; i < codonRange.size(); i++, j++)
		{
			currentMutationParameter[category][codonRange[i]] = mutationValues[j];
		}
	}
}



void ROCParameter::initMutationCategories(std::vector<std::string> files, unsigned numCategories)
{
	for (unsigned category = 0; category < numCategories; category++)
	{
		//Open the file for the category
		std::ifstream currentFile;
		currentFile.open(files[category].c_str());
		if (currentFile.fail())
		{
			std::cerr << "Error opening file " << category << " to initialize mutation values.\n";
			std::exit(1);
		}

		std::string tmp;
		currentFile >> tmp; //The first line is a header (Amino Acid, Codon, Value, Std_deviation)

		CodonTable *codonTable = CodonTable::getInstance();
		while (currentFile >> tmp)
		{
			//Get the Codon and Index
			std::size_t pos = tmp.find(",", 2); //Amino Acid and a comma will always be the first 2 characters
			std::string codon = tmp.substr(2, pos - 2);
			unsigned codonIndex = codonTable -> codonToIndex(codon, true);

			//get the value to store
			std::size_t pos2 = tmp.find(",", pos + 1);
			//std::cout << tmp.substr(pos + 1, pos2 - pos - 1 ) <<"\n";
			double value = std::atof(tmp.substr(pos + 1, pos2 - pos - 1).c_str());

			currentMutationParameter[category][codonIndex] = value;
			proposedMutationParameter[category][codonIndex] = value;
		}
		currentFile.close();
	} //END OF A CATEGORY/FILE
}


void ROCParameter::initSelectionCategories(std::vector<std::string> files, unsigned numCategories)
{
	for (unsigned category = 0; category < numCategories; category++)
	{
		//Open the file for the category
		std::ifstream currentFile;
		currentFile.open(files[category].c_str());
		if (currentFile.fail())
		{
			std::cerr << "Error opening file " << category << " to initialize mutation values.\n";
			std::exit(1);
		}

		std::string tmp;
		currentFile >> tmp; //The first line is a header (Amino Acid, Codon, Value, Std_deviation)

		while (currentFile >> tmp)
		{
			//Get the Codon and Index
			std::size_t pos = tmp.find(",", 2); //Amino Acid and a comma will always be the first 2 characters
			std::string codon = tmp.substr(2, pos - 2);
			CodonTable *codonTable = CodonTable::getInstance();
			unsigned codonIndex = codonTable -> codonToIndex(codon, true);

			//get the value to store
			std::size_t pos2 = tmp.find(",", pos + 1);
		//	std::cout << tmp.substr(pos + 1, pos2 - pos - 1 ) <<"\n";
			double value = std::atof(tmp.substr(pos + 1, pos2 - pos - 1).c_str());

			currentSelectionParameter[category][codonIndex] = value;
			proposedSelectionParameter[category][codonIndex] = value;
		}
		currentFile.close();
	} //END OF A CATEGORY/FILE
}

void ROCParameter::setNumPhiGroupings(unsigned _phiGroupings)
{
	phiGroupings = _phiGroupings;
	Aphi.resize(phiGroupings, 0.0);
	Aphi_proposed.resize(phiGroupings, 0.0);
	std_Aphi.resize(phiGroupings, 0.1);
	numAcceptForAphi.resize(phiGroupings, 0);
	Sepsilon.resize(phiGroupings, 0.0);
}

void ROCParameter::getParameterForCategory(unsigned category, unsigned paramType, std::string aa, bool proposal,
		double *returnSet)
{
	std::vector<double> *tempSet;
	if (paramType == ROCParameter::dM)
	{
		tempSet = (proposal ? &proposedMutationParameter[category] : &currentMutationParameter[category]);
	}
	else if (paramType == ROCParameter::dEta)
	{
		tempSet = (proposal ? &proposedSelectionParameter[category] : &currentSelectionParameter[category]);
	}
	else
	{
		std::cerr << "Warning in ROCParameter::getParameterForCategory: Unkown parameter type: " << paramType << "\n";
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

double ROCParameter::getCurrentCodonSpecificProposalWidth(unsigned aa)
{
	CodonTable *codonTable = CodonTable::getInstance();
	std::vector <unsigned> codonRange = codonTable -> AAIndexToCodonRange(aa, true); //checked
	return std_csp[codonRange[0]];
}

void ROCParameter::adaptSphiProposalWidth(unsigned adaptationWidth)
{
	double acceptanceLevel = (double) numAcceptForSphi / (double) adaptationWidth;
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

void ROCParameter::adaptAphiProposalWidth(unsigned adaptationWidth)
{
	for (unsigned i = 0; i < getNumPhiGroupings(); i++) {
		double acceptanceLevel = numAcceptForAphi[i] / (double)adaptationWidth;
		traces.updateAphiAcceptanceRatioTrace(i, acceptanceLevel);
		if (acceptanceLevel < 0.2)
		{
			std_Aphi[i] *= 0.8;
		}
		if (acceptanceLevel > 0.3)
		{
			std_Aphi[i] *= 1.2;
		}
		numAcceptForAphi[i] = 0u;
	}
}

void ROCParameter::adaptSynthesisRateProposalWidth(unsigned adaptationWidth)
{
	unsigned acceptanceUnder = 0u;
	unsigned acceptanceOver = 0u;

	for (unsigned cat = 0u; cat < numSelectionCategories; cat++)
	{
		unsigned numGenes = (unsigned)numAcceptForSynthesisRate[cat].size();
		for (unsigned i = 0; i < numGenes; i++)
		{
			double acceptanceLevel = (double) numAcceptForSynthesisRate[cat][i] / (double) adaptationWidth;
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

void ROCParameter::adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth)
{
	std::cout << "acceptance ratio for amino acid:\n\t";
	for (unsigned i = 0; i < groupList.size(); i++)
	{
		std::cout << groupList[i] << "\t";
	}
	std::cout << "\n\t";

	CodonTable *codonTable = CodonTable::getInstance();
	for (unsigned i = 0; i < groupList.size(); i++)
	{
		std::string aa = groupList[i];
		unsigned aaIndex = codonTable -> AAToAAIndex(aa);
		double acceptanceLevel = (double) numAcceptForMutationAndSelection[aaIndex] / (double) adaptationWidth;
		std::cout << acceptanceLevel << "\t";
		traces.updateCspAcceptanceRatioTrace(aaIndex, acceptanceLevel);
		std::vector <unsigned> codonRange = codonTable -> AAIndexToCodonRange(aaIndex, true); //checked
		for (unsigned k = 0; k < 8; k++)
		{
			if (codonRange[k] == 100) break;
			if (acceptanceLevel < 0.2)
			{
				covarianceMatrix[aaIndex] *= 0.8;
				covarianceMatrix[aaIndex].choleskiDecomposition();
				std_csp[codonRange[k]] *= 0.8;
			}
			if (acceptanceLevel > 0.3)
			{
				covarianceMatrix[aaIndex] *= 1.2;
				covarianceMatrix[aaIndex].choleskiDecomposition();
				std_csp[codonRange[k]] *= 1.2;
			}
		}
		numAcceptForMutationAndSelection[aaIndex] = 0u;
	}
	std::cout << "\n";
}

double ROCParameter::getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement)
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
	unsigned category;
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
	return posteriorMean / (double) usedSamples;
}

double ROCParameter::getSphiPosteriorMean(unsigned samples)
{
	double posteriorMean = 0.0;
	std::vector<double> sPhiTrace = traces.getSphiTrace();
	unsigned traceLength = (unsigned)sPhiTrace.size();

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
	return posteriorMean / (double) samples;
}

double ROCParameter::getAphiPosteriorMean(unsigned index, unsigned samples)
{
	double posteriorMean = 0.0;
	std::vector<double> aPhiTrace = traces.getAphiTrace(index);
	unsigned traceLength = (unsigned)aPhiTrace.size();

	if (samples > traceLength)
	{
		std::cerr << "Warning in Parameter::getAphiPosteriorMean throws: Number of anticipated samples (" << samples
			<< ") is greater than the length of the available trace (" << traceLength << ")."
			<< "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}
	unsigned start = traceLength - samples;
	for (unsigned i = start; i < traceLength; i++)
	{
		posteriorMean += aPhiTrace[i];
	}
	return posteriorMean / (double)samples;
}

std::vector<double> ROCParameter::getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex)
{
	std::vector<unsigned> mixtureAssignmentTrace = traces.getMixtureAssignmentTraceForGene(geneIndex);
	std::vector<double> probabilities(numMixtures, 0.0);
	unsigned traceLength = (unsigned)mixtureAssignmentTrace.size();

	if (samples > traceLength)
	{
		std::cerr
				<< "Warning in ROCParameter::getMixtureAssignmentPosteriorMean throws: Number of anticipated samples ("
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
		probabilities[i] /= (double) samples;
	}
	return probabilities;
}

double ROCParameter::getSphiVariance(unsigned samples, bool unbiased)
{
	std::vector<double> sPhiTrace = traces.getSphiTrace();
	unsigned traceLength = (unsigned)sPhiTrace.size();
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
	double normalizationTerm = unbiased ? (1 / ((double) samples - 1.0)) : (1 / (double) samples);
	return normalizationTerm * posteriorVariance;
}

double ROCParameter::getAphiVariance(unsigned index, unsigned samples, bool unbiased)
{
	std::vector<double> aPhiTrace = traces.getAphiTrace(index);
	unsigned traceLength = (unsigned)aPhiTrace.size();
	if (samples > traceLength)
	{
		std::cerr << "Warning in Parameter::getSphiVariance throws: Number of anticipated samples (" << samples
			<< ") is greater than the length of the available trace (" << traceLength << ")."
			<< "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}
	double posteriorMean = getAphiPosteriorMean(index, samples);

	double posteriorVariance = 0.0;

	unsigned start = traceLength - samples;
	for (unsigned i = start; i < traceLength; i++)
	{
		double difference = aPhiTrace[i] - posteriorMean;
		posteriorVariance += difference * difference;
	}
	double normalizationTerm = unbiased ? (1 / ((double)samples - 1.0)) : (1 / (double)samples);
	return normalizationTerm * posteriorVariance;
}

double ROCParameter::getSynthesisRateVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement,
		bool unbiased)
{
	std::vector<double> synthesisRateTrace = traces.getSynthesisRateTraceByMixtureElementForGene(mixtureElement,
			geneIndex);
	unsigned traceLength = (unsigned)synthesisRateTrace.size();
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
		double difference;
		for (unsigned i = start; i < traceLength; i++)
		{
			difference = synthesisRateTrace[i] - posteriorMean;
			posteriorVariance += difference * difference;
		}
	}
	double normalizationTerm = unbiased ? (1 / ((double) samples - 1.0)) : (1 / (double) samples);
	return normalizationTerm * posteriorVariance;
}

double ROCParameter::getMutationPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon)
{
	double posteriorMean = 0.0;
	std::vector<double> mutationParameterTrace = traces.getMutationParameterTraceByMixtureElementForCodon(
			mixtureElement, codon);
	unsigned traceLength = (unsigned)mutationParameterTrace.size();

	if (samples > traceLength)
	{
		std::cerr << "Warning in ROCParameter::getMutationPosteriorMean throws: Number of anticipated samples ("
				<< samples << ") is greater than the length of the available trace (" << traceLength << ")."
				<< "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}
	unsigned start = traceLength - samples;
	for (unsigned i = start; i < traceLength; i++)
	{
		posteriorMean += mutationParameterTrace[i];
	}
	return posteriorMean / (double) samples;
}

double ROCParameter::getSelectionPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon)
{
	double posteriorMean = 0.0;
	std::vector<double> selectionParameterTrace = traces.getSelectionParameterTraceByMixtureElementForCodon(
			mixtureElement, codon);
	unsigned traceLength = (unsigned)selectionParameterTrace.size();

	if (samples > traceLength)
	{
		std::cerr << "Warning in ROCParameter::getSelectionPosteriorMean throws: Number of anticipated samples ("
				<< samples << ") is greater than the length of the available trace (" << traceLength << ")."
				<< "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}
	unsigned start = traceLength - samples;
	for (unsigned i = start; i < traceLength; i++)
	{
		posteriorMean += selectionParameterTrace[i];
	}
	return posteriorMean / (double) samples;
}

double ROCParameter::getSelectionVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased)
{
	std::vector<double> selectionParameterTrace = traces.getSelectionParameterTraceByMixtureElementForCodon(
			mixtureElement, codon);
	unsigned traceLength = (unsigned)selectionParameterTrace.size();
	if (samples > traceLength)
	{
		std::cerr << "Warning in ROCParameter::getSelectionVariance throws: Number of anticipated samples (" << samples
				<< ") is greater than the length of the available trace (" << traceLength << ")."
				<< "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}

	double posteriorMean = getSelectionPosteriorMean(mixtureElement, samples, codon);

	double posteriorVariance = 0.0;

	unsigned start = traceLength - samples;
	double difference;
	for (unsigned i = start; i < traceLength; i++)
	{
		difference = selectionParameterTrace[i] - posteriorMean;
		posteriorVariance += difference * difference;
	}
	double normalizationTerm = unbiased ? (1 / ((double) samples - 1.0)) : (1 / (double) samples);
	return normalizationTerm * posteriorVariance;
}

double ROCParameter::getMutationVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased)
{
	std::vector<double> mutationParameterTrace = traces.getMutationParameterTraceByMixtureElementForCodon(
			mixtureElement, codon);
	unsigned traceLength = (unsigned)mutationParameterTrace.size();
	if (samples > traceLength)
	{
		std::cerr << "Warning in ROCParameter::getMutationVariance throws: Number of anticipated samples (" << samples
				<< ") is greater than the length of the available trace (" << traceLength << ")."
				<< "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}

	double posteriorMean = getMutationPosteriorMean(mixtureElement, samples, codon);

	double posteriorVariance = 0.0;

	unsigned start = traceLength - samples;
	double difference;
	for (unsigned i = start; i < traceLength; i++)
	{
		difference = mutationParameterTrace[i] - posteriorMean;
		posteriorVariance += difference * difference;
	}
	double normalizationTerm = unbiased ? (1 / ((double) samples - 1.0)) : (1 / (double) samples);
	return normalizationTerm * posteriorVariance;
}

// Cedric: I decided to use a normal distribution to propose Sphi and phi instead of a lognormal because:
// 1. It is a symmetric distribution and you therefore do not have to account for the unsymmetry in jump probabilities
// 2. The one log and exp operation that have to be performed per parameter are cheaper than the operations necessary to draw from a lognormal
// 3. phi has to be on a non log scale for the likelihood evaluation thus it does not help to keep phi on th elog scale all the time
// 4. the adjusment of the likelihood by the jacobian that arises from this transformation is cheap and by grouping everything in one class it takes place more or less at the same place
void ROCParameter::proposeCodonSpecificParameter()
{
	CodonTable *codonTable = CodonTable::getInstance();
	for (unsigned k = 0; k < getGroupListSize(); k++)
	{
		std::vector<double> iidProposed;
		std::string aa = getGrouping(k);
		std::vector <unsigned> codonRange = codonTable -> AAToCodonRange(aa, true); //checked
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

void ROCParameter::proposeAphi()
{
	for (unsigned i = 0; i < getNumPhiGroupings(); i++) {
		Aphi_proposed[i] = randNorm(Aphi[i], std_Aphi[i]);
	}
}

std::vector<double> ROCParameter::propose(std::vector<double> currentParam, double (*proposal)(double a, double b),
		double A, std::vector<double> B)
{
	unsigned _numParam = (unsigned)currentParam.size();
	std::vector<double> proposedParam(_numParam, 0.0);
	for (unsigned i = 0; i < _numParam; i++)
	{
		proposedParam[i] = (*proposal)(A + currentParam[i], B[i]);
	}
	return proposedParam;
}

void ROCParameter::updateCodonSpecificParameter(std::string grouping)
{
	CodonTable *codonTable = CodonTable::getInstance();
	std::vector <unsigned> codonRange = codonTable -> AAToCodonRange(grouping, true); //checked
	unsigned aaIndex = codonTable -> AAToAAIndex(grouping);
	numAcceptForMutationAndSelection[aaIndex]++;

	for (unsigned k = 0u; k < numMutationCategories; k++)
	{
		for (unsigned i = 0; i < codonRange.size(); i++)
		{
			currentMutationParameter[k][codonRange[i]] = proposedMutationParameter[k][codonRange[i]];
		}
	}
	for (unsigned k = 0u; k < numSelectionCategories; k++)
	{
		for (unsigned i = 0; i < codonRange.size(); i++)
		{
			currentSelectionParameter[k][i] = proposedSelectionParameter[k][i];
		}
	}
}



