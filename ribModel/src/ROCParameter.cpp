#include "include/ROCParameter.h"

#include <math.h>
#include <ctime>
#include <iostream>
#include <set>
#include <fstream>

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

ROCParameter::ROCParameter()
{
	std::vector<unsigned> empty(100, 1);
	std::vector<std::vector <unsigned>> empty2;
	initParameterSet(2, 1, empty, empty2, true, "allUnique");
}
#ifndef STANDALONE
ROCParameter::ROCParameter(double sphi, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix, bool splitSer)
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

	initParameterSet( sphi, _numMixtures, geneAssignment, thetaKMatrix, splitSer);

}

ROCParameter::ROCParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, bool splitSer, std::string _mutationSelectionState)
{
	std::vector<std::vector<unsigned>> thetaKMatrix;
	initParameterSet( sphi, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
}
#endif

ROCParameter::ROCParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
		std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer, std::string _mutationSelectionState)
{
	initParameterSet(sphi, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
}
ROCParameter::~ROCParameter()
{
	//dtor
}

ROCParameter::ROCParameter(const ROCParameter& other)
{
	numParam = other.numParam;

	covarianceMatrix = other.covarianceMatrix;
	Sphi = other.Sphi;
	Aphi = other.Aphi;
	Sphi_proposed = other.Sphi_proposed;
	Aphi_proposed = other.Aphi_proposed;
	numAcceptForSphi = other.numAcceptForSphi;
	categories = other.categories;

	// proposal bias and std for phi values
	bias_sphi = other.bias_sphi;
	std_sphi = other.std_sphi;
	prev_std_sphi = other.prev_std_sphi;

	// proposal bias and std for phi values
	bias_phi = other.bias_phi;
	std_phi = other.std_phi;
	prev_std_phi = other.prev_std_phi;

	// proposal bias and std for codon specific parameter
	bias_csp = other.bias_csp;
	std_csp = other.std_csp;
	prev_std_csp = other.prev_std_csp;

	priorA = other.priorA;
	priorB = other.priorB;

	currentExpressionLevel = other.currentExpressionLevel;
	proposedExpressionLevel = other.proposedExpressionLevel;
	numAcceptForExpression = other.numAcceptForExpression;

	currentMutationParameter = other.currentMutationParameter;
	proposedMutationParameter = other.proposedMutationParameter;

	currentSelectionParameter = other.currentSelectionParameter;
	proposedSelectionParameter = other.proposedSelectionParameter;

	numMutationCategories = other.numMutationCategories;
	numSelectionCategories = other.numSelectionCategories;

	phiEpsilon = other.phiEpsilon;
	phiEpsilon_proposed = other.phiEpsilon_proposed;

	numMixtures = other.numMixtures;
}
ROCParameter& ROCParameter::operator=(const ROCParameter& rhs)
{
	if (this == &rhs) return *this; // handle self assignment
	numParam = rhs.numParam;

	covarianceMatrix = rhs.covarianceMatrix;
	Sphi = rhs.Sphi;
	Aphi = rhs.Aphi;
	Sphi_proposed = rhs.Sphi_proposed;
	Aphi_proposed = rhs.Aphi_proposed;
	numAcceptForSphi = rhs.numAcceptForSphi;
	categories = rhs.categories;

	// proposal bias and std for phi values
	bias_sphi = rhs.bias_sphi;
	std_sphi = rhs.std_sphi;
	prev_std_sphi = rhs.prev_std_sphi;

	// proposal bias and std for phi values
	bias_phi = rhs.bias_phi;
	std_phi = rhs.std_phi;
	prev_std_phi = rhs.prev_std_phi;

	// proposal bias and std for codon specific parameter
	bias_csp = rhs.bias_csp;
	std_csp = rhs.std_csp;
	prev_std_csp = rhs.prev_std_csp;

	priorA = rhs.priorA;
	priorB = rhs.priorB;

	currentExpressionLevel = rhs.currentExpressionLevel;
	proposedExpressionLevel = rhs.proposedExpressionLevel;
	numAcceptForExpression = rhs.numAcceptForExpression;

	currentMutationParameter = rhs.currentMutationParameter;
	proposedMutationParameter = rhs.proposedMutationParameter;

	currentSelectionParameter = rhs.currentSelectionParameter;
	proposedSelectionParameter = rhs.proposedSelectionParameter;

	numMutationCategories = rhs.numMutationCategories;
	numSelectionCategories = rhs.numSelectionCategories;

	phiEpsilon = rhs.phiEpsilon;
	phiEpsilon_proposed = rhs.phiEpsilon_proposed;

	numMixtures = rhs.numMixtures;

	return *this;
}

void ROCParameter::initParameterSet(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer, std::string _mutationSelectionState)
{
	// assign genes to mixture element
	unsigned numGenes = geneAssignment.size();
	mixtureAssignment.resize(numGenes, 0);
#ifndef STANDALONE
	//TODO:need to check index are correct, consecutive, and don't exceed numMixtures
	//possibly just use a set?
	for(unsigned i = 0u; i < numGenes; i++)
	{
		mixtureAssignment[i] = geneAssignment[i] - 1;
	}
#else
	for(unsigned i = 0u; i < numGenes; i++)
	{
		mixtureAssignment[i] = geneAssignment[i];
	}
#endif	
	mutationSelectionState = _mutationSelectionState;
	numAcceptForMutationAndSelection.resize(22, 0u);

	numParam = ((splitSer) ? 40 : 41);
	Sphi = sphi;
	Sphi_proposed = sphi;
	bias_sphi = 0;
	std_sphi = 0.1;
	prev_std_sphi = 0.1;

	phiEpsilon = 0.1;
	phiEpsilon_proposed = 0.1;

	numAcceptForSphi = 0u;
	// proposal bias and std for phi values
	bias_phi = 0;

	// proposal bias and std for codon specific parameter
	bias_csp = 0;
	std_csp.resize(numParam, 0.1); //TODO hase to be initialized with 1 when switched to the covariance matrix!!!
	prev_std_csp.resize(numParam, 0.1); //TODO hase to be initialized with 1 when switched to the covariance matrix!!!

	priorA = 0;
	priorB = 1;

	numMixtures = _numMixtures;
	setNumMutationSelectionValues(_mutationSelectionState, thetaKMatrix);
	mutationIsInMixture.resize(numMutationCategories);
	selectionIsInMixture.resize(numSelectionCategories);
	initCategoryDefinitions(_mutationSelectionState, thetaKMatrix);



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
	categoryProbabilities.resize(numMixtures, 1.0/(double)numMixtures);

	currentExpressionLevel.resize(numSelectionCategories);
	proposedExpressionLevel.resize(numSelectionCategories);

	numAcceptForExpression.resize(numSelectionCategories);
	std_phi.resize(numSelectionCategories);
	prev_std_phi.resize(numSelectionCategories);
	for (unsigned i = 0u; i < numSelectionCategories; i++)
	{
		std::vector<double> tmp(numParam, 0.0);
		proposedSelectionParameter[i] = tmp;
		currentSelectionParameter[i] = tmp;

		std::vector<double> tempExpr(numGenes, 0.0);
		currentExpressionLevel[i] = tempExpr;
		proposedExpressionLevel[i] = tempExpr;

		std::vector<unsigned> tempAccExpr(numGenes, 0u);
		numAcceptForExpression[i] = tempAccExpr;

		std::vector<double> tempStdPhi(numGenes, 0.1);
		std_phi[i] = tempStdPhi;
		prev_std_phi[i] = tempStdPhi;
	}
	proposediidSum.resize(22);
	currentiidSum.resize(22);
	for (unsigned i = 0; i < 22; i++)
	{
		char aa = SequenceSummary::AminoAcidArray[i];
		unsigned numCodons = SequenceSummary::GetNumCodonsForAA(aa, true);
		CovarianceMatrix m((numMutationCategories + numSelectionCategories) * numCodons);
		m.choleskiDecomposition();
		covarianceMatrix.push_back(m);
	}
}

//CWO
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
			double phi = getExpressionPosteriorMean(sample, i, mixture);
			for (unsigned k = 0; k < tmp.size(); k++)
			{
				tmp[k] -= minValue;
				selectionCoefficients[i].push_back(phi * tmp[k]);
			}
		}
	}
	return selectionCoefficients;
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
void ROCParameter::initSelection(std::vector<double> selectionValues, unsigned mixtureElement, char aa)
{
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

void ROCParameter::initMutation(std::vector<double> mutationValues, unsigned mixtureElement, char aa)
{
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

#ifndef STANDALONE
using namespace Rcpp;
void ROCParameter::initCovarianceMatrix(SEXP _matrix, char aa)
{
	std::vector<double> tmp;
	NumericMatrix matrix(_matrix);
	aa = std::toupper(aa);
	unsigned aaIndex = SequenceSummary::aaToIndex.find(aa) -> second;
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

CovarianceMatrix& ROCParameter::getCovarianceMatrixForAA(char aa)
{
	aa = std::toupper(aa);
	unsigned aaIndex = SequenceSummary::aaToIndex.find(aa) -> second;
	return covarianceMatrix[aaIndex];
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

void ROCParameter::setNumMutationSelectionValues(std::string mutationSelectionState, std::vector<std::vector<unsigned>> thetaKMatrix)
{
	if (!thetaKMatrix.empty())
	{
		//sets allow only the unique numbers to be added.
		//at the end, the size of the set is equal to the number
		//of unique categories.
		std::set<unsigned> delMCounter;
		std::set<unsigned> delEtaCounter;

		for (unsigned i = 0u; i < numMixtures; i++)
		{
			delMCounter.insert(thetaKMatrix[i][0] - 1);
			delEtaCounter.insert(thetaKMatrix[i][1] - 1);
		}
		numMutationCategories = delMCounter.size();
		numSelectionCategories = delEtaCounter.size();
	}
	else if (mutationSelectionState == selectionShared)
	{
		numMutationCategories = numMixtures;
		numSelectionCategories = 1u;
	}
	else if (mutationSelectionState == mutationShared)
	{
		numMutationCategories = 1u;
		numSelectionCategories = numMixtures;
	}
	else //assuming the default of allUnique
	{
		numMutationCategories = numMixtures;
		numSelectionCategories = numMixtures;
	}
}

void ROCParameter::initCategoryDefinitions(std::string mutationSelectionState, std::vector<std::vector<unsigned>> thetaKMatrix)
{
	std::set<unsigned> delMCounter;
	std::set<unsigned> delEtaCounter;

	for (unsigned i = 0; i < numMixtures; i++)
	{
		categories.push_back(thetaK()); //push a blank thetaK on the vector, then alter.
		if (!thetaKMatrix.empty())
		{
			categories[i].delM = thetaKMatrix[i][0] - 1;
			categories[i].delEta = thetaKMatrix[i][1] - 1; //need check for negative and consecutive checks
			mutationIsInMixture[thetaKMatrix[i][0] - 1].push_back(i);
			selectionIsInMixture[thetaKMatrix[i][1] - 1].push_back(i);
		}
		else if (mutationSelectionState == selectionShared)
		{
			categories[i].delM = i;
			categories[i].delEta = 0;
			mutationIsInMixture[i].push_back(i);
			selectionIsInMixture[0].push_back(i);
		}
		else if (mutationSelectionState == mutationShared)
		{
			categories[i].delM = 0;
			categories[i].delEta = i;
			mutationIsInMixture[0].push_back(i);
			selectionIsInMixture[i].push_back(i);
		}
		else //assuming the default of allUnique
		{
			categories[i].delM = i;
			categories[i].delEta = i;
			mutationIsInMixture[i].push_back(i);
			selectionIsInMixture[i].push_back(i);
		}
		delMCounter.insert(categories[i].delM);
		delEtaCounter.insert(categories[i].delEta);
	}
	//	numMutationCategories = delMCounter.size();
	//	numSelectionCategories = delEtaCounter.size();
	//sets allow only the unique numbers to be added.
	//at the end, the size of the set is equal to the number
	//of unique categories.
}



// calculate SCUO values according to
// Wan et al. CodonO: a new informatics method for measuring synonymous codon usage bias within and across genomes
// International Journal of General Systems, Vol. 35, No. 1, February 2006, 109–125
// http://www.tandfonline.com/doi/pdf/10.1080/03081070500502967
double ROCParameter::calculateSCUO(Gene& gene)
{
	SequenceSummary seqsum = gene.getSequenceSummary();

	double totalDegenerateAACount = 0.0;
	for(int i = 0; i < 22; i++)
	{
		char curAA = seqsum.AminoAcidArray[i];
		// skip amino acids with only one codon or stop codons
		if(curAA == 'X' || curAA == 'M' || curAA == 'W') continue;
		totalDegenerateAACount += (double)seqsum.getAAcountForAA(i);
	}

	double scuoValue = 0.0;
	for(int i = 0; i < 22; i++)
	{
		char curAA = seqsum.AminoAcidArray[i];
		// skip amino acids with only one codon or stop codons
		if(curAA == 'X' || curAA == 'M' || curAA == 'W') continue;
		double numDegenerateCodons = SequenceSummary::GetNumCodonsForAA(curAA);

		double aaCount = (double)seqsum.getAAcountForAA(i);
		if(aaCount == 0) continue;

		unsigned codonRange[2];
		SequenceSummary::AAindexToCodonRange(i, false, codonRange);

		// calculate -sum(pij log(pij))
		double aaEntropy = 0.0;
		unsigned start = codonRange[0];
		unsigned endd = codonRange[1];
		for(unsigned k = start; k < endd; k++)
		{
			int currCodonCount = seqsum.getCodonCountForCodon(k);
			if(currCodonCount == 0) continue;
			double codonProportion = (double)currCodonCount / aaCount;
			aaEntropy += codonProportion*std::log(codonProportion);
		}
		aaEntropy = -aaEntropy;
		// calculate max entropy -log(1/n_i)
		double maxEntropyForAA = -std::log(1.0 / numDegenerateCodons);
		// get normalized difference in entropy O_i
		double normalizedEntropyDiff = (maxEntropyForAA - aaEntropy) / maxEntropyForAA;

		// calculate the composition ratio F_i
		double compositionRatio = aaCount / totalDegenerateAACount;
		// SCUO is the sum(F_i * O_i) over all aa
		scuoValue += compositionRatio * normalizedEntropyDiff;
	}
	return scuoValue;
}

void ROCParameter::InitializeExpression(Genome& genome, double sd_phi)
{
	unsigned genomeSize = genome.getGenomeSize();
	double* scuoValues = new double[genomeSize]();	
	double* expression = new double[genomeSize]();	
	int* index = new int[genomeSize]();	
	//double scuoValues[genomeSize];
	//double expression[genomeSize];
	//int index[genomeSize];

	for(unsigned i = 0u; i < genomeSize; i++)
	{
		index[i] = i;
		scuoValues[i] = calculateSCUO( genome.getGene(i) );
		expression[i] = ROCParameter::randLogNorm(-(sd_phi * sd_phi) / 2, sd_phi);
	}
	quickSortPair(scuoValues, index, 0, genomeSize);
	quickSort(expression, 0, genomeSize);

	for(unsigned category = 0u; category < numSelectionCategories; category++)
	{
		for(unsigned j = 0u; j < genomeSize; j++)
		{
			currentExpressionLevel[category][j] = expression[index[j]];
			//std::cout << currentExpressionLevel[category][j] <<"\n";
			std_phi[category][j] = 0.1;
			prev_std_phi[category][j] = 0.1;
			numAcceptForExpression[category][j] = 0u;
		}
	}

	delete [] scuoValues;
	delete [] expression;
	delete [] index;
}
void ROCParameter::InitializeExpression(double sd_phi)
{
	unsigned numGenes = currentExpressionLevel[1].size();
	for(unsigned category = 0u; category < numSelectionCategories; category++)
	{
		for(unsigned i = 0u; i < numGenes; i++)
		{
			currentExpressionLevel[category][i] = ROCParameter::randLogNorm(-(sd_phi * sd_phi) / 2, sd_phi);
			std_phi[category][i] = 0.1;
			prev_std_phi[category][i] = 0.1;
			numAcceptForExpression[category][i] = 0u;
		}
	}
}
void ROCParameter::InitializeExpression(std::vector<double> expression)
{
	unsigned numGenes = currentExpressionLevel[0].size();
	for(unsigned category = 0u; category < numSelectionCategories; category++)
	{
		for(unsigned i = 0u; i < numGenes; i++)
		{
			currentExpressionLevel[category][i] = expression[i];
			std_phi[category][i] = 0.1;
			prev_std_phi[category][i] = 0.1;
			numAcceptForExpression[category][i] = 0u;
		}
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


std::vector <double> ROCParameter::getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex)
{
	std::vector<unsigned> mixtureAssignmentTrace = traces.getMixtureAssignmentTraceForGene(geneIndex);
	std::vector<double> probabilities(numMixtures, 0.0);
	unsigned traceLength = mixtureAssignmentTrace.size();

	if(samples > traceLength)
	{
		std::cerr << "Warning in ROCParameter::getMixtureAssignmentPosteriorMean throws: Number of anticipated samples (" <<
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

unsigned ROCParameter::getEstimatedMixtureAssignment(unsigned samples, unsigned geneIndex)
{
	unsigned  rv;
	double value = -1.0;
	std::vector <double> probabilities;
	probabilities = getEstimatedMixtureAssignmentProbabilities(samples, geneIndex);

	for (unsigned i = 0; i < probabilities.size(); i++)
	{
		if (value < probabilities[i])
		{
			value = probabilities[i];
			rv = i;
		}
	}

	return rv;
}


double ROCParameter::getExpressionPosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement)
{
	unsigned expressionCategory = getExpressionCategory(mixtureElement);
	double posteriorMean = 0.0;
	std::vector<double> synthesisRateTrace = traces.getSynthesisRateTraceByMixtureElementForGene(mixtureElement, geneIndex);
	unsigned traceLength = synthesisRateTrace.size();


	if(samples > traceLength)
	{
		std::cerr << "Warning in ROCParameter::getExpressionPosteriorMean throws: Number of anticipated samples (" <<
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
		category = getExpressionCategory(category);
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
		std::cerr << "Warning in ROCParameter::getSphiPosteriorMean throws: Number of anticipated samples (" <<
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


double ROCParameter::getExpressionVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement, bool unbiased)
{
	std::vector<double> synthesisRateTrace = traces.getSynthesisRateTraceByMixtureElementForGene(mixtureElement, geneIndex);
	unsigned traceLength = synthesisRateTrace.size();
	if(samples > traceLength)
	{
		std::cerr << "Warning in ROCParameter::getExpressionVariance throws: Number of anticipated samples (" <<
			samples << ") is greater than the length of the available trace (" << traceLength << ")." << "Whole trace is used for posterior estimate! \n";
		samples = traceLength;
	}

	double posteriorMean = getExpressionPosteriorMean(samples, geneIndex, mixtureElement);

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
			category = getExpressionCategory(category);
			difference = synthesisRateTrace[i] - posteriorMean;
			posteriorVariance += difference * difference;
		}
	}
	double normalizationTerm = unbiased ? (1/((double)samples-1.0)) : (1/(double)samples);
	return normalizationTerm * posteriorVariance;
}


double ROCParameter::getSphiVariance(unsigned samples, bool unbiased)
{
	std::vector<double> sPhiTrace = traces.getSPhiTrace();
	unsigned traceLength = sPhiTrace.size();
	if(samples > traceLength)
	{
		std::cerr << "Warning in ROCParameter::getSphiVariance throws: Number of anticipated samples (" <<
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



void ROCParameter::adaptSphiProposalWidth(unsigned adaptationWidth)
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


void ROCParameter::adaptExpressionProposalWidth(unsigned adaptationWidth)
{
	for(unsigned cat = 0u; cat < numSelectionCategories; cat ++)
	{
		unsigned numGenes = numAcceptForExpression[cat].size();
		for(unsigned i = 0; i < numGenes; i++)
		{
			double acceptanceLevel = (double)numAcceptForExpression[cat][i] / (double)adaptationWidth;
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
			numAcceptForExpression[cat][i] = 0u;
		}
	}
}

void ROCParameter::adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth)
{
	unsigned numCSPsets = numAcceptForMutationAndSelection.size();
	for(unsigned i = 0; i < numCSPsets; i++)
	{
		if(i == 21 || i == 10 || i == 18) continue;
		double acceptanceLevel = (double)numAcceptForMutationAndSelection[i] / (double)adaptationWidth;
		std::cout << SequenceSummary::AminoAcidArray[i] << " acceptance ratio: " << acceptanceLevel << "\n";
		traces.updateCspAcceptanceRatioTrace(i, acceptanceLevel);
		unsigned codonRange[2];
		SequenceSummary::AAindexToCodonRange(i, true, codonRange);
		for(unsigned k = codonRange[0]; k < codonRange[1]; k++)
		{
			if(acceptanceLevel < 0.2)
			{
				prev_std_csp[k] = std_csp[k];
				std_csp[k] = std::max(0.01, std_csp[k] * 0.8);
			}
			if(acceptanceLevel > 0.3)
			{
				prev_std_csp[k] = std_csp[k];
				std_csp[k] = std::min(10.0, std_csp[k] * 1.2);
			}
		}
		if(prev_std_csp[codonRange[0]] != std_csp[codonRange[0]])
		{
			// rescale covariance matrix and decompose again
			covarianceMatrix[i] *= std_csp[codonRange[0]];
			covarianceMatrix[i].choleskiDecomposition();
		}
		numAcceptForMutationAndSelection[i] = 0u;
	}
}

// Cedric: I decided to use a normal distribution to propose Sphi and phi instead of a lognormal because:
// 1. It is a symmetric distribution and you therefore do not have to account for the unsymmetry in jump probabilities
// 2. The one log and exp operation that have to be performed per parameter are cheaper than the operations necessary to draw from a lognormal
// 3. phi has to be on a non log scale for the likelihood evaluation thus it does not help to keep phi on th elog scale all the time
// 4. the adjusment of the likelihood by the jacobian that arises from this transformation is cheap and by grouping everything in one class it takes place more or less at the same place
void ROCParameter::proposeSPhi()
{
	//Sphi_proposed = randLogNorm(Sphi, std_sphi);
	// avoid adjusting probabilities for asymmetry of distribution
	Sphi_proposed = std::exp(randNorm(std::log(Sphi), std_sphi));
}
void ROCParameter::proposeExpressionLevels()
{
	unsigned numExpressionLevels = currentExpressionLevel[0].size();
	for(unsigned category = 0; category < numSelectionCategories; category++)
	{
		for(unsigned i = 0u; i < numExpressionLevels; i++)
		{
			// avoid adjusting probabilities for asymmetry of distribution
			proposedExpressionLevel[category][i] = std::exp( randNorm( std::log(currentExpressionLevel[category][i]) , std_phi[category][i]) );
		}
	}
}

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
		for(unsigned i = 0u; i < 2*numCodons*numMixtures; i++)
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
				proposediidSum[k] += covaryingNums[j] * covaryingNums[j];
			}
		}
		for(unsigned i = 0; i < numSelectionCategories; i++)
		{
			for(unsigned j = i * numCodons, l = aaRange[0]; j < (i * numCodons) + numCodons; j++, l++)
			{
				proposedSelectionParameter[i][l] = currentSelectionParameter[i][l] + covaryingNums[(numMutationCategories * numCodons) + j];
				proposediidSum[k] += covaryingNums[(numMutationCategories * numCodons) + j] * covaryingNums[(numMutationCategories * numCodons) + j];
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



void ROCParameter::updateCodonSpecificParameter(char aa)
{
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
	currentiidSum = proposediidSum;
}

std::vector <double> ROCParameter::readPhiValues(std::string filename)
{
	std::size_t pos, pos2;
	std::ifstream currentFile;
	std::string tmpString;
	std::vector<double> RV;

	currentFile.open(filename);
	if (currentFile.fail())
	{
		std::cerr <<"Error opening file\n";
		std::exit(1);
	}

	currentFile >> tmpString; //trash the first line, no info given.


	while (currentFile >> tmpString)
	{
		pos = tmpString.find(",");
		pos2 = tmpString.find(",", pos + 1);
		if (pos != std::string::npos && pos2 != std::string::npos)
		{
			std::string val = tmpString.substr(pos + 1, pos2 - (pos + 1));
			//RV.push_back(std::stod(val));
			RV.push_back(std::atof(val.c_str()));
		}
	}

	return RV;
}


void ROCParameter::setMixtureAssignmentForGene(unsigned geneIndex, unsigned value) 
{
	bool check = checkIndex(geneIndex, 1, mixtureAssignment.size());
	if (check)
	{
		mixtureAssignment[geneIndex - 1] = value;
	}
}
unsigned ROCParameter::getMutationCategoryForMixture(unsigned mixtureElement)
{
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	return check ? categories[mixtureElement - 1].delM + 1 : 0;
}
unsigned ROCParameter::getSelectionCategoryForMixture(unsigned mixtureElement)
{
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	return check ? categories[mixtureElement - 1].delEta + 1 : 0;
}
unsigned ROCParameter::getExpressionCategoryForMixture(unsigned mixtureElement)
{
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	return check ? categories[mixtureElement - 1].delEta + 1 : 0;
}
bool ROCParameter::checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound)
{
	bool check = false;
	if (lowerbound <= index && index <= upperbound)
	{
		check = true;
	}
	else
	{
		std::cerr <<"Error with Index\nINDEX: " << index <<"\n";
		std::cerr <<"MUST BE BETWEEN " << lowerbound << " & " << upperbound <<"\n";
	}

	return check;

}


/*
 * STATIC FUNCTIONS
 */

const unsigned ROCParameter::dM = 0;
const unsigned ROCParameter::dEta = 1;
const std::string ROCParameter::allUnique = "allUnique";
const std::string ROCParameter::selectionShared = "selectionShared";
const std::string ROCParameter::mutationShared = "mutationShared";
std::default_random_engine ROCParameter::generator(time(NULL));

void ROCParameter::drawIidRandomVector(unsigned draws, double mean, double sd, double (*proposal)(double a, double b), double* randomNumbers)
{
	for(unsigned i = 0u; i < draws; i++)
	{
		randomNumbers[i] = (*proposal)(mean, sd);
	}
}
void ROCParameter::drawIidRandomVector(unsigned draws, double r, double (*proposal)(double r), double* randomNumbers)
{
	for(unsigned i = 0u; i < draws; i++)
	{
		randomNumbers[i] = (*proposal)(r);
	}
}
double ROCParameter::randNorm(double mean, double sd)
{
	double rv;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = rnorm(1, mean, sd);
	rv = xx[0];
#else
	std::normal_distribution<double> distribution(mean, sd);
	rv = distribution(generator);
#endif
	return rv;
}


double ROCParameter::randLogNorm(double m, double s)
{
	double rv;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = rlnorm(1, m, s);
	rv = xx[0];
#else
	std::lognormal_distribution<double> distribution(m, s);
	rv = distribution(generator);
#endif
	return rv;
}


double ROCParameter::randExp(double r)
{
	double rv;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = rexp(1, r);
	rv = xx[0];
#else
	std::exponential_distribution<double> distribution(r);
	rv = distribution(generator);
#endif
	return rv;
}

void ROCParameter::randDirichlet(double* input, unsigned numElements, double* output)
{
	// draw y_i from Gamma(a_i, 1)
	// normalize y_i such that x_i = y_i / sum(y_i)

	double sumTotal = 0.0;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	for(unsigned i = 0; i < numElements; i++)
	{
		xx = rgamma(1, input[i], 1);
		output[i] = xx[0];
		sumTotal += xx[0];
	}
#else
	for(unsigned i = 0; i < numElements; i++)
	{
		std::gamma_distribution<double> distribution(input[i], 1);
		output[i] = distribution(generator);
		sumTotal += output[i];
	}
#endif
	for(unsigned i = 0; i < numElements; i++)
	{
		output[i] = output[i] / sumTotal;
	}
}

double ROCParameter::randUnif(double minVal, double maxVal)
{
	double rv;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = runif(1, minVal, maxVal);
	rv = xx[0];
#else
	std::uniform_real_distribution<double> distribution(minVal, maxVal);
	rv = distribution(generator);	
#endif
	return rv;
}

unsigned ROCParameter::randMultinom(double* probabilities, unsigned mixtureElements)
{
	// calculate cummulative sum to determine group boundaries
	double* cumsum = new double[mixtureElements]();
	//std::vector<double> cumsum(groups);
	cumsum[0] = probabilities[0];

	for(unsigned i = 1u; i < mixtureElements; i++)
	{
		cumsum[i] = cumsum[i-1u] + probabilities[i];
	}
	// draw random number from U(0,1)
	double referenceValue;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = runif(1, 0, 1);
	referenceValue = xx[0];
#else
	std::uniform_real_distribution<double> distribution(0, 1);
	referenceValue = distribution(generator);
#endif
	// check in which category the element falls
	unsigned returnValue = 0u;
	for (unsigned i = 0u; i < mixtureElements; i++)
	{
		if (referenceValue <= cumsum[i])
		{
			returnValue = i;
			break;
		}
	}
	delete [] cumsum;
	return returnValue;
}

double ROCParameter::densityNorm(double x, double mean, double sd)
{
	const double inv_sqrt_2pi = 0.3989422804014327;
	double a = (x - mean) / sd;

	return (inv_sqrt_2pi / sd) * std::exp(-0.5 * a * a);
}
double ROCParameter::densityLogNorm(double x, double mean, double sd)
{
	double returnValue = 0.0;
	// logN is only defined for x > 0 => all values less or equal to zero have probability 0
	if(x > 0.0)
	{
		const double inv_sqrt_2pi = 0.3989422804014327;
		double a = (std::log(x) - mean) / sd;
		returnValue = (inv_sqrt_2pi / (x * sd)) * std::exp(-0.5 * a * a);
	}
	return returnValue;
}

// sort array interval from first (included) to last (excluded)!!
// quick sort, sorting arrays a and b by a.
// Elements in b corespond to a, a will be sorted and it will be assured that b will be sorted by a
void ROCParameter::quickSortPair(double a[], int b[], int first, int last)
{
	int pivotElement;

	if(first < last)
	{
		pivotElement = pivotPair(a, b, first, last);
		quickSortPair(a, b, first, pivotElement - 1);
		quickSortPair(a, b, pivotElement + 1, last);
	}
}
// sort array interval from first (included) to last (excluded)!!
void ROCParameter::quickSort(double a[], int first, int last)
{
	int pivotElement;

	if(first < last)
	{
		pivotElement = pivot(a, first, last);
		quickSort(a, first, pivotElement - 1);
		quickSort(a, pivotElement + 1, last);
	}
}
int ROCParameter::pivot(double a[], int first, int last)
{
	int p = first;
	double pivotElement = a[first];

	for(int i = (first + 1) ; i < last ; i++)
	{
		/* If you want to sort the list in the other order, change "<=" to ">" */
		if(a[i] <= pivotElement)
		{
			p++;
			swap(a[i], a[p]);
		}
	}
	swap(a[p], a[first]);

	return p;
}
int ROCParameter::pivotPair(double a[], int b[], int first, int last)
{
	int p = first;
	double pivotElement = a[first];

	for(int i = (first + 1) ; i < last ; i++)
	{
		/* If you want to sort the list in the other order, change "<=" to ">" */
		if(a[i] <= pivotElement)
		{
			p++;
			swap(a[i], a[p]);
			swap(b[i], b[p]);
		}
	}
	swap(a[p], a[first]);
	swap(b[p], b[first]);

	return p;
}
void ROCParameter::swap(double& a, double& b)
{
	double temp = a;
	a = b;
	b = temp;
}
void ROCParameter::swap(int& a, int& b)
{
	double temp = a;
	a = b;
	b = temp;
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



#ifndef STANDALONE

	RCPP_EXPOSED_CLASS(ROCTrace)
	RCPP_EXPOSED_CLASS(Genome)
RCPP_EXPOSED_CLASS(CovarianceMatrix)

RCPP_MODULE(ROCParameter_mod)
{
	class_<ROCParameter>( "ROCParameter" )
		.constructor("empty constructor")
		.constructor <double, std::vector<unsigned>, std::vector<unsigned>, bool>()
		.constructor <double, unsigned, std::vector<unsigned>, bool, std::string>()
		.method("initMutationSelectionCategories", &ROCParameter::initMutationSelectionCategoriesR)
		.method("readPhiValues", &ROCParameter::readPhiValues)
		.method("setMixtureAssignmentForGene", &ROCParameter::setMixtureAssignmentForGene)
		.method("getMutationCategoryForMixture", &ROCParameter::getMutationCategoryForMixture)
		.method("getSelectionCategoryForMixture", &ROCParameter::getSelectionCategoryForMixture)
		.method("getExpressionCategoryForMixture", &ROCParameter::getExpressionCategoryForMixture)

		.method("initSelection", &ROCParameter::initSelection)
		.method("initMutation", &ROCParameter::initMutation)
		.method("initCovarianceMatrix", &ROCParameter::initCovarianceMatrix)
		.method("getCovarianceMatrixForAA", &ROCParameter::getCovarianceMatrixForAA)
		.method("getTraceObject", &ROCParameter::getTraceObject)

		//R wrapper functions
		.method("initializeExpressionByGenome", &ROCParameter::initializeExpressionByGenome)
		.method("initializeExpressionByList", &ROCParameter::initializeExpressionByList)
		.method("initializeExpressionByRandom", &ROCParameter::initializeExpressionByRandom)
		.method("getEstimatedMixtureAssignmentForGene", &ROCParameter::getEstimatedMixtureAssignmentForGene, "returns the mixture assignment for a given gene")
		.method("getEstimatedMixtureAssignmentProbabilitiesForGene", &ROCParameter::getEstimatedMixtureAssignmentProbabilitiesForGene, "returns the probabilities assignment for a given gene")	
		.method("calculateSelectionCoefficients", &ROCParameter::calculateSelectionCoefficientsR)
		// Posterior functions
		.method("getSphiPosteriorMean", &ROCParameter::getSphiPosteriorMean)
		//.method("getMixtureAssignmentPosteriorMean", &ROCParameter::getMixtureAssignmentPosteriorMeanR)

		//R wrapper functions
		.method("getExpressionPosteriorMeanByMixtureElementForGene", &ROCParameter::getExpressionPosteriorMeanByMixtureElementForGene)
		.method("getMutationPosteriorMeanForCodon", &ROCParameter::getMutationPosteriorMeanForCodon)
		.method("getSelectionPosteriorMeanForCodon", &ROCParameter::getSelectionPosteriorMeanForCodon)

		// Variance functions
		.method("getSphiVariance", &ROCParameter::getSphiVariance)

		//R wrapper functions
		.method("getMutationVarianceForCodon", &ROCParameter::getMutationVarianceForCodon)
		.method("getSelectionVarianceForCodon", &ROCParameter::getSelectionVarianceForCodon)
		.method("getExpressionVarianceByMixtureElementForGene", &ROCParameter::getExpressionVarianceByMixtureElementForGene)
		.method("getCurrentExpressionForMixture", &ROCParameter::getCurrentExpressionForMixture)

		.property("numMutationCategories", &ROCParameter::getNumMutationCategories)
		.property("numSelectionCategories", &ROCParameter::getNumSelectionCategories)
		.property("numMixtures", &ROCParameter::getNumMixtureElements)
		;
	//	function("randUnif", &ROCParameter::randUnif);
	//	function("randExp", &ROCParameter::randExp);
	//	function("randLogNorm", &ROCParameter::randLogNorm);
	//	function("randNorm", &ROCParameter::randNorm);
}
#endif
