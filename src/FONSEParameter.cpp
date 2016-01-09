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

//--------------------------------------------------//
// ---------- Constructors & Destructors ---------- //
//--------------------------------------------------//


FONSEParameter::FONSEParameter() : Parameter()
{
	//CTOR
	bias_csp = 0;
	mutation_prior_sd = 0.35;
}


FONSEParameter::FONSEParameter(std::string filename) : Parameter(22)
{
	initFromRestartFile(filename);
}


FONSEParameter::FONSEParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
	std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer, std::string _mutationSelectionState) :
	Parameter(22)
{
	initParameterSet(stdDevSynthesisRate, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
	initFONSEParameterSet();
}


FONSEParameter& FONSEParameter::operator=(const FONSEParameter& rhs)
{
	if (this == &rhs)
		return *this; // handle self assignment

	Parameter::operator=(rhs);

	// proposal bias and std for codon specific parameter
	bias_csp = rhs.bias_csp;
	std_csp = rhs.std_csp;

	mutation_prior_sd = rhs.mutation_prior_sd;

	currentMutationParameter = rhs.currentMutationParameter;
	proposedMutationParameter = rhs.proposedMutationParameter;

	currentSelectionParameter = rhs.currentSelectionParameter;
	proposedSelectionParameter = rhs.proposedSelectionParameter;

	return *this;
}


FONSEParameter::FONSEParameter(const FONSEParameter &other) : Parameter(other)
{
	bias_csp = other.bias_csp;
	std_csp = other.std_csp;

	mutation_prior_sd = other.mutation_prior_sd;

	currentMutationParameter = other.currentMutationParameter;
	proposedMutationParameter = other.proposedMutationParameter;

	currentSelectionParameter = other.currentSelectionParameter;
	proposedSelectionParameter = other.proposedSelectionParameter;


}


FONSEParameter::~FONSEParameter()
{
	// destructor
}





//---------------------------------------------------------------//
// ---------- Initialization, Restart, Index Checking -----------//
//---------------------------------------------------------------//


void FONSEParameter::initFONSEParameterSet()
{
	mutation_prior_sd = 0.35;
	groupList = { "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R", "S", "T", "V", "Y", "Z" };
	// proposal bias and std for codon specific parameter
	bias_csp = 0;
	std_csp.resize(numParam, 0.1);
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

	for (unsigned i = 0; i < maxGrouping; i++)
	{
		std::string aa = SequenceSummary::AminoAcidArray[i];
		unsigned numCodons = SequenceSummary::GetNumCodonsForAA(aa, true);
		CovarianceMatrix m((numMutationCategories + numSelectionCategories) * numCodons);
		m.choleskiDecomposition();
		covarianceMatrix.push_back(m);
	}

}


void FONSEParameter::initFONSEValuesFromFile(std::string filename)
{
	std::ifstream input;
	std::vector <double> mat;
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
			mat.clear();
			cat = 0;
			variableName = tmp.substr(1, tmp.size() - 2);
			if (variableName == "covarianceMatrix")
			{
				getline(input, tmp);
				//char aa = tmp[0];
				cat = SequenceSummary::AAToAAIndex(tmp); // ????
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
			else if (variableName == "std_csp")
			{
				double val;
				iss.str(tmp);
				while (iss >> val)
				{
					std_csp.push_back(val);
				}
			}
			else if (variableName == "mutation_prior_sd")
			{
				iss.str(tmp);
				iss >> mutation_prior_sd;
			}
		}
	}
	input.close();

	//init other values
	bias_csp = 0;
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
	oss << ">mutation_prior_sd:\n" << mutation_prior_sd << "\n";
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
	for (unsigned i = 0; i < groupList.size(); i++)
	{
		std::string aa = groupList[i];
		oss << ">covarianceMatrix:\n" << aa << "\n";
		CovarianceMatrix m = covarianceMatrix[SequenceSummary::AAToAAIndex(aa)];
		std::vector<double>* tmp = m.getCovMatrix();
		int size = m.getNumVariates();
		for (unsigned k = 0; k < size * size; k++)
		{
			if (k % size == 0 && k != 0) { oss << "\n"; }
			oss << tmp->at(k) << "\t";
		}
		oss << "\n***\n";
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


void FONSEParameter::initAllTraces(unsigned samples, unsigned num_genes)
{
    traces.initializeFONSETrace(samples, num_genes, numMutationCategories, numSelectionCategories, numParam,
                         numMixtures, categories, maxGrouping);
}

void FONSEParameter::initMutationCategories(std::vector<std::string> files, unsigned numCategories)
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
            unsigned codonIndex = SequenceSummary::codonToIndex(codon, true);
            
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


void FONSEParameter::initSelectionCategories(std::vector<std::string> files, unsigned numCategories)
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
            unsigned codonIndex = SequenceSummary::codonToIndex(codon, true);
            
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






// --------------------------------------//
// ---------- Trace Functions -----------//
// --------------------------------------//


void FONSEParameter::updateCodonSpecificParameterTrace(unsigned sample, std::string grouping)
{
	traces.updateCodonSpecificParameterTraceForAA(sample, grouping, currentMutationParameter, dM);
	traces.updateCodonSpecificParameterTraceForAA(sample, grouping, currentSelectionParameter, dOmega);
}


// ------------------------------------------//
// ---------- Covariance Functions ----------//
// ------------------------------------------//


CovarianceMatrix& FONSEParameter::getCovarianceMatrixForAA(std::string aa)
{
    aa[0] = (char)std::toupper(aa[0]);
    unsigned aaIndex = SequenceSummary::aaToIndex.find(aa)->second;
    return covarianceMatrix[aaIndex];
}


// -----------------------------------//
// ---------- CSP Functions ----------//
// -----------------------------------//


double FONSEParameter::getCurrentCodonSpecificProposalWidth(unsigned aa)
{
	unsigned aaStart;
	unsigned aaEnd;
	SequenceSummary::AAIndexToCodonRange(aa, aaStart, aaEnd, false);
    return std_csp[aaStart];
}


void FONSEParameter::proposeCodonSpecificParameter()
{
    for (unsigned k = 0; k < getGroupListSize(); k++)
    {
        std::vector<double> iidProposed;
        std::string aa = getGrouping(k);
		unsigned aaStart;
		unsigned aaEnd;
		SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);
        unsigned numCodons = aaEnd - aaStart;
        for (unsigned i = 0u; i < numCodons * (numMutationCategories + numSelectionCategories); i++)
        {
            iidProposed.push_back(randNorm(0.0, 1.0));
        }
        
        std::vector<double> covaryingNums;
        covaryingNums = covarianceMatrix[SequenceSummary::AAToAAIndex(aa)].transformIidNumersIntoCovaryingNumbers(iidProposed);
        for (unsigned i = 0; i < numMutationCategories; i++)
        {
            for (unsigned j = i * numCodons, l = aaStart; j < (i * numCodons) + numCodons; j++, l++)
            {
                proposedMutationParameter[i][l] = currentMutationParameter[i][l] + covaryingNums[j];
            }
        }
        for (unsigned i = 0; i < numSelectionCategories; i++)
        {
            for (unsigned j = i * numCodons, l = aaStart; j < (i * numCodons) + numCodons; j++, l++)
            {
                proposedSelectionParameter[i][l] = currentSelectionParameter[i][l]
                + covaryingNums[(numMutationCategories * numCodons) + j];
            }
        }
    }
}


void FONSEParameter::updateCodonSpecificParameter(std::string grouping)
{
	unsigned aaStart;
	unsigned aaEnd;
	SequenceSummary::AAToCodonRange(grouping, aaStart, aaEnd, true);
    unsigned aaIndex = SequenceSummary::aaToIndex.find(grouping)->second;
	numAcceptForCodonSpecificParameters[aaIndex]++;
    
    for (unsigned k = 0u; k < numMutationCategories; k++)
    {
        for (unsigned i = aaStart; i < aaEnd; i++)
        {
            currentMutationParameter[k][i] = proposedMutationParameter[k][i];
        }
    }
    for (unsigned k = 0u; k < numSelectionCategories; k++)
    {
        for (unsigned i = aaStart; i < aaEnd; i++)
        {
            currentSelectionParameter[k][i] = proposedSelectionParameter[k][i];
        }
    }
}

// -------------------------------------//
// ---------- Other Functions ----------//
// -------------------------------------//

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
	unsigned aaStart;
	unsigned aaEnd;
	SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);
    
    unsigned j = 0u;
    for (unsigned i = aaStart; i < aaEnd; i++, j++)
    {
        returnSet[j] = tempSet->at(i);
    }
}


void FONSEParameter::proposeHyperParameters()
{
    for(unsigned i = 0u; i < numSelectionCategories; i++)
    {
		stdDevSynthesisRate_proposed[i] = std::exp(randNorm(std::log(stdDevSynthesisRate[i]), std_stdDevSynthesisRate));
    }
}








// -----------------------------------------------------------------------------------------------------//
// ---------------------------------------- R SECTION --------------------------------------------------//
// -----------------------------------------------------------------------------------------------------//


#ifndef STANDALONE


//--------------------------------------------------//
// ---------- Constructors & Destructors ---------- //
//--------------------------------------------------//


FONSEParameter::FONSEParameter(std::vector<double> stdDevSynthesisRate, std::vector<unsigned> geneAssignment,
                               std::vector<unsigned> _matrix, bool splitSer) : Parameter(22)
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
    initParameterSet(stdDevSynthesisRate, _matrix.size() / 2, geneAssignment, thetaKMatrix, splitSer);
    initFONSEParameterSet();
    
}


FONSEParameter::FONSEParameter(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
                               bool splitSer, std::string _mutationSelectionState) : Parameter(22)
{
    std::vector<std::vector<unsigned>> thetaKMatrix;
    initParameterSet(stdDevSynthesisRate, _numMixtures, geneAssignment, thetaKMatrix, splitSer, _mutationSelectionState);
    initFONSEParameterSet();
}





//---------------------------------------------------------------//
// ---------- Initialization, Restart, Index Checking ---------- //
//---------------------------------------------------------------//


void FONSEParameter::initCovarianceMatrix(SEXP _matrix, std::string aa)
{
    std::vector<double> tmp;
    NumericMatrix matrix(_matrix);
    
    for (unsigned i = 0u; i < aa.length(); i++)	aa[i] = (char)std::toupper(aa[i]);
    
    unsigned aaIndex = SequenceSummary::aaToIndex.find(aa)->second;
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

void FONSEParameter::initMutation(std::vector<double> mutationValues, unsigned mixtureElement, std::string aa)
{
    //TODO: seperate out the R wrapper functionality and make the wrapper
    //currentMutationParameter
    bool check = checkIndex(mixtureElement, 1, numMixtures);
    if (check)
    {
        mixtureElement--;
        
        
        unsigned category = getMutationCategory(mixtureElement);
        aa[0] = (char)std::toupper(aa[0]);
	unsigned aaStart;
	unsigned aaEnd;
	SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);
        for (unsigned i = aaStart, j = 0; i < aaEnd; i++, j++)
        {
            currentMutationParameter[category][i] = mutationValues[j];
        }
    }
}


void FONSEParameter::initSelection(std::vector<double> selectionValues, unsigned mixtureElement, std::string aa)
{
    //TODO: seperate out the R wrapper functionality and make the wrapper
    //currentSelectionParameter
    bool check = checkIndex(mixtureElement, 1, numMixtures);
    if (check)
    {
        mixtureElement--;
        
        int category = getSelectionCategory(mixtureElement);
        
        aa[0] = (char)std::toupper(aa[0]);
	unsigned aaStart;
	unsigned aaEnd;
	SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);
        for (unsigned i = aaStart, j = 0; i < aaEnd; i++, j++)
        {
            currentSelectionParameter[category][i] = selectionValues[j];
        }
    }
}


// -----------------------------------//
// ---------- CSP Functions ----------//
// -----------------------------------//


std::vector< std::vector <double> > FONSEParameter::getCurrentMutationParameter()
{
    return currentMutationParameter;
}


std::vector< std::vector <double> > FONSEParameter::getCurrentSelectionParameter()
{
    return currentSelectionParameter;
}

#endif
