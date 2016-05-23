#include "include/CovarianceMatrix.h"
#include "include/SequenceSummary.h"

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif



//--------------------------------------------------//
// ---------- Constructors & Destructors ---------- //
//--------------------------------------------------//


CovarianceMatrix::CovarianceMatrix()
{
	initCovarianceMatrix(2); //TODO: should not do this
}


CovarianceMatrix::CovarianceMatrix(int _numVariates)
{
	initCovarianceMatrix(_numVariates);
}


CovarianceMatrix::CovarianceMatrix(std::vector <double> &matrix)
{
    numVariates = (int)std::sqrt(matrix.size());
    covMatrix = matrix;
    choleskiMatrix.resize(matrix.size(), 0.0);
}


CovarianceMatrix::CovarianceMatrix(const CovarianceMatrix& other)
{
    numVariates = other.numVariates;
    covMatrix = other.covMatrix;
    choleskiMatrix = other.choleskiMatrix;
}


CovarianceMatrix& CovarianceMatrix::operator=(const CovarianceMatrix& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    numVariates = rhs.numVariates;
    covMatrix = rhs.covMatrix;
	choleskiMatrix = rhs.choleskiMatrix;
    return *this;
}


CovarianceMatrix& CovarianceMatrix::operator+(const CovarianceMatrix& rhs)
{
	std::vector<double> cov = rhs.covMatrix;
	for (unsigned i = 0; i < covMatrix.size(); i++)
	{
		covMatrix[i] += cov[i];
	}
	return *this;
}


CovarianceMatrix& CovarianceMatrix::operator*(const double &value)
{
	for (unsigned i = 0; i < covMatrix.size(); i++)
	{
		covMatrix[i] *= value;
	}
	return *this;
}


void CovarianceMatrix::operator*=(const double &value)
{
  for (unsigned i = 0; i < covMatrix.size(); i++)
  {
    covMatrix[i] *= value;
  }
}


bool CovarianceMatrix::operator==(const CovarianceMatrix& other) const
{
    bool match = true;

    if(this->covMatrix != other.covMatrix) { match = false; }
    if(this->choleskiMatrix != other.choleskiMatrix) { match = false; }
    if(this->numVariates != other.numVariates) { match = false; }

    return match;
}


CovarianceMatrix::~CovarianceMatrix()
{
    //dtor
}





//--------------------------------------//
//---------- Matrix Functions ----------//
//--------------------------------------//


void CovarianceMatrix::initCovarianceMatrix(unsigned _numVariates)
{
    numVariates = _numVariates;
    unsigned vectorLength = numVariates * numVariates;
    covMatrix.resize(vectorLength);
    choleskiMatrix.resize(vectorLength);

	double diag_const = 0.01 / (double)numVariates;
    for (unsigned i = 0u; i < vectorLength; i++)
    {
        covMatrix[i] = (i % (numVariates + 1) ? 0.0 : diag_const);
        choleskiMatrix[i] = 0.0;
    }
}


void CovarianceMatrix::setDiag(double val)
{
	for (unsigned i = 0u; i < covMatrix.size(); i++)
	{
		covMatrix[i] = (i % (numVariates + 1) ? covMatrix[i] : val);
	}
}


// adaptation of http://en.wikipedia.org/wiki/Cholesky_decomposition
// http://rosettacode.org/wiki/Cholesky_decomposition#C
void CovarianceMatrix::choleskiDecomposition()
{
    for (int i = 0; i < numVariates; i++)
    {
        for (int j = 0; j < (i + 1); j++)
        {
            double LsubstractSum = 0.0;
            for (int k = 0; k < j; k++)
            {
                LsubstractSum += choleskiMatrix[i * numVariates + k] * choleskiMatrix[j * numVariates + k];
            }
            choleskiMatrix[i * numVariates + j] = (i == j) ? std::sqrt(covMatrix[i * numVariates + i] - LsubstractSum) :
                (1.0 / choleskiMatrix[j * numVariates + j]) * (covMatrix[i * numVariates + j] - LsubstractSum);
        }
    }
}


void CovarianceMatrix::printCovarianceMatrix()
{

    for (int i = 0; i < numVariates * numVariates; i++)
    {
        if (i % numVariates == 0 && i != 0)
            my_print("\n");
        my_print("%\t", covMatrix[i]);
    }

    my_print("\n");
}


void CovarianceMatrix::printCholeskiMatrix()
{
    for (int i = 0; i < numVariates * numVariates; i++)
    {
        if (i % numVariates == 0 && i != 0)
            my_print("\n");
        my_print("%\t", choleskiMatrix[i]);
    }

    my_print("\n");
}


std::vector<double>* CovarianceMatrix::getCovMatrix()
{
    std::vector<double> *ptr = &covMatrix;
    return ptr;
}


int CovarianceMatrix::getNumVariates()
{
    return numVariates;
}


std::vector<double> CovarianceMatrix::transformIidNumersIntoCovaryingNumbers(std::vector <double> iidnumbers)
{
    std::vector<double> covnumbers;
    for (int i = 0; i < numVariates; i++)
    {
        double sum = 0.0;
        for (int k = 0; k < numVariates; k++)
        {
			// testing if [i * numVariates + k] or [k * numVariates + i], first option was default
            sum += choleskiMatrix[k * numVariates + i] * iidnumbers[k];
        }

        covnumbers.push_back(sum);
    }
    return covnumbers;
}


void CovarianceMatrix::calculateSampleCovariance(std::vector<std::vector<std::vector<std::vector<double>>>> codonSpecificParameterTrace, std::string aa, unsigned samples, unsigned lastIteration)
{
	//order of codonSpecificParameterTrace: paramType, category, numparam, samples
	unsigned numParamTypesInModel = codonSpecificParameterTrace.size();
	unsigned numCategoriesInModel = codonSpecificParameterTrace[0].size();

	unsigned start = lastIteration - samples;
	
	unsigned aaStart;
	unsigned aaEnd;
	SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, true);

	unsigned IDX = 0;
	for (unsigned paramType1 = 0; paramType1 < numParamTypesInModel; paramType1++)
	{
		for (unsigned category1 = 0; category1 < numCategoriesInModel; category1++)
		{
			for (unsigned param1 = aaStart; param1 < aaEnd; param1++)
			{
				double mean1 = sampleMean(codonSpecificParameterTrace[paramType1][category1][param1], samples, lastIteration);
				for (unsigned paramType2 = 0; paramType2 < numParamTypesInModel; paramType2++)
				{
					for (unsigned category2 = 0; category2 < numCategoriesInModel; category2++)
					{
						for (unsigned param2 = aaStart; param2 < aaEnd; param2++)
						{
							double mean2 = sampleMean(codonSpecificParameterTrace[paramType2][category2][param2], samples, lastIteration);
							double unscaledSampleCov = 0.0;
							for (unsigned i = start; i < lastIteration; i++)
							{
								unscaledSampleCov += (codonSpecificParameterTrace[paramType1][category1][param1][i] - mean1) * (codonSpecificParameterTrace[paramType2][category2][param2][i] - mean2);
							}
							covMatrix[IDX] = unscaledSampleCov / (samples - 1);
							IDX++;
						}
					}
				}
			}
		}
	}
}


double CovarianceMatrix::sampleMean(std::vector<double> sampleVector, unsigned samples, unsigned lastIteration)
{
	double posteriorMean = 0.0;
	unsigned start = lastIteration - samples;
	for (unsigned i = start; i < lastIteration; i++)
	{
		posteriorMean += sampleVector[i];
	}
	return posteriorMean / (double)samples;
}





// -----------------------------------------------------------------------------------------------------//
// ---------------------------------------- R SECTION --------------------------------------------------//
// -----------------------------------------------------------------------------------------------------//



#ifndef STANDALONE

void CovarianceMatrix::setCovarianceMatrix(SEXP _matrix)
{
  std::vector<double> tmp;
  NumericMatrix matrix(_matrix);
  unsigned numRows = matrix.nrow();
  covMatrix.resize(numRows * numRows, 0.0);
  numVariates = numRows;
 
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
}





//----------------------------------//
//---------- RCPP Module -----------//
//----------------------------------//


RCPP_MODULE(CovarianceMatrix_mod)
{
  class_<CovarianceMatrix>( "CovarianceMatrix" )

        //Constructors & Destructors:
		.constructor("Empty Constructor")



		//Matrix Functions:
		.method("choleskiDecomposition", &CovarianceMatrix::choleskiDecomposition)
		.method("printCovarianceMatrix", &CovarianceMatrix::printCovarianceMatrix)
		.method("printCholeskiMatrix", &CovarianceMatrix::printCholeskiMatrix)
		.method("setCovarianceMatrix", &CovarianceMatrix::setCovarianceMatrix)
		;
}
#endif
