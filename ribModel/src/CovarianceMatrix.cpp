#include "include/CovarianceMatrix.h"
#include <iostream>
#include <cmath>


CovarianceMatrix::CovarianceMatrix()
{
	initCovarianceMatrix(2);

}
CovarianceMatrix::CovarianceMatrix(int _numVariates)
{
	initCovarianceMatrix(_numVariates);
}
CovarianceMatrix::CovarianceMatrix(std::vector <double> &matrix)
{
    numVariates = std::sqrt(matrix.size());
    covMatrix = matrix;
    choleskiMatrix.resize(matrix.size(), 0.0);
}


CovarianceMatrix::~CovarianceMatrix()
{
    //dtor
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

void CovarianceMatrix::operator*(const double &value)
{
	for (unsigned i = 0; i < covMatrix.size(); i++)
	{
		covMatrix[i] *= value;
	}
}


void CovarianceMatrix::operator*=(const double &value)
{
  for (unsigned i = 0; i < covMatrix.size(); i++)
  {
    covMatrix[i] *= value;
  }
}


void CovarianceMatrix::initCovarianceMatrix(unsigned _numVariates)
{
    numVariates = _numVariates;
    unsigned vectorLength = numVariates * numVariates;
    covMatrix.resize(vectorLength);
    choleskiMatrix.resize(vectorLength);

    for(unsigned i = 0u; i < vectorLength; i++)
    {
        covMatrix[i] = (i % (numVariates + 1) ? 0.0 : 1.0);
        choleskiMatrix[i] = 0.0;
    }
}


// addaptatoin of http://en.wikipedia.org/wiki/Cholesky_decomposition
// http://rosettacode.org/wiki/Cholesky_decomposition#C
void CovarianceMatrix::choleskiDecomposition()
{
    std::cout <<"choleskiDecomposition called\n";
    for(int i = 0; i < numVariates; i++)
    {
        for(int j = 0; j < (i + 1); j++)
        {
            double LsubstractSum = 0.0;
            for(int k = 0; k < j; k++)
            {
                LsubstractSum += choleskiMatrix[i * numVariates + k] * choleskiMatrix[j * numVariates + k];
            }
            choleskiMatrix[i * numVariates + j] = (i == j) ? std::sqrt(covMatrix[i * numVariates + i] - LsubstractSum) :
                (1.0 / choleskiMatrix[j * numVariates + j]) * (covMatrix[i * numVariates + j] - LsubstractSum);
        }
    }
}


void CovarianceMatrix::calculateCovarianceMatrixFromTraces(std::vector<std::vector <std::vector<double>>> trace, unsigned geneIndex, unsigned curSample, unsigned adaptiveWidth)
{
    // calculate all means
    unsigned numVariates = trace.size(); // <- number of mixture elements or number of selection categories
    double* means = new double[numVariates]();
    unsigned start = curSample - adaptiveWidth;
    // calculate all means from trace
    for(unsigned i = 0u; i < numVariates; i++)
    {
        for(unsigned j = start; j < curSample; j++)
        {
            means[i] += trace[i][j][geneIndex];
        }
        means[i] /= adaptiveWidth;
    }


    for(unsigned i = 0u; i < numVariates; i++)
    {
        for (unsigned k = 0u; k < numVariates; k++)
        {
            double nonNormalizedCovariance = 0.0; // missing term 1/(n-1)
            for(unsigned j = start; j < curSample; j++)
            {
                nonNormalizedCovariance += (trace[i][j][geneIndex] - means[i]) * (trace[k][j][geneIndex] - means[k]);
            }
            covMatrix[i * numVariates + k] = (1.0/(adaptiveWidth - 1.0)) * nonNormalizedCovariance;
        }

    }
}
std::vector<double> CovarianceMatrix::transformIidNumersIntoCovaryingNumbers(std::vector <double> iidnumbers)
{
    std::vector<double> covnumbers;
    for(int i = 0; i < numVariates; i++)
    {
        double sum = 0.0;
        for (int k = 0; k < numVariates; k++)
        {
            sum += choleskiMatrix[i * numVariates + k] * iidnumbers[k];
        }

        covnumbers.push_back(sum);
    }
    return covnumbers;
}

void CovarianceMatrix::printCovarianceMatrix()
{
    for(int i = 0; i < numVariates * numVariates; i++)
    {
        if (i % numVariates == 0 && i != 0) { std::cout << std::endl; }
        std::cout << covMatrix[i]<< "\t";
    }
    std::cout <<"\n";


}
void CovarianceMatrix::printCholeskiMatrix()
{
    for(int i = 0; i < numVariates * numVariates; i++)
    {
        if (i % numVariates == 0 && i != 0) { std::cout << std::endl; }
        std::cout << choleskiMatrix[i]<< "\t";
    }
    std::cout <<"\n";
}

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
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
#endif
#ifndef STANDALONE

RCPP_MODULE(CovarianceMatrix_mod)
{
  class_<CovarianceMatrix>( "CovarianceMatrix" )
		.constructor("Empty Constructor")
		.method("choleskiDecomposition", &CovarianceMatrix::choleskiDecomposition)
		.method("printCovarianceMatrix", &CovarianceMatrix::printCovarianceMatrix)
		.method("printCholeskiMatrix", &CovarianceMatrix::printCholeskiMatrix)
		.method("setCovarianceMatrix", &CovarianceMatrix::setCovarianceMatrix)
		;
}
#endif
