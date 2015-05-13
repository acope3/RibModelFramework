#include "../include/My1DCovarianceMatrix.h"
#include <iostream>
#include <cmath>


CovarianceMatrix::CovarianceMatrix()
{
		std::cout <<"0 arg constructor called\n";
		numVariates = 4; //Square 2x2 matrix
		covMatrix.resize(numVariates * numVariates);
		choleskiMatrix.resize(numVariates * numVariates);

    for(int i = 0; i < numVariates * numVariates; i++)
    {
            covMatrix[i] = (i % (numVariates + 1) ? 0.0 : 1.0);
            choleskiMatrix[i] = 0.0;
		}

}
CovarianceMatrix::CovarianceMatrix(int _numVariates)
{

		std::cout <<"1 arg constructor called\n";
		numVariates = _numVariates; 
		covMatrix.resize(numVariates * numVariates);
		choleskiMatrix.resize(numVariates * numVariates);

    for(int i = 0; i < numVariates * numVariates; i++)
    {
            covMatrix[i] = (i % (numVariates + 1) ? 0.0 : 1.0);
            choleskiMatrix[i] = 0.0;
		}

}

CovarianceMatrix::CovarianceMatrix(std::vector <double> &matrix)
{
		std::cout <<"Matrix constructor called\n";
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
		std::cout <<"copy constructor called\n";
    numVariates = other.numVariates;
		covMatrix = other.covMatrix;
		choleskiMatrix = other.choleskiMatrix;
}
/*
CovarianceMatrix& CovarianceMatrix::operator=(const CovarianceMatrix& rhs)
{
		std::cout <<"op over constructor called\n";
    if (this == &rhs) return *this; // handle self assignment
    covMatrix[rhs.numVariates][rhs.numVariates];
    numVariates = rhs.numVariates;
    for(int i = 0; i < rhs.numVariates; i++)
    {
        for(int j = 0; j < rhs.numVariates; j++)
        {
            covMatrix[i][j] = rhs.covMatrix[i][j];
        }
    }
    return *this;
}
*/
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

void CovarianceMatrix::transformIidNumersIntoCovaryingNumbers(double* iidnumbers, double* covnumbers)
{
    //double covnumbers[numVariates];
    for(int i = 0; i < numVariates; i++)
    {
        double sum = 0.0;
        for (int k = 0; k < numVariates; k++)
        {
            sum += choleskiMatrix[i * numVariates + k] * iidnumbers[k];
        }

        covnumbers[i] = sum;
    }
    //return covnumbers;
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

