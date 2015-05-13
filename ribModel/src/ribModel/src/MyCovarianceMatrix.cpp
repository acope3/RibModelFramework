#include "../include/MyCovarianceMatrix.h"

#include <iostream>
#include <cmath>


CovarianceMatrix::CovarianceMatrix()
{
		std::cout <<"0 arg constructor called\n";
		numVariates = 2;
		covMatrix.resize(numVariates);
		choleskiMatrix.resize(numVariates);

    for(int i = 0; i < numVariates; i++)
    {
				covMatrix[i].resize(numVariates);
				choleskiMatrix[i].resize(numVariates);
        for(int j = 0; j < numVariates; j++)
        {
            covMatrix[i][j] = (i == j ? 1.0 : 0.0);
            choleskiMatrix[i][j] = 0.0;
				}
		}
}

CovarianceMatrix::CovarianceMatrix(int _numVariates)
{

		std::cout <<"1 arg constructor called\n";
		numVariates = _numVariates;
		covMatrix.resize(numVariates);
		choleskiMatrix.resize(numVariates);

    for(int i = 0; i < numVariates; i++)
    {
				covMatrix[i].resize(numVariates);
				choleskiMatrix[i].resize(numVariates);
        for(int j = 0; j < numVariates; j++)
        {
            covMatrix[i][j] = (i == j ? 1.0 : 0.0);
            choleskiMatrix[i][j] = 0.0;
        }
    }
}

CovarianceMatrix::CovarianceMatrix(std::vector <std::vector <double>> &matrix)
{
		std::cout <<"Matrix constructor called\n";
		numVariates = matrix.size();   
		covMatrix.resize(numVariates);
		choleskiMatrix.resize(numVariates);
		
		for(int i = 0; i < numVariates; i++)
    {
				covMatrix[i].resize(numVariates);
				choleskiMatrix[i].resize(numVariates);
        for(int j = 0; j < numVariates; j++)
        {
            covMatrix[i][j] = matrix[i][j];
            choleskiMatrix[i][j] = 0.0;
        }
    }
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
            double a = covMatrix[i][j];
            double LsubstractSum = 0.0;
            for(int k = 0; k < j; k++)
            {
                double b = covMatrix[i][k];
                double c = covMatrix[j][k];
                LsubstractSum += choleskiMatrix[i][k] * choleskiMatrix[j][k];
            }
            choleskiMatrix[i][j] = (i == j) ? std::sqrt(covMatrix[i][j] - LsubstractSum) :
                (1.0 / choleskiMatrix[j][j]) * (covMatrix[i][j] - LsubstractSum);
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
            sum += choleskiMatrix[i][k] * iidnumbers[k];
        }

        covnumbers[i] = sum;
    }
    //return covnumbers;
}

void CovarianceMatrix::printCovarianceMatrix()
{
    for(int i = 0; i < numVariates; i++)
    {
        for(int j = 0; j < numVariates; j++)
        {
            std::cout << covMatrix[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}
void CovarianceMatrix::printCholeskiMatrix()
{
    for(int i = 0; i < numVariates; i++)
    {
        for(int j = 0; j < numVariates; j++)
        {
            std::cout << choleskiMatrix[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

