#include "CovarianceMatrix.h"

#include <iostream>
#include <cmath>

/*
    !!!!IMPORTANT NOTE AT THE BOTTOM OF FILE!!!!
*/

template<unsigned _numVariates>
CovarianceMatrix<_numVariates>::CovarianceMatrix()
{
    for(int i = 0; i < numVariates; i++)
    {
        for(int j = 0; j < numVariates; j++)
        {
            covMatrix[i][j] = (i == j ? 1.0 : 0.0);
            choleskiMatrix[i][j] = 0.0;
        }
    }
}

template<unsigned _numVariates>
CovarianceMatrix<_numVariates>::CovarianceMatrix(double (&matrix)[_numVariates][_numVariates])
{
    for(int i = 0; i < numVariates; i++)
    {
        for(int j = 0; j < numVariates; j++)
        {
            covMatrix[i][j] = matrix[i][j];
            choleskiMatrix[i][j] = 0.0;
        }
    }
}

template<unsigned _numVariates>
CovarianceMatrix<_numVariates>::~CovarianceMatrix()
{
    //dtor
}

template<unsigned _numVariates>
CovarianceMatrix<_numVariates>::CovarianceMatrix(const CovarianceMatrix& other)
{
    covMatrix[other.numVariates][other.numVariates];
    numVariates = other.numVariates;
    for(int i = 0; i < other.numVariates; i++)
    {
        for(int j = 0; j < other.numVariates; j++)
        {
            covMatrix[i][j] = other.covMatrix[i][j];
        }
    }
}

template<unsigned L> template<unsigned R>
CovarianceMatrix<L>& CovarianceMatrix<L>::operator=(const CovarianceMatrix<R>& rhs)
{
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

// addaptatoin of http://en.wikipedia.org/wiki/Cholesky_decomposition
// http://rosettacode.org/wiki/Cholesky_decomposition#C
template<unsigned _numVariates>
void CovarianceMatrix<_numVariates>::choleskiDecomposition()
{
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

template<unsigned _numVariates>
void CovarianceMatrix<_numVariates>::transformIidNumersIntoCovaryingNumbers(double* iidnumbers, double* covnumbers)
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

template<unsigned _numVariates>
void CovarianceMatrix<_numVariates>::printCovarianceMatrix()
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
template<unsigned _numVariates>
void CovarianceMatrix<_numVariates>::printCholeskiMatrix()
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

// To avoid this you can move all the code in this CPP file (exluding the next few template lines) into the header file
// I do this to keep the header file readable and since I know I only need covariance matrices of size 2,3,4,5 for the codons
// http://stackoverflow.com/questions/8752837/undefined-reference-to-template-class-constructor
template class CovarianceMatrix<2>;
template class CovarianceMatrix<3>;
template class CovarianceMatrix<4>;
template class CovarianceMatrix<5>;


