#ifndef COVARIANCEMATRIX_H
#define COVARIANCEMATRIX_H

#include <vector>

template<unsigned _numVariates>
class CovarianceMatrix
{
    private:
    //std::vector<std::vector<double>> covMatrix;
    //std::vector<std::vector<double>> choleskiMatrix;
    double covMatrix[_numVariates][_numVariates];
    double choleskiMatrix[_numVariates][_numVariates];
    static const int numVariates = _numVariates;


    public:
        CovarianceMatrix();
        CovarianceMatrix(double (&covMatrix)[_numVariates][_numVariates]);
        virtual ~CovarianceMatrix();
        CovarianceMatrix(const CovarianceMatrix& other);
        template<unsigned R> CovarianceMatrix<_numVariates>& operator=(const CovarianceMatrix<R>& other);

        void choleskiDecomposition();
        void printCovarianceMatrix();
        void printCholeskiMatrix();
        void transformIidNumersIntoCovaryingNumbers(double* iidnumbers, double* covnumbers);

    protected:
};

#endif // COVARIANCEMATRIX_H
