#ifndef COVARIANCEMATRIX_H
#define COVARIANCEMATRIX_H

#include <vector>

class CovarianceMatrix
{
    private:
    std::vector<double> covMatrix;
    std::vector<double> choleskiMatrix;
    int numVariates; //make static const again


    public:
        CovarianceMatrix();
        CovarianceMatrix(int _numVariates);
				CovarianceMatrix(std::vector <double> &matrix);
        virtual ~CovarianceMatrix();
        CovarianceMatrix(const CovarianceMatrix& other);
        //CovarianceMatrix& operator=(const CovarianceMatrix& other);
        void choleskiDecomposition();
        void calculateCovarianceMatrixFromTraces(std::vector <std::vector <double>> trace);
				void printCovarianceMatrix();
        void printCholeskiMatrix();
        void transformIidNumersIntoCovaryingNumbers(double* iidnumbers, double* covnumbers);

    protected:
};

#endif // COVARIANCEMATRIX_H
