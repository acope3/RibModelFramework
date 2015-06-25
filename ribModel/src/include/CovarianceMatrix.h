#ifndef COVARIANCEMATRIX_H
#define COVARIANCEMATRIX_H

#include <vector>
#ifndef STANDALONE
#include <Rcpp.h>
#endif
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
        CovarianceMatrix& operator=(const CovarianceMatrix& other);
        void operator*(const double &value);
        void operator*=(const double &value);
				void initCovarianceMatrix(unsigned _numVariates);
        void choleskiDecomposition();
        void calculateCovarianceMatrixFromTraces(std::vector<std::vector <std::vector<double>>> trace, unsigned geneIndex, unsigned curSample, unsigned adaptiveWidte);
				void printCovarianceMatrix();
        void printCholeskiMatrix();
        std::vector<double> transformIidNumersIntoCovaryingNumbers(std::vector<double> iidnumbers);
				std::vector<double>* getCovMatrix();
				int getNumVariates() {return numVariates;}

				#ifndef STANDALONE
				void setCovarianceMatrix(SEXP _matrix);
				#endif
    protected:
};

#endif // COVARIANCEMATRIX_H
