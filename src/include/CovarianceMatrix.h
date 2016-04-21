#ifndef COVARIANCEMATRIX_H
#define COVARIANCEMATRIX_H

#include <vector>
#include <iostream>
#include <cmath>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

class CovarianceMatrix
{
    private:
        std::vector<double> covMatrix;
        std::vector<double> choleskiMatrix;
        int numVariates; //make static const again

		double sampleMean(std::vector<double> sampleVector, unsigned samples, unsigned lastIteration);

    public:
        //Constructors & Destructors:
        CovarianceMatrix();
        CovarianceMatrix(int _numVariates);
        CovarianceMatrix(std::vector <double> &matrix);
        CovarianceMatrix(const CovarianceMatrix& other);
        CovarianceMatrix& operator=(const CovarianceMatrix& other);
		CovarianceMatrix& operator+(const CovarianceMatrix& rhs);
		CovarianceMatrix& operator*(const double &value);
        void operator*=(const double &value);
        virtual ~CovarianceMatrix();



        //Matrix Functions:
	    void initCovarianceMatrix(unsigned _numVariates);
		void setDiag(double val);
        void choleskiDecomposition();
	    void printCovarianceMatrix();
        void printCholeskiMatrix();
        std::vector<double>* getCovMatrix();
        int getNumVariates();
        std::vector<double> transformIidNumersIntoCovaryingNumbers(std::vector<double> iidnumbers);
		void calculateSampleCovariance(std::vector<std::vector<std::vector<std::vector<double>>>> codonSpecificParameterTrace, std::string aa, unsigned samples, unsigned lastIteration);

#ifndef STANDALONE
    void setCovarianceMatrix(SEXP _matrix);
#endif


    protected:
};

#endif // COVARIANCEMATRIX_H
