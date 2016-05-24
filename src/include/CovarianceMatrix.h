#ifndef COVARIANCEMATRIX_H
#define COVARIANCEMATRIX_H


#include <vector>
#include <cmath>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

class CovarianceMatrix
{
    private:
        std::vector<double> covMatrix;
        std::vector<double> choleskyMatrix;
        int numVariates; //make static const again

		double sampleMean(std::vector<double> sampleVector, unsigned samples, unsigned lastIteration);

    public:
        //Constructors & Destructors:
        CovarianceMatrix(); //TODO: Currently defaults to initCovarianceMatrix(2)
        CovarianceMatrix(int _numVariates);
        CovarianceMatrix(std::vector <double> &matrix);
        CovarianceMatrix(const CovarianceMatrix& other);
        CovarianceMatrix& operator=(const CovarianceMatrix& other);
		CovarianceMatrix& operator+(const CovarianceMatrix& rhs);
		CovarianceMatrix& operator*(const double &value);
        void operator*=(const double &value);
		bool operator==(const CovarianceMatrix& other) const;
        virtual ~CovarianceMatrix();



        //Matrix Functions:
	    void initCovarianceMatrix(unsigned _numVariates);
		void setDiag(double val);
        void choleskyDecomposition();
	    void printCovarianceMatrix();
        void printCholeskyMatrix();
        std::vector<double>* getCovMatrix();
		std::vector<double>* getCholeskyMatrix(); //Only for unit testing.
        int getNumVariates();
        std::vector<double> transformIidNumersIntoCovaryingNumbers(std::vector<double> iidnumbers);
		void calculateSampleCovariance(std::vector<std::vector<std::vector<std::vector<double>>>> codonSpecificParameterTrace, std::string aa, unsigned samples, unsigned lastIteration);

#ifndef STANDALONE
    void setCovarianceMatrix(SEXP _matrix);
#endif


    protected:
};

#endif // COVARIANCEMATRIX_H
