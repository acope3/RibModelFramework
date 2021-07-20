#ifndef COVARIANCEMATRIX_H
#define COVARIANCEMATRIX_H


#include <vector>
#include <cmath>
#include <string>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

class CovarianceMatrix
{
    private:
        std::vector<double> covMatrix; //sample covariance matrix
        std::vector<double> choleskyMatrix; //Cholesky factorization of proposal covariance matrix
        unsigned numVariates; //make static const again
        double sampleMean(std::vector<float> sampleVector, unsigned samples, unsigned latestSample,bool log_scale=false);
                double lambda = 1; // Scales proposal matrix
                // Can be fixed or dynamic
                double gamma = 0.6 // Scales adjustment of lambda; \in (1/2, 1) but 1 is not advised by Vihola2012, p. 999 remark 3
                // From equation 29 in LuengoEtAl2020
                // log(lambda_t) = log(lambda_{t-1}) + gamma_t(alpha_t - alphaStar)
    public:
        //Constructors & Destructors:
        CovarianceMatrix(); //Defaults to initCovarianceMatrix(2)
        CovarianceMatrix(unsigned _numVariates);
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
        std::vector<double> transformIidNumbersIntoCovaryingNumbers(std::vector<double> iidNumbers);
	void calculateSampleCovariance(std::vector<std::vector<std::vector<std::vector<float>>>> codonSpecificParameterTrace, std::string aa, unsigned samples, unsigned latestSample);
        void calculateSampleCovarianceForPANSE(std::vector<std::vector<std::vector<std::vector<float>>>> codonSpecificParameterTrace, std::string codon, unsigned samples, unsigned latestSample);
	void updateSampleCovariance(std::vector<std::vector<std::vector<std::vector<float>>>> codonSpecificParameterTrace, std::string aa, unsigned samples, unsigned latestSample);
    	void updateCholeskyCovariance(std::vector<std::vector<std::vector<std::vector<float>>>> codonSpecificParameterTrace, std::string aa, unsigned samples, unsigned latestSample);

#ifndef STANDALONE
    void setCovarianceMatrix(SEXP _matrix);
#endif


    protected:
};

#endif // COVARIANCEMATRIX_H
