#ifndef MYCOVARIANCEMATRIX_H
#define MYCOVARIANCEMATRIX_H

#include <vector>

class CovarianceMatrix
{
    private:
    std::vector<std::vector<double>> covMatrix;
    std::vector<std::vector<double>> choleskiMatrix;
    int numVariates; //make static const again


    public:
        CovarianceMatrix();
        CovarianceMatrix(int _numVariates);
				CovarianceMatrix(std::vector <std::vector<double>> &matrix);
        virtual ~CovarianceMatrix();
        CovarianceMatrix(const CovarianceMatrix& other);
        //CovarianceMatrix& operator=(const CovarianceMatrix& other);
        void choleskiDecomposition();
        void printCovarianceMatrix();
        void printCholeskiMatrix();
        void transformIidNumersIntoCovaryingNumbers(double* iidnumbers, double* covnumbers);

    protected:
};

#endif // MYCOVARIANCEMATRIX_H
