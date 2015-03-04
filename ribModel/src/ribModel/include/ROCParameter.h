#ifndef MODELPARAMETER_H
#define MODELPARAMETER_H

#include <vector>
#include <random>
#include <string>


#include "../include/SequenceSummary.h"
#include "../include/Genome.h"

class ROCParameter
{
    private:
        //members
        unsigned numParam = 40;

        double Sphi;
        double Aphi;
        double Sphi_proposed;
        double Aphi_proposed;
        std::vector<double> sPhiTrace;
        std::vector<double> aPhiTrace;
        // proposal bias and std for phi values
        double bias_sphi;
        double std_sphi;

        // proposal bias and std for phi values
        double bias_phi;
        double std_phi;

        // proposal bias and std for codon specific parameter
        double bias_csp;
        double std_csp;

        double priorA;
        double priorB;

        std::vector<double> currentExpressionLevel;
        std::vector<double> proposedExpressionLevel;
        std::vector<std::vector<double>> expressionTrace;

        std::vector<std::vector<double>> currentMutationParameter;
        std::vector<std::vector<double>> proposedMutationParameter;

        std::vector<std::vector<double>> currentSelectionParameter;
        std::vector<std::vector<double>> proposedSelectionParameter;

        //static members



        // functions
        std::vector<double> propose(std::vector<double> currentParam, double (*proposal)(double a, double b), double A, double B);
        double calculateSCUO(Gene& gene);

        static void quickSortPair(double a[], int b[], int first, int last);
        static void quickSort(double a[], int first, int last);
        static int pivotPair(double a[], int b[], int first, int last);
        static int pivot(double a[], int first, int last);
        static void swap(double& a, double& b);
        static void swap(int& a, int& b);

    public:
        //static const members
        static const unsigned dM;
        static const unsigned dEta;
        static std::default_random_engine generator; // static to make sure that the same generator is during the runtime.

        explicit ROCParameter();
        ROCParameter(unsigned numMutationCategories, unsigned numSelectionCategories, unsigned numGenes, double sphi);
        virtual ~ROCParameter();
        //ROCParameter(const ROCParameter& other);

        double getExpression(int geneIndex, bool proposed = false) {return (proposed ? proposedExpressionLevel[geneIndex] : currentExpressionLevel[geneIndex]);}
        void setExpression(double phi, int geneIndex) {currentExpressionLevel[geneIndex] = phi;}
        void updateExpression(int geneIndex) {currentExpressionLevel[geneIndex] = proposedExpressionLevel[geneIndex];}

        double getSphi(bool proposed = false) {return (proposed ? Sphi_proposed : Sphi);}
        double setSphi(double sPhi) {Sphi = sPhi;}
        double updateSphi() {Sphi = Sphi_proposed;}
        std::vector<double> getSPhiTrace() {return sPhiTrace;}

        void proposeSPhi();
        void proposeExpressionLevels();
        void proposeCodonSpecificParameter();

        void InitializeExpression(Genome& genome, double sd_phi);
        void InitializeExpression(double sd_phi);
        void InitializeExpression(double* expression);

        void getParameterForCategory(unsigned category, unsigned parameter, char aa, bool proposal, double* returnValue);

        void initAllTraces(int samples, int num_genes);
        void initExpressionTrace(int samples, int num_genes);
        void initSphiTrace(int sample) {sPhiTrace.resize(sample);}
        void updateExpressionTrace(int sample, int geneIndex) {expressionTrace[sample][geneIndex] = currentExpressionLevel[geneIndex];}
        void updateSphiTrace(int sample) {sPhiTrace[sample] = Sphi;}

        double getExpressionPosteriorMean(int samples, int geneIndex);
        double getSphiPosteriorMean(int samples);


        //static functions
        static double randNorm(double mean, double sd);
        static double randLogNorm(double m, double s);
        static double randExp(double r);

        static double densityNorm(double x, double mean, double sd);
        static double densityLogNorm(double x, double mean, double sd);

    protected:

};

#endif // MODELPARAMETER_H
