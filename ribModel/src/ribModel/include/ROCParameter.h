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
        unsigned numAcceptForSphi;

        double phiEpsilon;
        double phiEpsilon_proposed;

        // proposal bias and std for phi values
        double bias_sphi;
        double std_sphi;

        // proposal bias and std for phi values
        double bias_phi;
        std::vector<double> std_phi;

        // proposal bias and std for codon specific parameter
        double bias_csp;
        double std_csp;

        double priorA;
        double priorB;

        std::vector<std::vector<double>> currentExpressionLevel;
        std::vector<std::vector<double>> proposedExpressionLevel;
        std::vector<std::vector<std::vector<double>>> expressionTrace;
        std::vector<unsigned> numAcceptForExpression;

        std::vector<std::vector<double>> currentMutationParameter;
        std::vector<std::vector<double>> proposedMutationParameter;

        std::vector<std::vector<double>> currentSelectionParameter;
        std::vector<std::vector<double>> proposedSelectionParameter;


        unsigned numMixtures;
        std::vector<unsigned[2]> categories;
        std::vector<std::vector<double>> categoryProbabilities;
        //static members



        // functions
        std::vector<double> propose(std::vector<double> currentParam, double (*proposal)(double a, double b), double A, double B);


        // sorting functions
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
        ROCParameter(unsigned numMutationCategories, unsigned numSelectionCategories, unsigned numGenes, double sphi, unsigned _numMixtures);
        virtual ~ROCParameter();
        ROCParameter(const ROCParameter& other);
        ROCParameter& operator=(const ROCParameter& rhs);

        unsigned getNumMixtureElements() {return numMixtures;}
        unsigned getMutationCategory(unsigned group) {return categories[group][0];}
        unsigned getSelectionCategory(unsigned group) {return categories[group][1];}
        unsigned getExpressionCategory(unsigned group) {return categories[group][1];}
        double getCategoryProbability(unsigned group, unsigned gene) {return categoryProbabilities[group][gene];}
        void setCategoryProbability(unsigned group, unsigned gene, double value) {categoryProbabilities[group][gene] = value;}

        // Phi epsilon functions
        double getPhiEpsilon() {return phiEpsilon;}

        // functions to manage adaptive step
        void adaptSphiProposalWidth(unsigned adaptationWidth);
        void adaptExpressionProposalWidth(unsigned adaptationWidth);
        void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth);

        // functions to manage Expression
        double getExpression(unsigned geneIndex, unsigned category, bool proposed = false)
        {
            return (proposed ? proposedExpressionLevel[category][geneIndex] : currentExpressionLevel[category][geneIndex]);
        }
        void setExpression(double phi, unsigned geneIndex, unsigned category) {currentExpressionLevel[category][geneIndex] = phi;}
        void updateExpression(unsigned geneIndex)
        {
            for(unsigned category = 0; category < numMixtures; category++)
            {
                currentExpressionLevel[category][geneIndex] = proposedExpressionLevel[category][geneIndex]; numAcceptForExpression[geneIndex]++;
            }
        }
        std::vector<std::vector<std::vector<double>>> getExpressionTrace() {return expressionTrace;}
        void InitializeExpression(Genome& genome, double sd_phi);
        void InitializeExpression(double sd_phi);
        void InitializeExpression(double* expression);

        // functions to manage codon specific parameter
        void updateCodonSpecificParameter(unsigned category, unsigned paramType, char aa);

        // functions to manage Sphi
        double getSphi(bool proposed = false) {return (proposed ? Sphi_proposed : Sphi);}
        void setSphi(double sPhi) {Sphi = sPhi;}
        void updateSphi() {Sphi = Sphi_proposed; numAcceptForSphi++;}
        std::vector<double> getSPhiTrace() {return sPhiTrace;}
        double getSphiProposalWidth() {return std_sphi;}

        // poposal functions
        void proposeSPhi();
        void proposeExpressionLevels();
        void proposeCodonSpecificParameter();

        void getParameterForCategory(unsigned category, unsigned parameter, char aa, bool proposal, double* returnValue);
        unsigned getNumMutationCategories() {return currentMutationParameter.size();}
        unsigned getNumSelectionCategories() {return currentSelectionParameter.size();}

        // functions to manage traces
        void initAllTraces(unsigned samples, unsigned num_genes);
        void initExpressionTrace(unsigned samples, unsigned num_genes);
        void initSphiTrace(unsigned sample) {sPhiTrace.resize(sample);}
        void updateExpressionTrace(unsigned sample, unsigned geneIndex)
        {
            for(unsigned category = 0; category < numMixtures; category++)
            {
                expressionTrace[category][sample][geneIndex] = currentExpressionLevel[category][geneIndex];
            }
        }
        void updateSphiTrace(unsigned sample) {sPhiTrace[sample] = Sphi;}

        // functions to return estimates
        double getExpressionPosteriorMean(unsigned samples, unsigned geneIndex, unsigned category);
        double getSphiPosteriorMean(unsigned samples);


        // static functions
        static double calculateSCUO(Gene& gene);

        static double* drawIidRandomVector(unsigned draws, double mean, double sd, double (*proposal)(double a, double b));
        static double* drawIidRandomVector(unsigned draws, double r, double (*proposal)(double r));
        static double randNorm(double mean, double sd);
        static double randLogNorm(double m, double s);
        static double randExp(double r);
        static void randDirichlet(double* input, unsigned numElements, double* output);
        static void randMultinom(double* probabilities, unsigned groups, unsigned* output);

        static double densityNorm(double x, double mean, double sd);
        static double densityLogNorm(double x, double mean, double sd);

    protected:

};

#endif // MODELPARAMETER_H
