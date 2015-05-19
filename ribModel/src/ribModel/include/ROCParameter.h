#ifndef MODELPARAMETER_H
#define MODELPARAMETER_H

#include <vector>
#include <random>
#include <string>
#include <iostream>

#include "../include/SequenceSummary.h"
#include "../include/Genome.h"

class ROCParameter
{
    private:

        struct thetaK
        {
            unsigned delM;
            unsigned delEta;
        };
        //members
        unsigned int numParam;
				
        double Sphi;
        double Aphi;
        double Sphi_proposed;
        double Aphi_proposed;
        std::vector<double> sPhiTrace;
        std::vector<double> aPhiTrace;
        unsigned numAcceptForSphi;

        double phiEpsilon;
        double phiEpsilon_proposed;

				//Keywords
				static const std::string allUnique;
				static const std::string selectionShared;
				static const std::string mutationShared;

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
        std::vector<unsigned> numAcceptForMutationAndSelection;

        unsigned numMixtures;
				unsigned numMutationCategories;
				unsigned numSelectionCategories;
        std::vector<thetaK> categories;
        std::vector<unsigned> mixtureAssignment;

        std::vector<std::vector<unsigned>> mixtureAssignmentTrace;
        std::vector<double> categoryProbabilities;
        std::vector<std::vector<double>> categoryProbabilitiesTrace;
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
				ROCParameter(unsigned numGenes, double sphi, unsigned _numMixtures, double* geneAssignment = nullptr, bool splitSer = true, std::string mutationSelectionState = "allUnique", unsigned thetaKMatrix[][2] = nullptr, std::string files[] = nullptr);
        virtual ~ROCParameter();
        ROCParameter(const ROCParameter& other);
        ROCParameter& operator=(const ROCParameter& rhs);

        void initCategoryDefinitions(std::string mutationSelectionState, unsigned thetaKMatrix[][2]);
				void initMutationSelectionCategories(std::string files[], int numCategories, unsigned paramType);
				void printThetaKMatrix()
				{
					for (int i = 0; i < numMixtures; i++)
					{
						std::cout << categories[i].delM <<"\t" << categories[i].delEta <<"\n";
					}
				}

        unsigned getNumMixtureElements() {return numMixtures;}
        unsigned getNumMutationCategoriesVar() {return numMutationCategories;}
				unsigned getNumSelectionCategoriesVar() {return numSelectionCategories;}
				unsigned getMutationCategory(unsigned group) {return categories[group].delM;}
        unsigned getSelectionCategory(unsigned group) {return categories[group].delEta;}
        unsigned getExpressionCategory(unsigned group) {return categories[group].delEta;}
        double getCategoryProbability(unsigned group) {return categoryProbabilities[group];}
        void setCategoryProbability(unsigned group, double value) {categoryProbabilities[group]= value;}
        void setMixtureAssignment(unsigned gene, unsigned value) {mixtureAssignment[gene] = value;}
        unsigned getMixtureAssignment(unsigned gene) {return mixtureAssignment[gene];}
        void updateCategoryProbabilitiesTrace(int samples)
        {
            for(unsigned category = 0; category < numMixtures; category++)
            {
                categoryProbabilitiesTrace[category][samples] = categoryProbabilities[category];
            }
        }

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
            numAcceptForExpression[geneIndex]++;
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
        void updateCodonSpecificParameter(char aa);

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
        void initMixtureAssignmentTrace(unsigned samples, unsigned num_genes);
        void initCategoryProbabilitesTrace(int samples);
        void updateExpressionTrace(unsigned sample, unsigned geneIndex)
        {
            for(unsigned category = 0; category < numMixtures; category++)
            {
                expressionTrace[category][sample][geneIndex] = currentExpressionLevel[category][geneIndex];
            }
        }
        void updateSphiTrace(unsigned sample) {sPhiTrace[sample] = Sphi;}
        void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex) {mixtureAssignmentTrace[sample][geneIndex] = mixtureAssignment[geneIndex];}
        std::vector<double> getAPhiTrace() {return aPhiTrace;}
        std::vector<std::vector <unsigned>> getMixtureAssignmentTrace() {return mixtureAssignmentTrace;}

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
        static unsigned randMultinom(double* probabilities, unsigned groups);

        static double densityNorm(double x, double mean, double sd);
        static double densityLogNorm(double x, double mean, double sd);

    protected:

};

#endif // MODELPARAMETER_H
