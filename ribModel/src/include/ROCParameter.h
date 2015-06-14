#ifndef MODELPARAMETER_H
#define MODELPARAMETER_H

#include <vector>
#include <random>
#include <string>
#include <iostream>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

#include "Genome.h"
#include "CovarianceMatrix.h"
#include "ROCTrace.h"


class ROCParameter
{
	private:

		ROCTrace traces;
		//members
		unsigned int numParam;

		double Sphi;
		double Aphi;
		double Sphi_proposed;
		double Aphi_proposed;
		unsigned numAcceptForSphi;
		std::vector<CovarianceMatrix> covarianceMatrix;


		double phiEpsilon;
		double phiEpsilon_proposed;

		// proposal bias and std for phi values
		double bias_sphi;
		double std_sphi;
		double prev_std_sphi;

		// proposal bias and std for phi values
		double bias_phi;
		std::vector<std::vector<double>> std_phi;
		std::vector<std::vector<double>> prev_std_phi;

		// proposal bias and std for codon specific parameter
		double bias_csp;
		std::vector<double> std_csp;
		std::vector<double> prev_std_csp;

		double priorA;
		double priorB;

		std::vector<std::vector<double>> currentExpressionLevel;
		std::vector<std::vector<double>> proposedExpressionLevel;
		std::vector<std::vector<unsigned>> numAcceptForExpression;

		std::vector<std::vector<double>> currentMutationParameter;
		std::vector<std::vector<double>> proposedMutationParameter;

		std::vector<std::vector<double>> currentSelectionParameter;
		std::vector<std::vector<double>> proposedSelectionParameter;
		std::vector<unsigned> numAcceptForMutationAndSelection;

		std::string mutationSelectionState;
		unsigned numMixtures;
		unsigned numMutationCategories;
		unsigned numSelectionCategories;
		std::vector<thetaK> categories;
		std::vector<unsigned> mixtureAssignment;
		std::vector<std::vector<unsigned>> selectionIsInMixture;
		std::vector<std::vector<unsigned>> mutationIsInMixture;

		std::vector<double> categoryProbabilities;
		//static members



		// functions
		std::vector<double> propose(std::vector<double> currentParam, double (*proposal)(double a, double b), double A, std::vector<double> B);

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
		//Keywords
		static const std::string allUnique;
		static const std::string selectionShared;
		static const std::string mutationShared;

		static std::default_random_engine generator; // static to make sure that the same generator is during the runtime.

		bool checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound);
		explicit ROCParameter();
		ROCParameter(double sphi, unsigned _numMixtures,
				std::vector<unsigned> geneAssignment, std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer = true, std::string _mutationSelectionState = "allUnique");
		virtual ~ROCParameter();
		#ifndef STANDALONE
		ROCParameter(double sphi, std::vector<unsigned> geneAssignment, std::vector<unsigned> _matrix, bool splitSer = true);
		ROCParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, bool splitSer = true, std::string _mutationSelectionState = "allUnique");
		#endif
		ROCParameter(const ROCParameter& other);
		ROCParameter& operator=(const ROCParameter& rhs);
		void initParameterSet(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, 
				std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer = true, std::string _mutationSelectionState = allUnique);


		void initSelection(std::vector<double> selectionValues, unsigned mixtureElement, char aa);
		void initMutation(std::vector<double> mutationValues, unsigned mixtureElement, char aa);
		std::vector<std::vector<double>> calculateSelectionCoefficients(unsigned sample, unsigned mixture);
		#ifndef STANDALONE
		void initCovarianceMatrix(SEXP matrix, char aa);
		SEXP calculateSelectionCoefficientsR(unsigned sample, unsigned mixture);
		#endif
		CovarianceMatrix& getCovarianceMatrixForAA(char aa);
		void initAllTraces(unsigned samples, unsigned num_genes, unsigned adaptiveSamples) {traces.initAllTraces(samples, num_genes, adaptiveSamples, 
				numMutationCategories, numSelectionCategories, numParam, numMixtures, categories);}
		ROCTrace& getTraceObject() {return traces;}

		std::vector <double> readPhiValues(std::string filename);
		void setNumMutationSelectionValues(std::string mutationSelectionState, std::vector<std::vector<unsigned>> thetaKMatrix);
		void initCategoryDefinitions(std::string mutationSelectionState, std::vector<std::vector<unsigned>> thetaKMatrix);
		void initMutationSelectionCategories(std::vector<std::string> files, unsigned numCategories, unsigned paramType);
		void printThetaKMatrix()
		{
			for (unsigned i = 0u; i < numMixtures; i++)
			{
				std::cout << categories[i].delM <<"\t" << categories[i].delEta <<"\n";
			}
		}

		unsigned getNumMixtureElements() {return numMixtures;}
		unsigned getNumMutationCategories() {return numMutationCategories;}
		unsigned getNumSelectionCategories() {return numSelectionCategories;}
		unsigned getNumExpressionCategories() {return numSelectionCategories;}
		std::vector<unsigned> getMixtureElementsOfMutationCategory(unsigned category) {return mutationIsInMixture[category];}
		std::vector<unsigned> getMixtureElementsOfSelectionCategory(unsigned category) {return selectionIsInMixture[category];}
		std::string getMutationSelectionState() {return mutationSelectionState;}

		unsigned int getNumParam() {return numParam;}
		std::vector<std::vector<double>> getCurrentMutationParameter() {return currentMutationParameter;}
		std::vector<std::vector<double>> getCurrentSelectionParameter() {return currentSelectionParameter;}
		unsigned getMutationCategory(unsigned mixtureElement) {return categories[mixtureElement].delM;}
		unsigned getSelectionCategory(unsigned mixtureElement) {return categories[mixtureElement].delEta;}
		unsigned getExpressionCategory(unsigned mixtureElement) {return categories[mixtureElement].delEta;}
		double getCategoryProbability(unsigned mixtureElement) {return categoryProbabilities[mixtureElement];}
		void setCategoryProbability(unsigned mixtureElement, double value) {categoryProbabilities[mixtureElement]= value;}
		void setMixtureAssignment(unsigned gene, unsigned value) {mixtureAssignment[gene] = value;}
		unsigned getMixtureAssignment(unsigned gene) {return mixtureAssignment[gene];}

		// Phi epsilon functions
		double getPhiEpsilon() {return phiEpsilon;}

		// functions to manage adaptive step
		void adaptSphiProposalWidth(unsigned adaptationWidth);
		void adaptExpressionProposalWidth(unsigned adaptationWidth);
		void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth);

		// functions to manage Expression
		double getExpression(unsigned geneIndex, unsigned mixtureElement, bool proposed = false)
		{
			unsigned category = getSelectionCategory(mixtureElement);
			return (proposed ? proposedExpressionLevel[category][geneIndex] : currentExpressionLevel[category][geneIndex]);
		}
		void setExpression(double phi, unsigned geneIndex, unsigned mixtureElement) 
		{
			unsigned category = getSelectionCategory(mixtureElement);
			currentExpressionLevel[category][geneIndex] = phi;
		}
		double getExpressionProposalWidth(unsigned geneIndex, unsigned mixtureElement) 
		{
			unsigned category = getSelectionCategory(mixtureElement);
			return std_phi[category][geneIndex];
		}
		void updateExpression(unsigned geneIndex)
		{
			for(unsigned category = 0; category < numSelectionCategories; category++)
			{
				numAcceptForExpression[category][geneIndex]++;
				currentExpressionLevel[category][geneIndex] = proposedExpressionLevel[category][geneIndex];
			}
		}
		void updateExpression(unsigned geneIndex, unsigned mixtureElement)
		{
			unsigned category = getSelectionCategory(mixtureElement);
			numAcceptForExpression[category][geneIndex]++;
			currentExpressionLevel[category][geneIndex] = proposedExpressionLevel[category][geneIndex];
		}
		void InitializeExpression(Genome& genome, double sd_phi);
		void InitializeExpression(double sd_phi);
		void InitializeExpression(std::vector<double> expression);

		// functions to manage codon specific parameter
		void updateCodonSpecificParameter(char aa);
		double getCodonSpecificProposalWidth(unsigned aa)
		{
			unsigned codonRange[2];
			SequenceSummary::AAindexToCodonRange(aa, true, codonRange);
			return std_csp[codonRange[0]];
		}

		// functions to manage Sphi
		double getSphi(bool proposed = false) {return (proposed ? Sphi_proposed : Sphi);}
		void setSphi(double sPhi) {Sphi = sPhi;}
		void updateSphi() {Sphi = Sphi_proposed; numAcceptForSphi++;}
		double getSphiProposalWidth() {return std_sphi;}


		//update trace functions
		void updateSphiTrace(unsigned sample) {traces.updateSphiTrace(sample, Sphi);}
		void updateExpressionTrace(unsigned sample, unsigned geneIndex){traces.updateExpressionTrace(sample, geneIndex, currentExpressionLevel);}
		void updateCodonSpecificParameterTrace(unsigned sample, char aa) {traces.updateCodonSpecificParameterTrace(sample, aa, currentMutationParameter, 
				currentSelectionParameter);}
		void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex) {traces.updateMixtureAssignmentTrace(sample, geneIndex, mixtureAssignment[geneIndex]);} 
		void updateMixtureProbabilitiesTrace(unsigned samples) {traces.updateMixtureProbabilitiesTrace(samples, categoryProbabilities);}

		// poposal functions
		void proposeSPhi();
		void proposeExpressionLevels();
		void proposeCodonSpecificParameter();
		double getCurrentExpressionProposalWidth(unsigned expressionCategory, unsigned geneIndex) {return std_phi[expressionCategory][geneIndex];}
		double getPreviousExpressionProposalWidth(unsigned expressionCategory, unsigned geneIndex) {return prev_std_phi[expressionCategory][geneIndex];}
		double getCurrentSphiProposalWidth() {return std_sphi;}
		double getPreviousSphiProposalWidth() {return prev_std_sphi;}

		void getParameterForCategory(unsigned category, unsigned parameter, char aa, bool proposal, double* returnValue);


		// functions to return estimates
		double getExpressionPosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement);
		double getSphiPosteriorMean(unsigned samples);
		//double getMixtureAssignmentPosteriorMean(unsigned samples, unsigned geneIndex); // TODO: implement variance function, fix Mean function (won't work with 3 groups)
		unsigned getEstimatedMixtureAssignment(unsigned samples, unsigned geneIndex);
		std::vector <double> getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex);
		double getMutationPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);
		double getSelectionPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon);
		double getMutationVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);
		double getSelectionVariance(unsigned mixtureElement, unsigned samples, std::string &codon, bool unbiased = true);
		double getSphiVariance(unsigned samples, bool unbiased = true);
		double getExpressionVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement, bool unbiased = true);

		// static functions
		static double calculateSCUO(Gene& gene);

		static void drawIidRandomVector(unsigned draws, double mean, double sd, double (*proposal)(double a, double b), double* randomNumbers);
		static void drawIidRandomVector(unsigned draws, double r, double (*proposal)(double r), double* randomNumber);
		static double randNorm(double mean, double sd);
		static double randLogNorm(double m, double s);
		static double randExp(double r);
		static void randDirichlet(double* input, unsigned numElements, double* output);
		static double randUnif(double minVal, double maxVal);
		static unsigned randMultinom(double* probabilities, unsigned mixtureElements);

		static double densityNorm(double x, double mean, double sd);
		static double densityLogNorm(double x, double mean, double sd);


		//R wrapper functions
		void initMutationSelectionCategoriesR(std::vector<std::string> files, unsigned numCategories, std::string paramType);
		double getMutationPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon)
		{
			double rv = -1.0;
			bool check = checkIndex(mixtureElement, 1, numMixtures);
			if (check)
			{ 
				rv = getMutationPosteriorMean(mixtureElement - 1, samples, codon);
			}
			return rv;
		}
		double getSelectionPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon)
		{
			double rv = -1.0;
			bool check = checkIndex(mixtureElement, 1, numMixtures);
			if (check)
			{
				rv = getSelectionPosteriorMean(mixtureElement - 1, samples, codon);
			}
			return rv;
		}
		double getMutationVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased)
		{
			double rv = -1.0;
			bool check = checkIndex(mixtureElement, 1, numMixtures);
			if (check)
			{
				rv = getMutationVariance(mixtureElement - 1, samples, codon, unbiased);
			}
			return rv;
		}
		double getSelectionVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased)
		{
			double rv = -1.0;
			bool check = checkIndex(mixtureElement, 1, numMixtures);
			if (check)
			{
				rv = getSelectionVariance(mixtureElement - 1, samples, codon, unbiased);
			}
			return rv;
		}

		void initializeExpressionByGenome(Genome& genome, double sd_phi) {InitializeExpression(genome, sd_phi);}
		void initializeExpressionByList(double sd_phi) {InitializeExpression(sd_phi);}
		void initializeExpressionByRandom(std::vector<double> expression) {InitializeExpression(expression);}

		unsigned getEstimatedMixtureAssignmentForGene(unsigned samples, unsigned geneIndex)
		{
			bool check = checkIndex(geneIndex, 1, mixtureAssignment.size());
			return check ? getEstimatedMixtureAssignment(samples, geneIndex - 1) + 1 : 0;
		}

		std::vector<double> getEstimatedMixtureAssignmentProbabilitiesForGene(unsigned samples, unsigned geneIndex)
		{
			std::vector <double> probabilities;
			bool check = checkIndex(geneIndex, 1, mixtureAssignment.size());
			if (check) 
			{
				probabilities = getEstimatedMixtureAssignmentProbabilities(samples, geneIndex - 1);
			}
			return probabilities;
		}


		void setMixtureAssignmentForGene(unsigned geneIndex, unsigned value);

		double getExpressionPosteriorMeanByMixtureElementForGene(unsigned samples, unsigned geneIndex, unsigned mixtureElement)
		{
			double rv = -1.0;
			bool checkGene = checkIndex(geneIndex, 1, mixtureAssignment.size());
			bool checkMixtureElement = checkIndex(mixtureElement, 1, numMixtures); 
			if (checkGene && checkMixtureElement)
			{
				rv = getExpressionPosteriorMean(samples, geneIndex - 1, mixtureElement - 1);
			}
			return rv;
		}
		double getExpressionVarianceByMixtureElementForGene(unsigned samples, unsigned geneIndex, unsigned mixtureElement, bool unbiased)
		{
			double rv = -1.0;
			bool checkGene = checkIndex(geneIndex, 1, mixtureAssignment.size());
			bool checkMixtureElement = checkIndex(mixtureElement, 1, numMixtures);
			if (checkGene && checkMixtureElement)
			{
				rv = getExpressionVariance(samples, geneIndex - 1, mixtureElement - 1, unbiased);
			}
			return rv;
		}
		unsigned getMutationCategoryForMixture(unsigned mixtureElement);
		unsigned getSelectionCategoryForMixture(unsigned mixtureElement);
		unsigned getExpressionCategoryForMixture(unsigned mixtureElement);
		std::vector<double> getCurrentExpressionForMixture(unsigned mixture)
		{
			bool checkMixture = checkIndex(mixture, 1, numMixtures);
			unsigned exprCat = 0u;
			if(checkMixture)
			{
				exprCat = getExpressionCategory(mixture - 1);
			}else{
				std::cerr << "WARNING: Mixture element " << mixture << " NOT found. Mixture element 1 is returned instead. \n";
			}
			return currentExpressionLevel[exprCat];
		}
	protected:

};

#endif // MODELPARAMETER_H
