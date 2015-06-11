#ifndef MODELPARAMETER_H
#define MODELPARAMETER_H

#include <vector>
#include <random>
#include <string>
#include <iostream>

#ifndef STANDALONE
#include <Rcpp.h>
#endif
//#include "../include/SequenceSummary.h"
#include "Genome.h"
#include "CovarianceMatrix.h"

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
		std::vector<CovarianceMatrix> covarianceMatrix;

		std::vector<double> sphiAcceptanceRatioTrace; // sample
		std::vector<std::vector<std::vector<double>>> expressionAcceptanceRatioTrace; // order: category, gene, sample
		std::vector<std::vector<double>> cspAcceptanceRatioTrace; // order: codon, sample

		double phiEpsilon;
		double phiEpsilon_proposed;

		// proposal bias and std for phi values
		double bias_sphi;
		double std_sphi;

		// proposal bias and std for phi values
		double bias_phi;
		std::vector<std::vector<double>> std_phi;
		std::vector<std::vector<double>> prev_std_phi;

		// proposal bias and std for codon specific parameter
		double bias_csp;
		std::vector<double> std_csp;

		double priorA;
		double priorB;

		std::vector<std::vector<double>> currentExpressionLevel;
		std::vector<std::vector<double>> proposedExpressionLevel;
		std::vector<std::vector<std::vector<double>>> expressionTrace;
		std::vector<std::vector<std::vector<double>>> mutationParameterTrace;
		std::vector<std::vector<std::vector<double>>> selectionParameterTrace;
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

		std::vector<std::vector<unsigned>> mixtureAssignmentTrace;
		std::vector<double> categoryProbabilities;
		std::vector<std::vector<double>> mixtureProbabilitiesTrace;
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
		ROCParameter(const ROCParameter& other);
		ROCParameter& operator=(const ROCParameter& rhs);
		void initParameterSet(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, 
				std::vector<std::vector<unsigned>> thetaKMatrix, bool splitSer = true, std::string _mutationSelectionState = allUnique);

		#ifndef STANDALONE
		ROCParameter(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, bool splitSer = true, std::string _mutationSelectionState = allUnique, SEXP _matrix = NULL);
		#endif

		void initSelection(std::vector<double> selectionValues, unsigned mixtureElement, char aa);
		void initMutation(std::vector<double> mutationValues, unsigned mixtureElement, char aa);
		#ifndef STANDALONE
		void initCovarianceMatrix(SEXP matrix, char aa);
		#endif
		CovarianceMatrix& getCovarianceMatrixForAA(char aa);
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
		void updateMixtureProbabilitiesTrace(int samples)
		{
			for(unsigned category = 0; category < numMixtures; category++)
			{
				mixtureProbabilitiesTrace[category][samples] = categoryProbabilities[category];
			}
		}

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
		std::vector<std::vector<std::vector<double>>> getExpressionTrace() {return expressionTrace;}
		std::vector<double> getExpressionTrace(unsigned geneIndex);
		std::vector<std::vector<std::vector<double>>> getMutationParameterTrace() {return mutationParameterTrace;}
		std::vector<std::vector<std::vector<double>>> getSelectionParameterTrace() {return selectionParameterTrace;}
		std::vector<double> getMixtureProbabilitiesTrace(unsigned mixtureIndex) {return mixtureProbabilitiesTrace[mixtureIndex];}
		std::vector<double> getSphiAcceptanceRatioTrace() {return sphiAcceptanceRatioTrace;}
		std::vector<double> getExpressionAcceptanceRatioTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex) 
		{
			unsigned category = getExpressionCategory(mixtureElement);
			return expressionAcceptanceRatioTrace[category][geneIndex];
		}
		std::vector<double> getCspAcceptanceRatioTraceForAA(char aa)
		{
			aa = std::toupper(aa);
			unsigned aaIndex = SequenceSummary:: aaToIndex.find(aa) -> second;
			return cspAcceptanceRatioTrace[aaIndex];
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
		std::vector<double> getSPhiTrace() {return sPhiTrace;}
		double getSphiProposalWidth() {return std_sphi;}

		// poposal functions
		void proposeSPhi();
		void proposeExpressionLevels();
		void proposeCodonSpecificParameter();
		double getCurrentExpressionProposalWidth(unsigned expressionCategory, unsigned geneIndex) {return std_phi[expressionCategory][geneIndex];}
		double getPreviousExpressionProposalWidth(unsigned expressionCategory, unsigned geneIndex) {return prev_std_phi[expressionCategory][geneIndex];}

		void getParameterForCategory(unsigned category, unsigned parameter, char aa, bool proposal, double* returnValue);

		// functions to manage traces
		void initAllTraces(unsigned samples, unsigned num_genes, unsigned adaptiveSamples);
		void initExpressionTrace(unsigned samples, unsigned num_genes);
		void initSelectionParameterTrace(unsigned samples);
		void initMutationParameterTrace(unsigned samples);
		void initSphiTrace(unsigned sample) {sPhiTrace.resize(sample);}
		void initMixtureAssignmentTrace(unsigned samples, unsigned num_genes);
		void initMixtureProbabilitesTrace(int samples);
		void initExpressionAcceptanceRatioTrace(unsigned samples, unsigned num_genes);

		void updateExpressionTrace(unsigned sample, unsigned geneIndex)
		{
			for(unsigned category = 0; category < numSelectionCategories; category++)
			{
				expressionTrace[category][sample][geneIndex] = currentExpressionLevel[category][geneIndex];
			}
		}
		void updateCodonSpecificParameterTrace(unsigned sample, char aa);
		void updateSphiTrace(unsigned sample) {sPhiTrace[sample] = Sphi;}
		void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex) {mixtureAssignmentTrace[sample][geneIndex] = mixtureAssignment[geneIndex];}
		std::vector<double> getAPhiTrace() {return aPhiTrace;}
		std::vector<std::vector <unsigned>> getMixtureAssignmentTrace() {return mixtureAssignmentTrace;}
		std::vector<double> getExpectedPhiTrace();

		// functions to return estimates
		double getExpressionPosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement);
		double getSphiPosteriorMean(unsigned samples);
		//double getMixtureAssignmentPosteriorMean(unsigned samples, unsigned geneIndex); // TODO: implement variance function, fix Mean function (won't work with 3 groups)
		unsigned getEstimatedMixtureAssignment(unsigned samples, unsigned geneIndex);
		std::vector <double> getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex);
		double getMutationPosteriorMean(unsigned mixtureElement, unsigned samples, unsigned paramIndex);
		double getSelectionPosteriorMean(unsigned mixtureElement, unsigned samples, unsigned paramIndex);
		double getMutationVariance(unsigned mixtureElement, unsigned samples, unsigned paramIndex, bool unbiased = true);
		double getSelectionVariance(unsigned mixtureElement, unsigned samples, unsigned paramIndex, bool unbiased = true);
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
				unsigned codonIndex = SequenceSummary::CodonToIndex(codon, true);
				rv = getMutationPosteriorMean(mixtureElement - 1, samples, codonIndex);
			}
			return rv;
		}
		double getSelectionPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon)
		{
			double rv = -1.0;
			bool check = checkIndex(mixtureElement, 1, numMixtures);
			if (check)
			{
				unsigned codonIndex = SequenceSummary::CodonToIndex(codon, true);
				rv = getSelectionPosteriorMean(mixtureElement - 1, samples, codonIndex);
			}
			return rv;
		}
		double getMutationVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased)
		{
			double rv = -1.0;
			bool check = checkIndex(mixtureElement, 1, numMixtures);
			if (check)
			{
				unsigned codonIndex = SequenceSummary::CodonToIndex(codon, true);
				rv = getMutationVariance(mixtureElement - 1, samples, codonIndex, unbiased);
			}
			return rv;
		}
		double getSelectionVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, bool unbiased)
		{
			double rv = -1.0;
			bool check = checkIndex(mixtureElement, 1, numMixtures);
			if (check)
			{
				unsigned codonIndex = SequenceSummary::CodonToIndex(codon, true);
				rv = getSelectionVariance(mixtureElement - 1, samples, codonIndex, unbiased);
			}
			return rv;
		}

		void initializeExpressionByGenome(Genome& genome, double sd_phi) {InitializeExpression(genome, sd_phi);}
		void initializeExpressionByList(double sd_phi) {InitializeExpression(sd_phi);}
		void initializeExpressionByRandom(std::vector<double> expression) {InitializeExpression(expression);}
		std::vector<double> getExpressionTraceForGene(int geneIndex) 
		{
			std::vector<double> RV;
			bool check = checkIndex(geneIndex, 1, expressionTrace[0][0].size());
			if (check)
			{ 
				RV = getExpressionTrace(geneIndex - 1);
			}
			return RV;
		}
		std::vector<double> getExpressionTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex)
		{
			std::vector<double> RV;
			bool checkGene = checkIndex(geneIndex, 1, expressionTrace[0][0].size());
			bool checkMixtureElement = checkIndex(mixtureElement, 1, numMixtures);
			if (checkGene && checkMixtureElement)
			{
				unsigned category = getExpressionCategory(mixtureElement - 1);
				unsigned samples = expressionTrace[category].size();
				for (unsigned i = 0u; i < samples; i++)
				{
					RV.push_back(expressionTrace[category][i][geneIndex - 1]);
				}
			}
			return RV;
		}
		std::vector<double> getMutationParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon)
		{
			std::vector<double> RV;
			bool check = checkIndex(mixtureElement, 1, numMixtures);
			if (check)
			{	
				unsigned category = getMutationCategory(mixtureElement - 1);
				unsigned codonIndex = SequenceSummary::CodonToIndex(codon, true);
				unsigned samples = mutationParameterTrace[category].size();
				for (unsigned i = 0u; i < samples; i++)
				{
					RV.push_back(mutationParameterTrace[category][i][codonIndex]);
				}
			}
			return RV;
		}
		std::vector<double> getSelectionParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon)
		{
			std::vector<double> RV;
			bool check = checkIndex(mixtureElement, 1, numMixtures);
			if (check)
			{	
				unsigned category = getSelectionCategory(mixtureElement - 1);
				unsigned codonIndex = SequenceSummary::CodonToIndex(codon, true);
				unsigned samples = selectionParameterTrace[category].size();
				for (unsigned i = 0u; i < samples; i++)
				{
					RV.push_back(selectionParameterTrace[category][i][codonIndex]);
				}
			}
			return RV;
		}
		std::vector<unsigned> getMixtureAssignmentTraceForGene(int geneIndex)
		{
			std::vector <unsigned> RV;
			bool check = checkIndex(geneIndex, 1, mixtureAssignmentTrace[0].size());
			if (check)
			{
				unsigned samples = mixtureAssignmentTrace.size();
				for(unsigned i = 0u; i < samples; i++)
				{
					RV.push_back(mixtureAssignmentTrace[i][geneIndex - 1] + 1);
				}
			}
			return RV;
		}

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
		std::vector<double> getMixtureProbabilitiesTraceForCategory(unsigned categoryIndex)
		{
			std::vector<double> RV;
			bool check = checkIndex(categoryIndex, 1, mixtureProbabilitiesTrace.size());
			if (check)
			{
				RV = mixtureProbabilitiesTrace[categoryIndex - 1];
			}

			return RV;
		}

		std::vector<double> getExpressionAcceptanceRatioTraceByMixtureElementForGeneR(unsigned mixtureElement, unsigned geneIndex) 
		{
			std::vector<double> RV;
			bool checkGene = checkIndex(geneIndex, 1, expressionAcceptanceRatioTrace.size());
			bool checkMixtureElement = checkIndex(mixtureElement, 1, numMixtures); 
			if (checkGene && checkMixtureElement)
			{
				RV = getExpressionAcceptanceRatioTraceByMixtureElementForGene(mixtureElement - 1, geneIndex - 1);
			}
			return RV;
		}
		double getExpressionPosteriorMeanByMixtureElementForGene(unsigned samples, unsigned geneIndex, unsigned mixtureElement)
		{
			double rv = -1.0;
			bool checkGene = checkIndex(geneIndex, 1, expressionTrace[0][0].size());
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
			bool checkGene = checkIndex(geneIndex, 1, expressionTrace[0][0].size());
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
