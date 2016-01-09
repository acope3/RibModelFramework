#ifndef PARAMETER_H
#define PARAMETER_H

#include <vector>
#include <random>
#include <string>
#include <iostream>
#include <set>
#include <fstream>
#include <ctime>
#include <sstream>
#ifndef STANDALONE
#include <Rcpp.h>
#endif

#include "../Genome.h"
#include "../CovarianceMatrix.h"
#include "Trace.h"



class Parameter {
	private:

		//STATICS - Sorting Functions:
		void quickSortPair(double a[], int b[], int first, int last);
		void quickSort(double a[], int first, int last);
		static int pivotPair(double a[], int b[], int first, int last);
		static int pivot(double a[], int first, int last);
		static void swap(double& a, double& b);
		static void swap(int& a, int& b);

	public:

		static const std::string allUnique;
		static const std::string selectionShared;
		static const std::string mutationShared;

		static const unsigned dM;
		static const unsigned dEta;
		static const unsigned dOmega;
		static const unsigned alp;
		static const unsigned lmPri;

#ifdef STANDALONE
		static std::default_random_engine generator; // static to make sure that the same generator is during the runtime.
#endif



		//Constructors & Destructors:
		Parameter();
		Parameter(unsigned maxGrouping);
		Parameter& operator=(const Parameter& rhs);
		virtual ~Parameter();



		//Initialization and Restart Functions:
		void initParameterSet(std::vector<double> stdDevSynthesisRate, unsigned _numMixtures, std::vector<unsigned> geneAssignment,
							  std::vector<std::vector<unsigned>> mixtureDefinitionMatrix,
							  bool splitSer = true, std::string _mutationSelectionState = "allUnique");
		void initBaseValuesFromFile(std::string filename);
		void writeBasicRestartFile(std::string filename);
		void initCategoryDefinitions(std::string mutationSelectionState,
								 std::vector<std::vector<unsigned>> mixtureDefinitionMatrix);
		void InitializeSynthesisRate(Genome& genome, double sd_phi);
		void InitializeSynthesisRate(double sd_phi);
		void InitializeSynthesisRate(std::vector<double> expression);
		std::vector<double> readPhiValues(std::string filename); //General function, possibly move



		//Mixture Definition Matrix and Category Functions:
		void setNumMutationSelectionValues(std::string mutationSelectionState,
									   std::vector<std::vector<unsigned>> mixtureDefinitionMatrix);
		void printMixtureDefinitionMatrix();
		double getCategoryProbability(unsigned mixtureElement);
		void setCategoryProbability(unsigned mixtureElement, double value);
		unsigned getNumMutationCategories();
		unsigned getNumSelectionCategories();
		unsigned getNumSynthesisRateCategories();
		unsigned getMutationCategory(unsigned mixtureElement);
		unsigned getSelectionCategory(unsigned mixtureElement); //TODO: Add comments explaining reasonsing here for same function
		unsigned getSynthesisRateCategory(unsigned mixtureElement);
		std::vector<unsigned> getMixtureElementsOfMutationCategory(unsigned category);
		std::vector<unsigned> getMixtureElementsOfSelectionCategory(unsigned category);
		std::string getMutationSelectionState();



		//Group List Functions:
		void setGroupList(std::vector<std::string> gl);
		std::string getGrouping(unsigned index);
		std::vector<std::string> getGroupList();
		unsigned getGroupListSize();



		//stdDevSynthesisRate Functions:
		double getStdDevSynthesisRate(unsigned selectionCategory, bool proposed = false);
		virtual void proposeStdDevSynthesisRate();
		void setStdDevSynthesisRate(double stdDevSynthesisRate, unsigned selectionCategory);
		double getCurrentStdDevSynthesisRateProposalWidth();
		void updateStdDevSynthesisRate();



		//Synthesis Rate Functions:
		double getSynthesisRate(unsigned geneIndex, unsigned mixtureElement, bool proposed = false);
		double getCurrentSynthesisRateProposalWidth(unsigned expressionCategory, unsigned geneIndex);
		double getSynthesisRateProposalWidth(unsigned geneIndex, unsigned mixtureElement);
		void proposeSynthesisRateLevels();
		void setSynthesisRate(double phi, unsigned geneIndex, unsigned mixtureElement);
		void updateSynthesisRate(unsigned geneIndex);
		void updateSynthesisRate(unsigned geneIndex, unsigned mixtureElement);



		//Iteration Functions:
		unsigned getLastIteration();
		void setLastIteration(unsigned iteration);



		//Trace Functions:
		Trace& getTraceObject();
		void setTraceObject(Trace _trace);
		void updateStdDevSynthesisRateTrace(unsigned sample);
		void updateSynthesisRateTrace(unsigned sample, unsigned geneIndex);
		void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex);
		void updateMixtureProbabilitiesTrace(unsigned samples);


		//Adaptive Width Functions:
		void adaptStdDevSynthesisRateProposalWidth(unsigned adaptationWidth);
		void adaptSynthesisRateProposalWidth(unsigned adaptationWidth);
		virtual void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth);


		//Posterior, Variance, and Estimates Functions:
		double getStdDevSynthesisRatePosteriorMean(unsigned samples, unsigned mixture);
		double getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement);

		double getCodonSpecificPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon, unsigned paramType);
		double getStdDevSynthesisRateVariance(unsigned samples, unsigned mixture, bool unbiased);
		double getSynthesisRateVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement,
											bool unbiased = true);
		double getCodonSpecificVariance(unsigned mixtureElement, unsigned samples, std::string &codon, unsigned paramType, bool unbiased);
		unsigned getEstimatedMixtureAssignment(unsigned samples, unsigned geneIndex);
		std::vector<double> getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex);



		//Other Functions:
		unsigned getNumParam();
		unsigned getNumMixtureElements();
		unsigned getNumObservedPhiSets();
		void setMixtureAssignment(unsigned gene, unsigned value);
		unsigned getMixtureAssignment(unsigned gene);
		virtual void setNumObservedPhiSets(unsigned _phiGroupings);
		virtual std::vector <std::vector <double> > calculateSelectionCoefficients(unsigned sample, unsigned mixture);

		

		//Static Functions:
		static double calculateSCUO(Gene& gene, unsigned maxAA);
		static void drawIidRandomVector(unsigned draws, double mean, double sd, double (*proposal)(double a, double b),
				double* randomNumbers);
		static void drawIidRandomVector(unsigned draws, double r, double (*proposal)(double r), double* randomNumber);
		static double randNorm(double mean, double sd);
		static double randLogNorm(double m, double s);
		static double randExp(double r);
		static double randGamma(double shape, double rate);
		static void randDirichlet(double *input, unsigned numElements, double *output);
		static double randUnif(double minVal, double maxVal);
		static unsigned randMultinom(double* probabilities, unsigned mixtureElements);
		static double densityNorm(double x, double mean, double sd, bool log = false);
		static double densityLogNorm(double x, double mean, double sd, bool log = false);
		//double getMixtureAssignmentPosteriorMean(unsigned samples, unsigned geneIndex); // TODO: implement variance function, fix Mean function (won't work with 3 groups)







		//R Section:

#ifndef STANDALONE

		//Initialization and Restart Functions:
		void initializeSynthesisRateByGenome(Genome& genome, double sd_phi);
		void initializeSynthesisRateByRandom(double sd_phi);
		void initializeSynthesisRateByList(std::vector<double> expression);
		bool checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound);



		//Mixture Definition Matrix and Category Functions:
		unsigned getMutationCategoryForMixture(unsigned mixtureElement);
		unsigned getSelectionCategoryForMixture(unsigned mixtureElement);
		unsigned getSynthesisRateCategoryForMixture(unsigned mixtureElement);
		std::vector<std::vector<unsigned>> getCategories();
		void setCategories(std::vector<std::vector<unsigned>> _categories);
		void setNumMutationCategories(unsigned _numMutationCategories);
		void setNumSelectionCategories(unsigned _numSelectionCategories);



		//Synthesis Rate Functions:
		std::vector<std::vector<double>> getSynthesisRateR();
		std::vector<double> getCurrentSynthesisRateForMixture(unsigned mixture);



		//Posterior, Variance, and Estimates Functions:
		double getSynthesisRatePosteriorMeanByMixtureElementForGene(unsigned samples, unsigned geneIndex,
																	unsigned mixtureElement);
		double getSynthesisRateVarianceByMixtureElementForGene(unsigned samples, unsigned geneIndex,
														   unsigned mixtureElement, bool unbiased);
		unsigned getEstimatedMixtureAssignmentForGene(unsigned samples, unsigned geneIndex);
		std::vector<double> getEstimatedMixtureAssignmentProbabilitiesForGene(unsigned samples, unsigned geneIndex);
		double getCodonSpecificPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon, unsigned paramType);
		double getCodonSpecificVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, unsigned paramType, bool unbiased);


		//Other Functions:
		SEXP calculateSelectionCoefficientsR(unsigned sample, unsigned mixture);
		std::vector<unsigned> getMixtureAssignmentR();
		void setMixtureAssignmentR(std::vector<unsigned> _mixtureAssignment);
		unsigned getMixtureAssignmentForGeneR(unsigned geneIndex);
		void setMixtureAssignmentForGene(unsigned geneIndex, unsigned value);
		void setNumMixtureElements(unsigned _numMixtures);

#endif

	protected:
		Trace traces;

		std::vector<CovarianceMatrix> covarianceMatrix;
		std::vector<mixtureDefinition> categories;
		std::vector<double> categoryProbabilities;
		std::vector<std::vector<unsigned>> mutationIsInMixture;
		std::vector<std::vector<unsigned>> selectionIsInMixture;
		unsigned numMutationCategories; //TODO Probably needs to be renamed
		unsigned numSelectionCategories; //TODO Probably needs to be renamed
		std::vector<unsigned> numAcceptForCodonSpecificParameters;
		std::string mutationSelectionState; //TODO: Probably needs to be renamed


		std::vector<unsigned> mixtureAssignment;
		std::vector<std::string> groupList;
		unsigned maxGrouping;


		std::vector<double> stdDevSynthesisRate_proposed;
		std::vector<double> stdDevSynthesisRate;
		double bias_stdDevSynthesisRate;
		double std_stdDevSynthesisRate;
		unsigned numAcceptForStdDevSynthesisRate;
		std::vector<double> std_csp;



		std::vector<std::vector<double>> proposedSynthesisRateLevel;
		std::vector<std::vector<double>> currentSynthesisRateLevel;
		std::vector<std::vector<unsigned>> numAcceptForSynthesisRate;

		unsigned lastIteration;

		unsigned int numParam;
		unsigned numMixtures;
		unsigned obsPhiSets;

		double bias_phi;
		std::vector<std::vector<double>> std_phi;

};

#endif // PARAMETER_H
