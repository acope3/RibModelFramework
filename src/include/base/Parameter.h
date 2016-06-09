#ifndef PARAMETER_H
#define PARAMETER_H


#include "../Genome.h"
#include "../CovarianceMatrix.h"
#include "Trace.h"


#include <vector>
#include <random>
#include <string>
#include <set>
#include <fstream>
#include <ctime>
#include <sstream>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

/* Note 1) -- on getSelectionCategory and getSynthesisRateCategory
 * These two functions are technically the same for readability.
 * Selection and Synthesis Rate are directly related even if they are not known
 * and thus are represented by the same variable. By splitting this
 * into Selection and Synthesis, avoids confusing the two, however.
*/

class Parameter {
	private:

		//STATICS - Sorting Functions:
		void quickSortPair(double a[], int b[], int first, int last);
		void quickSort(double a[], int first, int last);
		static int pivotPair(double a[], int b[], int first, int last);
		static int pivot(double a[], int first, int last);
		static void swap(double& a, double& b);
		static void swap(int& a, int& b);

		unsigned adaptiveStepPrev;
		unsigned adaptiveStepCurr;


		std::vector<double> codonSpecificPrior;
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
			bool splitSer = true, std::string _mutationSelectionState = "allUnique"); //Mostly tested; TODO caveats
		void initBaseValuesFromFile(std::string filename);
		void writeBasicRestartFile(std::string filename);
		void initCategoryDefinitions(std::string mutationSelectionState,
			std::vector<std::vector<unsigned>> mixtureDefinitionMatrix);
		void InitializeSynthesisRate(Genome& genome, double sd_phi);
		void InitializeSynthesisRate(double sd_phi);
		void InitializeSynthesisRate(std::vector<double> expression);
		std::vector<double> readPhiValues(std::string filename); //General function, possibly move


		//prior functions
		double getCodonSpecificPriorStdDev(unsigned paramType);


		//Mixture Definition Matrix and Category Functions:
		void setNumMutationSelectionValues(std::string mutationSelectionState,
			std::vector<std::vector<unsigned>> mixtureDefinitionMatrix);
		void printMixtureDefinitionMatrix();
		double getCategoryProbability(unsigned mixtureElement); //Tested
		void setCategoryProbability(unsigned mixtureElement, double value); //Tested
		unsigned getNumMutationCategories(); //Tested; TODO caveat
		unsigned getNumSelectionCategories(); //Tested; TODO caveat
		unsigned getNumSynthesisRateCategories();
		unsigned getMutationCategory(unsigned mixtureElement); //Tested
		unsigned getSelectionCategory(unsigned mixtureElement); //Tested; see Note 1) at top of file.
		unsigned getSynthesisRateCategory(unsigned mixtureElement); //Tested; see Note 1) at top of file.
		std::vector<unsigned> getMixtureElementsOfMutationCategory(unsigned category); //Tested; TODO caveat
		std::vector<unsigned> getMixtureElementsOfSelectionCategory(unsigned category); //Tested; TODO caveat
		std::string getMutationSelectionState(); //Tested
		unsigned getNumAcceptForCspForIndex(unsigned i); //Tested; only for unit testing.


		//Group List Functions: All tested
		void setGroupList(std::vector<std::string> gl);
		std::string getGrouping(unsigned index);
		std::vector<std::string> getGroupList();
		unsigned getGroupListSize();


		//stdDevSynthesisRate Functions:
		double getStdDevSynthesisRate(unsigned selectionCategory, bool proposed = false); //Tested
		virtual void proposeStdDevSynthesisRate();
		void setStdDevSynthesisRate(double stdDevSynthesisRate, unsigned selectionCategory); //Tested
		double getCurrentStdDevSynthesisRateProposalWidth(); //Tested
		unsigned getNumAcceptForStdDevSynthesisRate(); //Tested; only for unit testing.
		void updateStdDevSynthesisRate();
		double getStdCspForIndex(unsigned i); //Tested; only for unit testing.


		//Synthesis Rate Functions:
		double getSynthesisRate(unsigned geneIndex, unsigned mixtureElement, bool proposed = false); //Tested
		double getCurrentSynthesisRateProposalWidth(unsigned expressionCategory, unsigned geneIndex); //Tested
		double getSynthesisRateProposalWidth(unsigned geneIndex, unsigned mixtureElement); //Tested
		void proposeSynthesisRateLevels();
		void setSynthesisRate(double phi, unsigned geneIndex, unsigned mixtureElement); //Tested
		void updateSynthesisRate(unsigned geneIndex);
		void updateSynthesisRate(unsigned geneIndex, unsigned mixtureElement);
		unsigned getNumAcceptForSynthesisRate(unsigned expressionCategory, unsigned geneIndex); //Tested; only for unit testing


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
		void adaptStdDevSynthesisRateProposalWidth(unsigned adaptationWidth, bool adapt);
		void adaptSynthesisRateProposalWidth(unsigned adaptationWidth, bool adapt);
		virtual void adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth, unsigned lastIteration, bool adapt);


		//Posterior, Variance, and Estimates Functions:
		double getStdDevSynthesisRatePosteriorMean(unsigned samples, unsigned mixture);
		double getSynthesisRatePosteriorMean(unsigned samples, unsigned geneIndex, unsigned mixtureElement);

		double getCodonSpecificPosteriorMean(unsigned mixtureElement, unsigned samples, std::string &codon, unsigned paramType,
			bool withoutReference = true);
		double getStdDevSynthesisRateVariance(unsigned samples, unsigned mixture, bool unbiased);
		double getSynthesisRateVariance(unsigned samples, unsigned geneIndex, unsigned mixtureElement,
			bool unbiased = true);
		double getCodonSpecificVariance(unsigned mixtureElement, unsigned samples, std::string &codon, unsigned paramType,
			bool unbiased, bool withoutReference = true);
        std::vector<double> getCodonSpecificQuantile(unsigned mixtureElement, unsigned samples, std::string &codon,
			unsigned paramType, std::vector<double> probs, bool withoutReference);
		unsigned getEstimatedMixtureAssignment(unsigned samples, unsigned geneIndex);
		std::vector<double> getEstimatedMixtureAssignmentProbabilities(unsigned samples, unsigned geneIndex);


		//Other Functions:
		unsigned getNumParam(); //Tested
		unsigned getNumMixtureElements(); //Tested //TODO style question: Why call it getNumMixtureElements instead of simply getNumMixtures? Alternatively, change name of variable.
		unsigned getNumObservedPhiSets();
		void setMixtureAssignment(unsigned gene, unsigned value); //Tested
		unsigned getMixtureAssignment(unsigned gene); //Tested
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
		void setCategoriesForTrace();
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
		double getCodonSpecificPosteriorMeanForCodon(unsigned mixtureElement, unsigned samples, std::string codon, unsigned paramType,
			bool withoutReference);
		double getCodonSpecificVarianceForCodon(unsigned mixtureElement, unsigned samples, std::string codon, unsigned paramType, bool unbiased,
			bool withoutReference);
        std::vector<double> getCodonSpecificQuantileForCodon(unsigned mixtureElement, unsigned samples, std::string &codon, unsigned paramType, std::vector<double> probs,
	       bool withoutReference);


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

		std::vector<std::vector<std::vector<double>>> proposedCodonSpecificParameter;
		std::vector<std::vector<std::vector<double>>> currentCodonSpecificParameter;

		std::vector<unsigned> mixtureAssignment;
		std::vector<std::string> groupList;
		unsigned maxGrouping;


		std::vector<double> stdDevSynthesisRate_proposed;
		std::vector<double> stdDevSynthesisRate;
		double bias_stdDevSynthesisRate; //NOTE: Currently, this value is always set to 0.0
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

		double bias_phi; //NOTE: Currently, this value is always set to 0.0
		std::vector<std::vector<double>> std_phi;

};

#endif // PARAMETER_H
