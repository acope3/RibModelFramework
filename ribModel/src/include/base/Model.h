#ifndef MODEL_H
#define MODEL_H


#include "../Genome.h"


class Model
{
    private:


    public:
		//Constructors & Destructors:
        explicit Model();
        virtual ~Model();


        //Likelihood ratio functions:
        virtual void calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k, double* logProbabilityRatio) = 0;
        virtual void calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome,
        		double& logAcceptanceRatioForAllMixtures) = 0;
		virtual void calculateLogLikelihoodRatioForHyperParameters(Genome &genome, unsigned iteration, std::vector <double> &logProbabilityRatio) = 0;


		//Other functions:
		virtual void simulateGenome(Genome &genome) =0;
		virtual void updateGibbsSampledHyperParameters(Genome &genome) = 0;

		//Parameter wrapper functions:
		virtual void initTraces(unsigned samples, unsigned num_genes) = 0;
		virtual void writeRestartFile(std::string filename) = 0;
		virtual double getSphi(bool proposed = false) = 0;
		virtual unsigned getNumPhiGroupings() = 0;
		virtual void setNumPhiGroupings(unsigned value) = 0;
		virtual unsigned getNumMixtureElements() = 0;
		virtual double getCategoryProbability(unsigned i) = 0;
		virtual void proposeCodonSpecificParameter() = 0;
		virtual void adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth) = 0;
		virtual void updateCodonSpecificParameter(std::string grouping) = 0;
		virtual void updateCodonSpecificParameterTrace(unsigned sample, std::string grouping) = 0;
		virtual void proposeHyperParameters() = 0;
		virtual unsigned getMixtureAssignment(unsigned index) = 0;
		virtual unsigned getSynthesisRateCategory(unsigned mixture) = 0;
		virtual unsigned getSelectionCategory(unsigned mixture) = 0;
		virtual unsigned getMutationCategory(unsigned mixture) = 0;
		virtual double getSynthesisRate(unsigned index, unsigned mixture, bool proposed = false) = 0;
		virtual double getCurrentSphiProposalWidth() = 0;
		virtual void updateHyperParameter(unsigned hp) = 0;
		virtual void updateSphi() = 0;
		virtual void updateSphiTrace(unsigned sample) = 0;
		virtual void updateHyperParameterTraces(unsigned sample) = 0;
		virtual void adaptSphiProposalWidth(unsigned adaptiveWidth) = 0;
		virtual void adaptHyperParameterProposalWidths(unsigned adaptiveWidth) = 0;
		virtual void proposeSynthesisRateLevels() = 0;
		virtual unsigned getNumSynthesisRateCategories() = 0;
		virtual std::vector<unsigned> getMixtureElementsOfSelectionCategory(unsigned k) = 0;
		virtual void updateSynthesisRate(unsigned i, unsigned k) = 0;
		virtual void setMixtureAssignment(unsigned i, unsigned catOfGene) = 0;
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned i) = 0;
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned i) = 0;
		virtual void setCategoryProbability(unsigned mixture, double value) = 0;
		virtual void updateMixtureProbabilitiesTrace(unsigned sample) = 0;
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptiveWidth) = 0;
		//virtual void getParameterForCategory(unsigned category, unsigned param, std::string aa, bool proposal, double* returnValue) = 0;
		virtual unsigned getGroupListSize() = 0;
		virtual std::string getGrouping(unsigned index) = 0;

	protected:
};

#endif // MODEL_H
