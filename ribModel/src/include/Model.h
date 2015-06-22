#ifndef MODEL_H
#define MODEL_H
#include "Genome.h"
class Model
{
    private:


        virtual void obtainCodonCount(SequenceSummary& seqsum, char curAA, int codonCount[]) =0;
        virtual double calculateLogLikelihoodPerAAPerGene(unsigned numCodons, int codonCount[], double mutation[], double selection[], double phiValue)=0;
    public:
        explicit Model();
        virtual ~Model();

        virtual void calculateCodonProbabilityVector(unsigned numCodons, double* mutation, double* selection, double phi, double* codonProb)=0;
        // Likelihood ratio functions
        virtual void calculateLogLiklihoodRatioPerGene(Gene& gene, int geneIndex, unsigned k, double* logProbabilityRatio)=0;
        virtual void calculateLogLikelihoodRatioPerAAPerCategory(char curAA, Genome& genome, double& logAcceptanceRatioForAllMixtures)=0;
				virtual void initTraces(unsigned samples, unsigned num_genes, unsigned adaptiveSamples) =0;
    		virtual double getSphi(bool proposed = false) =0;
				virtual double getSphiProposalWidth() =0;
				virtual unsigned getNumMixtureElements() =0;
				virtual double getCategoryProbability(unsigned i) =0;
				virtual void proposeCodonSpecificParameter() =0;
				virtual void adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth) =0;
				virtual void updateCodonSpecificParameter(char aa) =0;
				virtual void updateCodonSpecificParameterTrace(unsigned sample, char aa) =0;
				virtual void proposeSPhi() =0;
				virtual unsigned getMixtureAssignment(unsigned index) =0;
				virtual unsigned getSynthesisRateCategory(unsigned mixture) =0;
				virtual unsigned getSelectionCategory(unsigned mixture) =0;
				virtual unsigned getMutationCategory(unsigned mixture) =0;
				virtual double getSynthesisRate(unsigned index, unsigned mixture, bool proposed = false) =0; 
				virtual double getCurrentSphiProposalWidth() =0;
				virtual double getPreviousSphiProposalWidth() =0;
				virtual void updateSphi() =0;
				virtual void updateSphiTrace(unsigned sample) =0;
				virtual void adaptSphiProposalWidth(unsigned adaptiveWidth) =0;
				virtual void proposeSynthesisRateLevels() =0;
				virtual unsigned getNumSynthesisRateCategories() =0;
				virtual std::vector<unsigned> getMixtureElementsOfSelectionCategory(unsigned k) =0;
				virtual void updateSynthesisRate(unsigned i, unsigned k) =0;
				virtual void setMixtureAssignment(unsigned i, unsigned catOfGene) =0;
				virtual void updateSynthesisRateTrace(unsigned sample, unsigned i) =0;
				virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned i) =0;
				virtual void setCategoryProbability(unsigned mixture, double value) =0;
				virtual void updateMixtureProbabilitiesTrace(unsigned sample) =0;
				virtual void adaptSynthesisRateProposalWidth(unsigned adaptiveWidth) =0;
				virtual void getParameterForCategory(unsigned category, unsigned param, char aa, bool proposal, double* returnValue) =0;
		protected:
};

#endif // MODEL_H
