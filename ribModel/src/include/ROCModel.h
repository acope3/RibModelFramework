#ifndef ROCMODEL_H
#define ROCMODEL_H

#include "Model.h"
#include "ROCParameter.h"

class ROCModel : public Model
{
    private:

				ROCParameter *parameter;
        virtual void obtainCodonCount(SequenceSummary& seqsum, char curAA, int codonCount[]);
        virtual double calculateLogLikelihoodPerAAPerGene(unsigned numCodons, int codonCount[], double mutation[], double selection[], double phiValue);
    public:
        ROCModel();
        virtual ~ROCModel();
        ROCModel(const ROCModel& other);
				void setParameter(ROCParameter &_parameter);
        virtual void calculateCodonProbabilityVector(unsigned numCodons, double* mutation, double* selection, double phi, double* codonProb);
        // Likelihood ratio functions
        virtual void calculateLogLiklihoodRatioPerGene(Gene& gene, int geneIndex, unsigned k, double* logProbabilityRatio);
        virtual void calculateLogLikelihoodRatioPerAAPerCategory(char curAA, Genome& genome, double& logAcceptanceRatioForAllMixtures);
				virtual void initTraces(unsigned samples, unsigned num_genes, unsigned adaptiveSamples) {parameter -> initAllTraces(samples, num_genes, adaptiveSamples);}
       
				virtual double getSphi(bool proposed = false) {return parameter->getSphi(proposed);}
				virtual double getSphiProposalWidth() {return parameter->getSphiProposalWidth();}
				virtual unsigned getNumMixtureElements() {return parameter->getNumMixtureElements();}
				virtual unsigned getCategoryProbability(unsigned i) {return parameter->getCategoryProbability(i);}
				virtual void proposeCodonSpecificParameter() {parameter->proposeCodonSpecificParameter();}
				virtual void adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth) {parameter->adaptCodonSpecificParameterProposalWidth(adaptiveWidth);}
				virtual void updateCodonSpecificParameter(char aa) {parameter->updateCodonSpecificParameter(aa);}
				virtual void updateCodonSpecificParameterTrace(unsigned sample, char aa) {parameter->updateCodonSpecificParameterTrace(sample, aa);}
				virtual void proposeSPhi() {parameter->proposeSPhi();}
				virtual unsigned getMixtureAssignment(unsigned index) {return parameter->getMixtureAssignment(index);}
				virtual unsigned getSynthesisRateCategory(unsigned mixture) {return parameter->getSynthesisRateCategory(mixture);}
        virtual unsigned getSelectionCategory(unsigned mixture) {return parameter ->getSelectionCategory(mixture);}
        virtual unsigned getMutationCategory(unsigned mixture)  {return parameter ->getMutationCategory(mixture);} 				
				virtual double getSynthesisRate(unsigned index, unsigned mixture, bool proposed = false) {return parameter->getSynthesisRate(index, mixture, false);}
				virtual double getCurrentSphiProposalWidth() {return parameter->getCurrentSphiProposalWidth();}
				virtual double getPreviousSphiProposalWidth() {return parameter->getPreviousSphiProposalWidth();}
				virtual void updateSphi() {parameter->updateSphi();}
				virtual void updateSphiTrace(unsigned sample) {parameter->updateSphiTrace(sample);}
				virtual void adaptSphiProposalWidth(unsigned adaptiveWidth) {parameter->adaptSphiProposalWidth(adaptiveWidth);}
				virtual void proposeSynthesisRateLevels() {parameter->proposeSynthesisRateLevels();}
				virtual unsigned getNumSynthesisRateCategories() {return parameter->getNumSynthesisRateCategories();}
				virtual std::vector<unsigned> getMixtureElementsOfSelectionCategory(unsigned k) {return parameter->getMixtureElementsOfSelectionCategory(k);}
				virtual void updateSynthesisRate(unsigned i, unsigned k) {parameter->updateSynthesisRate(i,k);}
				virtual void setMixtureAssignment(unsigned i, unsigned catOfGene) {parameter->setMixtureAssignment(i, catOfGene);}
        virtual void updateSynthesisRateTrace(unsigned sample, unsigned i) {parameter->updateSynthesisRateTrace(sample, i);}
        virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned i) {parameter->updateMixtureAssignmentTrace(sample, i);}
        virtual void setCategoryProbability(unsigned mixture, double value) {parameter->setCategoryProbability(mixture, value);}
        virtual void updateMixtureProbabilitiesTrace(unsigned sample) {parameter->updateMixtureProbabilitiesTrace(sample);}
				virtual void adaptSynthesisRateProposalWidth(unsigned adaptiveWidth) {parameter->adaptSynthesisRateProposalWidth(adaptiveWidth);}
				virtual void getParameterForCategory(unsigned category, unsigned param, char aa, bool proposal, double* returnValue)
				{
					parameter -> getParameterForCategory(category, param, aa, proposal, returnValue);
				}
				 // R wrapper
        std::vector<double> CalculateProbabilitiesForCodons(std::vector<double> mutation, std::vector<double> selection, double phi);

    protected:
};

#endif // ROCMODEL_H
