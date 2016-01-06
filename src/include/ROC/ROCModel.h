#ifndef ROCMODEL_H
#define ROCMODEL_H

#include "../base/Model.h"
#include "ROCParameter.h"

class ROCModel : public Model
{
    private:

	ROCParameter *parameter;
	virtual void obtainCodonCount(SequenceSummary *seqsum, std::string curAA, int codonCount[]);
	double calculateMutationPrior(std::string grouping, bool proposed = false); // TODO add to FONSE as well? // cedric
	double calculateLogLikelihoodPerAAPerGene(unsigned numCodons, int codonCount[], double mutation[], double selection[], double phiValue);
	bool withPhi;

    public:

	ROCModel(bool _withPhi = false);
	virtual ~ROCModel();
	void setParameter(ROCParameter &_parameter);
	void calculateCodonProbabilityVector(unsigned numCodons, double mutation[], double selection[], double phi, double codonProb[]);
	// Likelihood ratio functions
	virtual void calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k, double* logProbabilityRatio);
	virtual void calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome, double& logAcceptanceRatioForAllMixtures);
	void simulateGenome(Genome &genome);
	virtual void calculateLogLikelihoodRatioForHyperParameters(Genome &genome, unsigned iteration, std::vector <double> &logProbabilityRatio);
	virtual double calculateAllPriors();

	virtual void updateGibbsSampledHyperParameters(Genome &genome);

	//Parameter wrapper functions:
	virtual void initTraces(unsigned samples, unsigned num_genes) {parameter -> initAllTraces(samples, num_genes);}
	virtual void writeRestartFile(std::string filename) {return parameter->writeEntireRestartFile(filename);}       
	virtual double getSphi(unsigned selectionCategory, bool proposed = false) {return parameter->getStdDevSynthesisRate(selectionCategory, proposed);}
	virtual unsigned getNumPhiGroupings() { return parameter->getNumObservedPhiSets(); }
	virtual void setNumPhiGroupings(unsigned value) { parameter->setNumObservedPhiSets(value); }
	double getAphi(unsigned index, bool proposed = false) { return parameter->getAphi(index, proposed); }
	double getSepsilon(unsigned index) { return parameter->getSepsilon(index); }
	virtual unsigned getNumMixtureElements() {return parameter->getNumMixtureElements();}
	virtual double getCategoryProbability(unsigned i) {return parameter->getCategoryProbability(i);}
	virtual void proposeCodonSpecificParameter() {parameter->proposeCodonSpecificParameter();}
	virtual void adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth) {parameter->adaptCodonSpecificParameterProposalWidth(adaptiveWidth);}
	virtual void updateCodonSpecificParameter(std::string grouping) {parameter->updateCodonSpecificParameter(grouping);}
	virtual void updateCodonSpecificParameterTrace(unsigned sample, std::string grouping) {parameter->updateCodonSpecificParameterTrace(sample,grouping);}
	virtual void proposeHyperParameters();
	virtual unsigned getMixtureAssignment(unsigned index) {return parameter->getMixtureAssignment(index);}
	virtual unsigned getSynthesisRateCategory(unsigned mixture) {return parameter->getSynthesisRateCategory(mixture);}
	virtual unsigned getSelectionCategory(unsigned mixture) {return parameter ->getSelectionCategory(mixture);}
	virtual unsigned getMutationCategory(unsigned mixture)  {return parameter ->getMutationCategory(mixture);}
	virtual double getSynthesisRate(unsigned index, unsigned mixture, bool proposed = false) {return parameter->getSynthesisRate(index, mixture, proposed);}
	virtual double getCurrentSphiProposalWidth() {return parameter->getCurrentStdDevSynthesisRateProposalWidth();}
	double getCurrentAphiProposalWidth(unsigned index) { return parameter->getCurrentAphiProposalWidth(index); }
	virtual void adaptHyperParameterProposalWidths(unsigned adaptiveWidth);
	virtual void updateAllHyperParameter();
	virtual void updateHyperParameter(unsigned hp);
	virtual void updateSphi() {parameter->updateStdDevSynthesisRate();}
	virtual void updateSphiTrace(unsigned sample) {parameter->updateStdDevSynthesisRateTrace(sample);}
	virtual void adaptSphiProposalWidth(unsigned adaptiveWidth) {parameter->adaptStdDevSynthesisRateProposalWidth(adaptiveWidth);}
	void updateAphi(unsigned index) { parameter->updateAphi(index); }
	void updateAphiTrace(unsigned index, unsigned sample) { parameter->updateAphiTraces(sample); }
	void updateSepsilonTrace(unsigned index, unsigned sample) { parameter->updateSepsilonTraces(sample); }
	void adaptAphiProposalWidth(unsigned adaptiveWidth) { parameter->adaptAphiProposalWidth(adaptiveWidth); }
	virtual void updateHyperParameterTraces(unsigned sample);
	virtual void updateTracesWithInitialValues(Genome &genome);
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
	virtual void getParameterForCategory(unsigned category, unsigned param, std::string aa, bool proposal, double* returnValue)
	{
		parameter -> getParameterForCategory(category, param, aa, proposal, returnValue);
	}
	virtual unsigned getGroupListSize() {return parameter->getGroupListSize();} 
	virtual std::string getGrouping(unsigned index) {return parameter -> getGrouping(index);}
	// R wrapper
	std::vector<double> CalculateProbabilitiesForCodons(std::vector<double> mutation, std::vector<double> selection, double phi);

	virtual void setLastIteration(unsigned iteration) { parameter->setLastIteration(iteration); }
	virtual unsigned getLastIteration() { return parameter->getLastIteration(); }
	virtual void updateTmp() {}
	virtual void printHyperParameters();
    protected:
};

#endif // ROCMODEL_H
