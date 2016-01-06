#ifndef PANSEMODEL_H
#define PANSEMODEL_H


#include "../base/Model.h"
#include "PANSEParameter.h"


class PANSEModel: public Model {
	private:
		PANSEParameter *parameter;

		virtual double calculateLogLikelihoodPerCodonPerGene(double currAlpha, double currLambdaPrime,
				unsigned currRFPObserved, unsigned currNumCodonsInMRNA, double phiValue);


	public:
		//Constructors & Destructors:
		explicit PANSEModel();
		virtual ~PANSEModel();


		//Likelihood ratio functions
		virtual void calculateLogLikelihoodRatioPerGene(Gene& gene, unsigned geneIndex, unsigned k,
				double* logProbabilityRatio);
		virtual void calculateLogLikelihoodRatioPerGroupingPerCategory(std::string grouping, Genome& genome,
				double& logAcceptanceRatioForAllMixtures);


		//Other functions:
		void setParameter(PANSEParameter &_parameter);
		virtual void simulateGenome(Genome &genome);


		//Parameter wrapper functions:
		virtual void initTraces(unsigned samples, unsigned num_genes)
		{
			parameter->initAllTraces(samples, num_genes);
		}
		virtual void writeRestartFile(std::string filename)
		{
			return parameter->writeEntireRestartFile(filename);
		}
		virtual double getSphi(bool proposed = false)
		{
			return parameter->getSphi(proposed);
		}
		virtual unsigned getNumMixtureElements()
		{
			return parameter->getNumMixtureElements();
		}
		virtual double getCategoryProbability(unsigned i)
		{
			return parameter->getCategoryProbability(i);
		}
		virtual void proposeCodonSpecificParameter()
		{
			parameter->proposeCodonSpecificParameter();
		}
		virtual void adaptCodonSpecificParameterProposalWidth(unsigned adaptiveWidth)
		{
			parameter->adaptCodonSpecificParameterProposalWidth(adaptiveWidth);
		}
		virtual void updateCodonSpecificParameter(std::string aa)
		{
			parameter->updateCodonSpecificParameter(aa);
		}
		virtual void updateCodonSpecificParameterTrace(unsigned sample, std::string codon)
		{
			parameter->updateCodonSpecificParameterTrace(sample, codon);
		}
		virtual void proposeSPhi()
		{
			parameter->proposeSphi();
		}
		virtual unsigned getMixtureAssignment(unsigned index)
		{
			return parameter->getMixtureAssignment(index);
		}
		virtual unsigned getSynthesisRateCategory(unsigned mixture)
		{
			return parameter->getSynthesisRateCategory(mixture);
		}
		virtual unsigned getSelectionCategory(unsigned mixture)
		{
			return parameter->getSelectionCategory(mixture);
		}
		virtual unsigned getMutationCategory(unsigned mixture)
		{
			return parameter->getMutationCategory(mixture);
		}
		virtual double getSynthesisRate(unsigned index, unsigned mixture, bool proposed = false)
		{
			return parameter->getSynthesisRate(index, mixture, proposed);
		}
		virtual double getCurrentSphiProposalWidth()
		{
			return parameter->getCurrentSphiProposalWidth();
		}
		virtual void updateSphi()
		{
			parameter->updateSphi();
		}
		virtual void updateSphiTrace(unsigned sample)
		{
			parameter->updateSphiTrace(sample);
		}
		virtual void adaptSphiProposalWidth(unsigned adaptiveWidth)
		{
			parameter->adaptSphiProposalWidth(adaptiveWidth);
		}
		virtual void proposeSynthesisRateLevels()
		{
			parameter->proposeSynthesisRateLevels();
		}
		virtual unsigned getNumSynthesisRateCategories()
		{
			return parameter->getNumSynthesisRateCategories();
		}
		virtual std::vector<unsigned> getMixtureElementsOfSelectionCategory(unsigned k)
		{
			return parameter->getMixtureElementsOfSelectionCategory(k);
		}
		virtual void updateSynthesisRate(unsigned i, unsigned k)
		{
			parameter->updateSynthesisRate(i, k);
		}
		virtual void setMixtureAssignment(unsigned i, unsigned catOfGene)
		{
			parameter->setMixtureAssignment(i, catOfGene);
		}
		virtual void updateSynthesisRateTrace(unsigned sample, unsigned i)
		{
			parameter->updateSynthesisRateTrace(sample, i);
		}
		virtual void updateMixtureAssignmentTrace(unsigned sample, unsigned i)
		{
			parameter->updateMixtureAssignmentTrace(sample, i);
		}
		virtual void setCategoryProbability(unsigned mixture, double value)
		{
			parameter->setCategoryProbability(mixture, value);
		}
		virtual void updateMixtureProbabilitiesTrace(unsigned sample)
		{
			parameter->updateMixtureProbabilitiesTrace(sample);
		}
		virtual void adaptSynthesisRateProposalWidth(unsigned adaptiveWidth)
		{
			parameter->adaptSynthesisRateProposalWidth(adaptiveWidth);
		}
		virtual double getParameterForCategory(unsigned category, unsigned param, std::string codon, bool proposal)
		{
			return parameter->getParameterForCategory(category, param, codon, proposal);
		}
		virtual unsigned getListSize() //TODO: should not be static
		{
			return 61;
		}
		virtual std::string getGrouping(unsigned index)
		{
			return parameter->getGrouping(index);
		}
		virtual unsigned getGroupListSize() {return parameter->getGroupListSize();}

		// R wrapper
		std::vector<double> CalculateProbabilitiesForCodons(std::vector<double> mutation, std::vector<double> selection,
				double phi);

	protected:
};

#endif //PANSEMODEL_H
