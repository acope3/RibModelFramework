#include <iostream>
#include <time.h>
#include <fstream>
#include <algorithm>

#include "include/Genome.h"
#include "include/Gene.h"
#include "include/SequenceSummary.h"
#include "include/ROCParameter.h"
#include "include/ROCModel.h"
#include "include/MCMCAlgorithm.h"
#include "include/CovarianceMatrix.h"

using namespace std;

void testNumCodonsPerAA()
{
	std::cout << "------------------ # CODONS PER AA ------------------" << std::endl;
	for(int i = 0; i < 22; i++)
	{
		char aa = SequenceSummary::AminoAcidArray[i];
		int ncodons = SequenceSummary::GetNumCodonsForAA(aa);
		std::cout << "# codons for " << aa << "\t" << ncodons << std::endl;
	}
	std::cout << "------------------ # CODONS PER AA ------------------" << std::endl;
}

void testCodonRangePerAA(bool forParamVector)
{
	std::cout << "------------------ CODON RANGE PER AA ------------------" << std::endl;
	for(int i = 0; i < 22; i++)
	{
		char aa = SequenceSummary::AminoAcidArray[i];
		unsigned codonRange[2];
		SequenceSummary::AAToCodonRange(aa, forParamVector, codonRange);
		std::cout << "codon range for " << aa << "\t" << codonRange[0] << " - " << codonRange[1] << std::endl;
	}
	std::cout << "------------------ CODON RANGE PER AA ------------------" << std::endl;
}

void testLogNormDensity()
{
	std::cout << "------------------ LOG NORM DENSITY ------------------" << std::endl;
	for(int i = 0; i < 5; i++)
	{
		double result = std::log(ROCParameter::densityLogNorm(i, -1, 1));
		std::cout << "logP of " << i << "\t" << result << std::endl;
	}
	std::cout << "------------------ LOG NORM DENSITY ------------------" << std::endl;
}

void testSCUO(Genome& genome)
{
	std::cout << "------------------ SCUO VALUES ------------------" << std::endl;
	for(int n = 0; n < genome.getGenomeSize(); n++)
	{
		std::cout << genome.getGene(n).getId() << "\t" << ROCParameter::calculateSCUO(genome.getGene(n)) << std::endl;
	}
	std::cout << "------------------ SCUO VALUES ------------------" << std::endl;
}

void testCovarianceMatrix()
{
	std::cout << "------------------ TEST COVARIANCE ROUTINES ------------------" << std::endl;
	//double arr[3][3] { {4,2,5}, {2,2,3}, {5,3,9} };
	std::vector <double> arr {25,15,-5, 15,18,0, -5,0,11};
	CovarianceMatrix covMat = arr;
	std::cout << "Choleski Matrix pre decomposition" << std::endl;
	covMat.printCholeskiMatrix();
	covMat.choleskiDecomposition();
	std::cout << "Covariance Matrix" << std::endl;
	covMat.printCovarianceMatrix();
	std::cout << "Choleski Matrix" << std::endl;
	covMat.printCholeskiMatrix();
	std::cout << "Matrix vector multiplication" << std::endl;
	double brr[3] {1, 0, -1};
	double res[3];
	covMat.transformIidNumersIntoCovaryingNumbers(brr, res);
	for(int i = 0; i < 3; i++)
	{
		std::cout << res[i] << std::endl;
	}
	if(res[0] == 5 && res[1] == 3 && res[2] == -4)
	{
		std::cout << "CHECK" << std::endl;
	}else{
		std::cout << "ERROR" << std::endl;
	}
	std::cout << "------------------ TEST COVARIANCE ROUTINES ------------------" << std::endl;

}

void testGetCountsForAA(Genome genome)
{

	unsigned codonCounts[genome.getGenomeSize()][5];
	genome.getCountsForAA('A', codonCounts);

	for(unsigned i = 0u; 0 < genome.getGenomeSize(); i++)
	{
		Gene gene = genome.getGene(i);
		std::cout << i << "| " << gene.getId() << ": " << codonCounts[i][0] << " " << codonCounts[i][1] << " " << codonCounts[i][2] << " " << codonCounts[i][3] << "\n";
	}
}


void testRandMultiNom(unsigned numCat)
{
	std::cout << "------------------ TEST RANDMULTINOMIAL ------------------" << std::endl;
	double assignmentCounts[3] = {0, 0, 0};
	double probabilities[3] = {0.5, 0.2, 0.3};
	unsigned tmp;
	for (int i = 0; i < 50; i++)
	{
		tmp = ROCParameter::randMultinom(probabilities, numCat);
		assignmentCounts[tmp] += 1;
	}
	std::cout <<"Printing dirichlet numbers:\n";
	for (int i = 0; i < numCat; i++)
	{
		std::cout <<"  " << i <<": " << assignmentCounts[i] <<"\n";
		std::cout <<"  " << i <<"/50: " << assignmentCounts[i]/50 <<"\n\n";
	}

	std::cout << "------------------ TEST RANDMULTINOMIAL ------------------" << std::endl;

}

void testThetaKMatrix()
{
	unsigned matrix[2][2] = { {2,1}, {1,1} };
	std::cout << "------------------ TEST THETAKMATRIX ------------------" << std::endl;
	ROCParameter R(100, 2, 2, nullptr, true, "allUnique");

	R.printThetaKMatrix();
	std::cout <<"numMutationCategories: " << R.getNumMutationCategories() <<"\n";
	std::cout <<"numSelectionCategories: " << R.getNumSelectionCategories() <<"\n";
	std::string files[] = {std::string("Skluyveri_CSP_ChrA.csv"), std::string("Skluyveri_CSP_ChrCleft.csv")};
	std::cout <<"files array good to go\n";
	R.initMutationSelectionCategories(files, R.getNumMutationCategories(), ROCParameter::dM);
	R.initMutationSelectionCategories(files, R.getNumSelectionCategories(), ROCParameter::dEta);


	std::cout <<"looping through all mutation params:\n";
	std::vector<std::vector <double>> temp;
	temp = R.getCurrentMutationParameter();
	for (int i = 0; i < temp.size(); i++)
	{
		std::cout <<"Mutation param " << i <<":\n";
		for (int j = 0; j < temp[i].size(); j++)
		{
			std::cout << temp[i][j] <<"\n";
		}
		std::cout <<"\n\n";
	}

	std::cout <<"looping through all selection params:\n";
	temp = R.getCurrentSelectionParameter();
	for (int i = 0; i < temp.size(); i++)
	{
		std::cout <<"Selection param " << i <<":\n";
		for (int j = 0; j < temp[i].size(); j++)
		{
			std::cout << temp[i][j] <<"\n";
		}
		std::cout <<"\n\n";
	}
	std::cout << "------------------ TEST THETAKMATRIX ------------------" << std::endl;

}

void testSimulateGenome(Genome& genome)
{
	std::cout << "------------------ TEST SIMULATEGENOME ------------------" << std::endl;
	ROCModel model;
	unsigned geneAssignment[genome.getGenomeSize()];
	for(int i = 0; i < genome.getGenomeSize(); i++)
	{
		if(i < 448) geneAssignment[i] = 0u;
		else geneAssignment[i] = 1u;
	}

	std::cout << "initialize ROCParameter object" << std::endl;
	double sphi_init = 2;
	double numMixtures = 2;
	std::string mixDef = ROCParameter::allUnique;
	std::cout << "\tSphi init: " << sphi_init << "\n";
	std::cout << "\t# mixtures: " << numMixtures << "\n";
	std::cout << "\tmixture definition: " << mixDef << "\n";
	ROCParameter parameter = ROCParameter(genome.getGenomeSize(), sphi_init, numMixtures, geneAssignment, true, mixDef);
	std::string files[] = {std::string("Skluyveri_CSP_ChrA.csv"), std::string("Skluyveri_CSP_ChrCleft.csv")};
	parameter.initMutationSelectionCategories(files, parameter.getNumMutationCategories(), ROCParameter::dM);
	parameter.initMutationSelectionCategories(files, parameter.getNumSelectionCategories(), ROCParameter::dEta);

	//parameter.InitializeExpression(genome, sphi_init);
	double phiVals[genome.getGenomeSize()];
	parameter.readStaticPhiValues("Skluyveri_ChrA_ChrCleft_phi_est.csv", phiVals);
	parameter.InitializeExpression(phiVals);
	std::cout << "done initialize ROCParameter object" << std::endl;

	genome.simulateGenome(parameter, model);

	std::vector <Gene> simGenes = genome.getSimulatedGenes();
	unsigned aaRange[2];

	std::cout <<"FREQUENCIES:\n";
	double phi;
	unsigned mixtureElement;
	unsigned mutationCategory;
	unsigned selectionCategory;
	unsigned expressionCategory;


	for (int i = 0; i < simGenes.size(); i++)
	{
		mixtureElement = parameter.getMixtureAssignment(i);
		mutationCategory = parameter.getMutationCategory(mixtureElement);
		selectionCategory = parameter.getSelectionCategory(mixtureElement);
		expressionCategory = parameter.getExpressionCategory(mixtureElement);
		phi= parameter.getExpression(i, expressionCategory, false);

		std::cout <<"phi = " << phi <<"\n";
		Gene gene = simGenes[i];
		SequenceSummary simSeqSum = gene.geneData;
		for (int j = 0; j < 22; j++)
		{
			SequenceSummary::AAindexToCodonRange(j, false, aaRange);
			char curAA = SequenceSummary::IndexToAA(j);
			unsigned numCodons = simSeqSum.GetNumCodonsForAA(curAA);
			double codonProb[numCodons];
			double mutation[numCodons - 1];
			double selection[numCodons - 1];
			parameter.getParameterForCategory(mutationCategory, ROCParameter::dM, curAA, false, mutation);
			parameter.getParameterForCategory(selectionCategory, ROCParameter::dEta, curAA, false, selection);
			model.calculateCodonProbabilityVector(numCodons, mutation, selection, phi, codonProb);
			std::vector <double> counts(numCodons, 0.0);
			double sum = 0;
			int a = 0;
			for (int k = aaRange[0]; k < aaRange[1]; k++)
			{
				counts[a] = simSeqSum.getCodonCountForCodon(k);
				sum += counts[a];
				a++;
			}
			int aaCount = simSeqSum.getAAcountForAA(j);
			std::cout <<"amino acid " << curAA <<": " << aaCount <<"\n";
			for (int k = 0; k < counts.size(); k++)
			{
				std::cout <<"\tval: " << counts[k]/sum <<" VS " << codonProb[k] <<" with count of " << counts[k] <<"\n";
			}
			std::cout <<"\n";
		}

	}
	genome.writeFasta("fastaTest.fasta", true);
	std::cout << "------------------ TEST SIMULATEGENOME ------------------" << std::endl;
}

int main()
{
	std::cout << "Hello world!" << std::endl << std::endl;

	Genome genome;
	std::cout << "reading fasta file" << std::endl;
	//genome.readFasta("../../inst/testGenome.fasta");
	genome.readFasta("Skluyveri_A_andCleft_simulated.fasta");
	//genome.readFasta("/home/clandere/CodonUsageBias/organisms/yeast/data/LKluyveri/Skluyveri.fasta");
	//genome.writeFasta("../../inst/resGenome.fasta
	//genome.readFasta("testchromosome.fasta");
	std::cout << "done reading fasta file" << std::endl;
	bool testing =  false;

	if(testing)
	{
		/*
			 testNumCodonsPerAA();
			 testCodonRangePerAA(false);
			 testCodonRangePerAA(true);
			 testLogNormDensity();
			 testSCUO(genome);
			 testCovarianceMatrix();
			 testRandMultiNom(3);
			 testThetaKMatrix();
		 */
		testSimulateGenome(genome);
	}else{

		ROCModel model = ROCModel();
		unsigned geneAssignment[genome.getGenomeSize()];
		for(int i = 0; i < genome.getGenomeSize(); i++)
		{
			if(i < 448) geneAssignment[i] = 0u;
			else geneAssignment[i] = 1u;
		}
		std::cout << "initialize ROCParameter object" << std::endl;
		double sphi_init = 2;
		double numMixtures = 2;
		std::string mixDef = ROCParameter::allUnique;
		std::cout << "\tSphi init: " << sphi_init << "\n";
		std::cout << "\t# mixtures: " << numMixtures << "\n";
		std::cout << "\tmixture definition: " << mixDef << "\n";
		ROCParameter parameter = ROCParameter(genome.getGenomeSize(), sphi_init, numMixtures, geneAssignment, true, mixDef);
		std::string files[] = {std::string("Skluyveri_CSP_ChrA.csv"), std::string("Skluyveri_CSP_ChrCleft.csv")};
		parameter.initMutationSelectionCategories(files, parameter.getNumMutationCategories(), ROCParameter::dM);
		parameter.initMutationSelectionCategories(files, parameter.getNumSelectionCategories(), ROCParameter::dEta);
		parameter.InitializeExpression(genome, sphi_init);
		std::cout << "done initialize ROCParameter object" << std::endl;
		//double phiVals[genome.getGenomeSize()];
		//parameter.readStaticPhiValues("Skluyveri_ChrA_phi_est.csv", phiVals);

		//std::cout <<"End of phi Vals\n";
		//parameter.InitializeExpression(phiVals);


		std::cout << "initialize MCMCAlgorithm object" << std::endl;
        int samples = 600;
		int thining = 10;
		int useSamples = 150;
		std::cout << "\t# samples: " << samples << "\n";
		std::cout << "\t thining: " << thining << "\n";
		std::cout << "\t # samples used: " << useSamples << "\n";
		MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thining, true, true, true);
		std::cout << "done initialize MCMCAlgorithm object" << std::endl;


		std::ofstream scuoout("results/scuo.csv");
		for(int n = 0; n < genome.getGenomeSize(); n++)
		{
			scuoout << genome.getGene(n).getId() << "," << parameter.calculateSCUO(genome.getGene(n)) << std::endl;
		}
		scuoout.close();

		std::cout << "starting MCMC" << std::endl;
		mcmc.run(genome, model, parameter);
		std::cout << std::endl << "Finish MCMC" << std::endl;


		std::cout << "Sphi posterior estimate: " << parameter.getSphiPosteriorMean(useSamples) << std::endl;
		std::cout << "Sphi proposal width: " << parameter.getSphiProposalWidth() << std::endl;
		std::cout << "CSP proposal width: \n";
		for(unsigned n = 0; n < 22; n++)
		{
			std::cout << SequenceSummary::AminoAcidArray[n] << ": " << parameter.getCodonSpecificProposalWidth(n) << "\n";
		}

		std::ofstream phiout("results/phiPosterior.csv");
		std::ofstream mixAssignment("results/mixAssignment.csv");
		//std::cout << "Phi proposal width: \n";
		for(int n = 0; n < genome.getGenomeSize(); n++)
		{
			//std::cout << parameter.getExpressionProposalWidth(n) << std::endl;
			double mixtureAssignmentPosteriorMean = parameter.getMixtureAssignmentPosteriorMean(useSamples, n);
			unsigned mixtureElement = std::round( mixtureAssignmentPosteriorMean );
			unsigned expressionCategory = parameter.getExpressionCategory(mixtureElement);
			phiout << genome.getGene(n).getId() << "," << parameter.getExpressionPosteriorMean(useSamples, n, expressionCategory) << std::endl;
			mixAssignment << genome.getGene(n).getId() << "," << mixtureAssignmentPosteriorMean << std::endl;
		}
		phiout.close();
		mixAssignment.close();
	}

	return 0;
}




