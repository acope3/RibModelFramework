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
		unsigned* codonRange = SequenceSummary::AAToCodonRange(aa, forParamVector);
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

void testRandMultiNom(unsigned numCat)
{
	std::cout << "------------------ TEST RANDMULTINOMIAL ------------------" << std::endl;
	double dirichletParameters[3] = {0, 0, 0};
	double probabilities[3] = {0.5, 0.2, 0.3};
	unsigned tmp;
	for (int i = 0; i < 50; i++)
	{
		tmp = ROCParameter::randMultinom(probabilities, numCat);
		dirichletParameters[tmp] += 1;
	}
	std::cout <<"Printing dirichlet numbers:\n";
	for (int i = 0; i < numCat; i++)
	{
		std::cout <<"  " << i <<": " << dirichletParameters[i] <<"\n";
		std::cout <<"  " << i <<"/50: " << dirichletParameters[i]/50 <<"\n\n";
	}

	std::cout << "------------------ TEST RANDMULTINOMIAL ------------------" << std::endl;

}//In progress

int main()
{
	std::cout << "Hello world!" << std::endl << std::endl;


	Genome genome;
	std::cout << "reading fasta file" << std::endl;
	//genome.readFasta("../../inst/testGenome.fasta");
	genome.readFasta("../Skluyveri_chromosomeA.fasta");
	//genome.readFasta("/home/clandere/CodonUsageBias/organisms/yeast/data/LKluyveri/Skluyveri.fasta");
	//genome.writeFasta("../../inst/resGenome.fasta
	bool testing = false;

	if(testing)
	{
		/*
			 testNumCodonsPerAA();
			 testCodonRangePerAA(false);
			 testCodonRangePerAA(true);
			 testLogNormDensity();
			 testSCUO(genome);
			 testCovarianceMatrix();*/
		testRandMultiNom(3);

	}else{
		int samples = 1000;
		int thining = 10;
		int useSamples = 150;

		ROCModel model = ROCModel();
		ROCParameter parameter = ROCParameter(1, 1, genome.getGenomeSize(), .8, 1, false);
		parameter.InitializeExpression(genome, .8);
		MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thining, true, false, true);

		std::ofstream scuoout("/results/test.scuo");
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

		std::ofstream phiout("../results/test.phi");
		for(int n = 0; n < genome.getGenomeSize(); n++)
		{
			phiout << genome.getGene(n).getId() << "," << parameter.getExpressionPosteriorMean(useSamples, n, 0) << std::endl;
		}
		phiout.close();
	}


	return 0;
}




