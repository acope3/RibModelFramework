#include <iostream>
#include <time.h>
#include <fstream>
#include <algorithm>
#include <sstream>

#include "include/MCMCAlgorithm.h"
#include "include/CovarianceMatrix.h"


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
	for(unsigned n = 0u; n < genome.getGenomeSize(); n++)
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
	// TODO I changed the function, have to reimplement test function
	//getCountsForAA
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
	for (unsigned i = 0u; i < numCat; i++)
	{
		std::cout <<"  " << i <<": " << assignmentCounts[i] <<"\n";
		std::cout <<"  " << i <<"/50: " << assignmentCounts[i]/50 <<"\n\n";
	}

	std::cout << "------------------ TEST RANDMULTINOMIAL ------------------" << std::endl;

}

void testThetaKMatrix()
{
	//unsigned matrix[2][2] = { {2,1}, {1,1} };
	std::cout << "------------------ TEST THETAKMATRIX ------------------" << std::endl;
	std::vector<unsigned> empty(100, 1);
	std::vector<std::vector<unsigned>> thetaKMatrix;
	ROCParameter R(2, 2, empty, thetaKMatrix, true, "allUnique");

	R.printThetaKMatrix();
	std::cout <<"numMutationCategories: " << R.getNumMutationCategories() <<"\n";
	std::cout <<"numSelectionCategories: " << R.getNumSelectionCategories() <<"\n";
	std::vector<std::string> files = {"Skluyveri_CSP_ChrA.csv", "Skluyveri_CSP_ChrCleft.csv"};
	std::cout <<"files array good to go\n";
	R.initMutationSelectionCategories(files, R.getNumMutationCategories(), ROCParameter::dM);
	R.initMutationSelectionCategories(files, R.getNumSelectionCategories(), ROCParameter::dEta);


	std::cout <<"looping through all mutation params:\n";
	std::vector<std::vector <double>> temp;
	temp = R.getCurrentMutationParameter();
	for (unsigned i = 0u; i < temp.size(); i++)
	{
		std::cout <<"Mutation param " << i <<":\n";
		for (unsigned j = 0u; j < temp[i].size(); j++)
		{
			std::cout << temp[i][j] <<"\n";
		}
		std::cout <<"\n\n";
	}

	std::cout <<"looping through all selection params:\n";
	temp = R.getCurrentSelectionParameter();
	for (unsigned i = 0u; i < temp.size(); i++)
	{
		std::cout <<"Selection param " << i <<":\n";
		for (unsigned j = 0u; j < temp[i].size(); j++)
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
	std::vector<unsigned> geneAssignment(genome.getGenomeSize());
	for(unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		if(i < 448) geneAssignment[i] = 0u;
		else geneAssignment[i] = 1u;
	}



	std::cout << "initialize ROCParameter object" << std::endl;
	double sphi_init = 2;
	unsigned numMixtures = 2;
	std::string mixDef = ROCParameter::mutationShared;
	std::cout << "\tSphi init: " << sphi_init << "\n";
	std::cout << "\t# mixtures: " << numMixtures << "\n";
	std::cout << "\tmixture definition: " << mixDef << "\n";
	std::vector<std::vector<unsigned>> thetaKMatrix;
	ROCParameter parameter(sphi_init, numMixtures, geneAssignment, thetaKMatrix, true, mixDef);
	std::vector<std::string> files = {"Skluyveri_CSP_ChrA.csv", "Skluyveri_CSP_ChrCleft.csv"};
	parameter.initMutationSelectionCategories(files, parameter.getNumMutationCategories(), ROCParameter::dM);
	parameter.initMutationSelectionCategories(files, parameter.getNumSelectionCategories(), ROCParameter::dEta);



	//parameter.InitializeExpression(genome, sphi_init);
	std::vector<double> phiVals;
	phiVals = parameter.readPhiValues("Skluyveri_ChrA_ChrCleft_phi_est.csv");
	parameter.InitializeExpression(phiVals);

	std::cout << "done initialize ROCParameter object" << std::endl;

	genome.simulateGenome(parameter, model);
	std::vector <Gene> simGenes = genome.getSimulatedGenes();
	unsigned aaRange[2];
	std::cout <<"FREQUENCIES:\n";


	for (unsigned i = 0u; i < simGenes.size(); i++)
	{
		unsigned mixtureElement = parameter.getMixtureAssignment(i);
		unsigned mutationCategory = parameter.getMutationCategory(mixtureElement);
		unsigned selectionCategory = parameter.getSelectionCategory(mixtureElement);
		unsigned expressionCategory = parameter.getExpressionCategory(mixtureElement);
		double phi= parameter.getExpression(i, expressionCategory, false);

		std::cout <<"phi = " << phi <<"\n";
		Gene gene = simGenes[i];
		SequenceSummary simSeqSum = gene.geneData;


		for (int j = 0; j < 22; j++)
		{
			SequenceSummary::AAindexToCodonRange(j, false, aaRange);
			char curAA = SequenceSummary::IndexToAA(j);
			unsigned numCodons = simSeqSum.GetNumCodonsForAA(curAA);
			double* codonProb = new double[numCodons];
			double* mutation = new double[numCodons - 1];
			double* selection = new double[numCodons - 1];
			if (curAA == 'X') std::cout << numCodons <<"\n";
			parameter.getParameterForCategory(mutationCategory, ROCParameter::dM, curAA, false, mutation);
			parameter.getParameterForCategory(selectionCategory, ROCParameter::dEta, curAA, false, selection);
			model.calculateCodonProbabilityVector(numCodons, mutation, selection, phi, codonProb);
			std::vector <double> counts(numCodons, 0.0);
			double sum = 0;
			int a = 0;
			for (unsigned k = aaRange[0]; k < aaRange[1]; k++)
			{
				counts[a] = simSeqSum.getCodonCountForCodon(k);
				sum += counts[a];
				a++;
			}
			std::cout <<"size of counts: " << counts.size() <<"\n";
			int aaCount = simSeqSum.getAAcountForAA(j);
			std::cout <<"amino acid " << curAA <<": " << aaCount <<"\n";
			for (unsigned k = 0u; k < counts.size(); k++)
			{
				std::cout <<"\tval: " << counts[k]/sum <<" VS " << codonProb[k] <<" with count of " << counts[k] <<"\n";
			}
			std::cout <<"\n";
			delete [] codonProb;
			delete [] mutation;
			delete [] selection;
		}

	}
	std::cout <<"Writing new fasta file\n";
	genome.writeFasta("SIMULATED_MUTATION_SHARED_Skluyveri_A_andCleft.fasta", true);
	std::cout << "------------------ TEST SIMULATEGENOME ------------------" << std::endl;
}

void testCovMatrixOverloading()
{
	std::cout << "------------------ TEST COVMATRIXOVERLOADING ------------------" << std::endl;
	CovarianceMatrix cm(3);

	cm.printCovarianceMatrix();
	cm * 4;
	std::cout <<"\n";
	cm.printCovarianceMatrix();

	std::cout << "------------------ TEST COVMATRIXOVERLOADING ------------------" << std::endl;
}

int main()
{
	bool cedric = true;
	std::cout << "Hello world!" << std::endl << std::endl;

	Genome genome;
	std::cout << "reading fasta file" << std::endl;
	if(cedric){
		genome.readFasta("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrB_andCleft.fasta");
		//genome.readFasta("C:/Users/Cedric/Documents/GitHub/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrB_andCleft.fasta");
	}else{
		genome.readFasta("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/Skluyveri_ChrA_andCleft.fasta");
	}
	std::cout << "done reading fasta file" << std::endl;
	bool testing =  false;

	if(testing)
	{
			 //testNumCodonsPerAA();
			 //testCodonRangePerAA(false);
			 //testCodonRangePerAA(true);
			 //testLogNormDensity();
			 //testSCUO(genome);
			 //testCovarianceMatrix();
			 //testRandMultiNom(3);
			 //testThetaKMatrix();
		//testSimulateGenome(genome);
		testCovMatrixOverloading();
	}else{
		ROCModel model;
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		for(unsigned i = 0u; i < genome.getGenomeSize(); i++)
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
		std::vector<std::vector<unsigned>> thetaKMatrix;
		ROCParameter parameter = ROCParameter(sphi_init, numMixtures, geneAssignment, thetaKMatrix, true, mixDef);

		std::vector<std::string> files(2);
		if(cedric)
		{
			files[0] = std::string("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_CSP_ChrA.csv");
			files[1] = std::string("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_CSP_ChrCleft.csv");
		}else{
			files[0] = std::string("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/Skluyveri_CSP_ChrA.csv");
			files[1] = std::string("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/Skluyveri_CSP_ChrCleft.csv");
		}
		parameter.initMutationSelectionCategories(files, parameter.getNumMutationCategories(), ROCParameter::dM);
		parameter.initMutationSelectionCategories(files, parameter.getNumSelectionCategories(), ROCParameter::dEta);
		parameter.InitializeExpression(genome, sphi_init);
		//std::vector<double> phiVals = parameter.readPhiValues("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrCleft_phi_est.csv");
		//parameter.InitializeExpression(phiVals);
		std::cout << "done initialize ROCParameter object" << std::endl;

		std::cout << "initialize MCMCAlgorithm object" << std::endl;
        int samples = 100;
		int thining = 10;
		int useSamples = 100;

		std::cout << "\t# samples: " << samples << "\n";
		std::cout << "\t thining: " << thining << "\n";
		std::cout << "\t # samples used: " << useSamples << "\n";
		MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thining, 10, true, true, true);
		std::cout << "done initialize MCMCAlgorithm object" << std::endl;

		std::ofstream scuoout("results/scuo.csv");
		for(unsigned n = 0u; n < genome.getGenomeSize(); n++)
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
            if(n == 21 || n == 10 || n == 18) continue;
			std::cout << SequenceSummary::AminoAcidArray[n] << ": " << parameter.getCodonSpecificProposalWidth(n) << "\n";
		}

        std::cout << "writing mutation posterior file" << "\n";
		//Get posterior estimates for mutation & selection
		int numMutationCategories = parameter.getNumMutationCategories();
		int numSelectionCategories = parameter.getNumSelectionCategories();
		for (int i = 0; i < numMutationCategories; i++)
		{
			std::ostringstream strstream;
			strstream << i;
			std::string file = "results/mutationPosterior_Cat" + strstream.str() + ".csv";
			std::ofstream mutout(file);
			for (int n = 0; n < 22; n++) //going over the amino acids
			{
				unsigned aaRange[2];
				char aa = SequenceSummary::AminoAcidArray[n];
				if (aa == 'X' || aa == 'M' || aa == 'W') continue;
				SequenceSummary::AAToCodonRange(aa, true, aaRange);
				for (unsigned a = aaRange[0]; a <= aaRange[1]; a++)
				{
					double estimate = 0.0;
					double variance = 0.0;
					if (a != aaRange[1])
					{
						estimate = parameter.getMutationPosteriorMean(i, useSamples, a);
						variance = parameter.getMutationVariance(i, useSamples, a);
					}
					std::string codon = SequenceSummary::IndexToCodon(a);
					mutout << aa <<"." << codon <<".deltaM," << estimate <<"," << variance <<"\n";
				}
			}
			mutout.close();
		}
        std::cout << "finished writing mutation posterior file" << "\n";
        std::cout << "writing selection posterior file" << "\n";
		for (int i = 0; i < numSelectionCategories; i++)
		{
			std::ostringstream strstream;
			strstream << i;
			std::string file = "results/selectionPosterior_Cat" + strstream.str() + ".csv";
			std::ofstream selectout(file);
			for (int n = 0; n < 22; n++)
			{
				unsigned aaRange[2];
				char aa = SequenceSummary::AminoAcidArray[n];
				if (aa == 'X' || aa == 'M' || aa == 'W') continue;
				SequenceSummary::AAToCodonRange(aa, true, aaRange);
				for (unsigned a = aaRange[0]; a <= aaRange[1]; a++)
				{
					std::string codon = SequenceSummary::IndexToCodon(a);
					double estimate = 0.0;
					double variance = 0.0;
					if (a != aaRange[1])
					{
						estimate = parameter.getSelectionPosteriorMean(i, useSamples, a);
						variance = parameter.getSelectionVariance(i, useSamples, a);
					}
					selectout << aa << "." << codon <<".deltaEta," << estimate <<"," << variance <<"\n";
				}
			}
			selectout.close();
		}
        std::cout << "finished writing selection posterior file" << "\n";



        for(unsigned k = 0u; k < parameter.getNumExpressionCategories(); k++)
        {
			std::ostringstream strstream;
			strstream << k;
            std::string file = "results/expressionPosterior_Cat" + strstream.str() + ".csv";
            std::ofstream phiout(file);
            for(unsigned n = 0u; n < genome.getGenomeSize(); n++)
            {
                phiout << genome.getGene(n).getId() << "," << parameter.getExpressionPosteriorMean(useSamples, n, k) << "," << parameter.getExpressionVariance(useSamples, n, k) << std::endl;
            }
            phiout.close();
        }

		std::ofstream mixAssignment("results/mixAssignment.csv");
		for(unsigned n = 0u; n < genome.getGenomeSize(); n++)
		{
			unsigned mixtureAssignment = parameter.getEstimatedMixtureAssignment(useSamples, n);
			mixAssignment << genome.getGene(n).getId() << "," << mixtureAssignment << std::endl;
		}
		mixAssignment.close();
	}

	return 0;
}




