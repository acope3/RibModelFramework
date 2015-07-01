#include <iostream>
#include <time.h>
#include <fstream>
#include <algorithm>
#include <sstream>

#include "include/MCMCAlgorithm.h"

void testNumCodonsPerAA()
{
	std::cout << "------------------ # CODONS PER AA ------------------" << std::endl;
	for (int i = 0; i < 22; i++)
	{
		std::string aa = SequenceSummary::AminoAcidArray[i];
		int ncodons = SequenceSummary::GetNumCodonsForAA(aa);
		std::cout << "# codons for " << aa << "\t" << ncodons << std::endl;
	}
	std::cout << "------------------ # CODONS PER AA ------------------" << std::endl;
}

void testCodonRangePerAA(bool forParamVector)
{
	std::cout << "------------------ CODON RANGE PER AA ------------------" << std::endl;
	for (int i = 0; i < 22; i++)
	{
		std::string aa = SequenceSummary::AminoAcidArray[i];
		unsigned codonRange[2];
		SequenceSummary::AAToCodonRange(aa, forParamVector, codonRange);
		std::cout << "codon range for " << aa << "\t" << codonRange[0] << " - " << codonRange[1] << std::endl;
	}
	std::cout << "------------------ CODON RANGE PER AA ------------------" << std::endl;
}

void testLogNormDensity()
{
	std::cout << "------------------ LOG NORM DENSITY ------------------" << std::endl;
	for (int i = 0; i < 5; i++)
	{
		double result = std::log(Parameter::densityLogNorm(i, -1, 1));
		std::cout << "logP of " << i << "\t" << result << std::endl;
	}
	std::cout << "------------------ LOG NORM DENSITY ------------------" << std::endl;
}

void testSCUO(Genome& genome)
{
	std::cout << "------------------ SCUO VALUES ------------------" << std::endl;
	for (unsigned n = 0u; n < genome.getGenomeSize(); n++)
	{
		std::cout << genome.getGene(n).getId() << "\t" << Parameter::calculateSCUO(genome.getGene(n), 22) << std::endl;
	}
	std::cout << "------------------ SCUO VALUES ------------------" << std::endl;
}

void testCovarianceMatrix()
{
	std::cout << "------------------ TEST COVARIANCE ROUTINES ------------------" << std::endl;
	//double arr[3][3] { {4,2,5}, {2,2,3}, {5,3,9} };
	std::vector <double> arr{ 25, 15, -5, 15, 18, 0, -5, 0, 11 };
	CovarianceMatrix covMat = arr;
	std::cout << "Choleski Matrix pre decomposition" << std::endl;
	covMat.printCholeskiMatrix();
	covMat.choleskiDecomposition();
	std::cout << "Covariance Matrix" << std::endl;
	covMat.printCovarianceMatrix();
	std::cout << "Choleski Matrix" << std::endl;
	covMat.printCholeskiMatrix();
	std::cout << "Matrix vector multiplication" << std::endl;
	//double brr[3] {1, 0, -1};
	double res[3];
	//covMat.transformIidNumersIntoCovaryingNumbers(brr, res);
	for (int i = 0; i < 3; i++)
	{
		std::cout << res[i] << std::endl;
	}
	if (res[0] == 5 && res[1] == 3 && res[2] == -4)
	{
		std::cout << "CHECK" << std::endl;
	}
	else{
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
	double assignmentCounts[3] = { 0, 0, 0 };
	double probabilities[3] = { 0.5, 0.2, 0.3 };
	unsigned tmp;
	for (int i = 0; i < 50; i++)
	{
		tmp = Parameter::randMultinom(probabilities, numCat);
		assignmentCounts[tmp] += 1;
	}
	std::cout << "Printing dirichlet numbers:\n";
	for (unsigned i = 0u; i < numCat; i++)
	{
		std::cout << "  " << i << ": " << assignmentCounts[i] << "\n";
		std::cout << "  " << i << "/50: " << assignmentCounts[i] / 50 << "\n\n";
	}

	std::cout << "------------------ TEST RANDMULTINOMIAL ------------------" << std::endl;

}

void testThetaKMatrix()
{
	//unsigned matrix[2][2] = { {2,1}, {1,1} };
	std::cout << "------------------ TEST THETAKMATRIX ------------------" << std::endl;
	std::vector<unsigned> empty(100, 1);
	std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
	ROCParameter R(2, 2, empty, mixtureDefinitionMatrix, true, "allUnique");

	//R.printThetaKMatrix();
	std::cout << "numMutationCategories: " << R.getNumMutationCategories() << "\n";
	std::cout << "numSelectionCategories: " << R.getNumSelectionCategories() << "\n";
	std::vector<std::string> files = { "Skluyveri_CSP_ChrA.csv", "Skluyveri_CSP_ChrCleft.csv" };
	std::cout << "files array good to go\n";
	R.initMutationSelectionCategories(files, R.getNumMutationCategories(), ROCParameter::dM);
	R.initMutationSelectionCategories(files, R.getNumSelectionCategories(), ROCParameter::dEta);


	std::cout << "looping through all mutation params:\n";
	std::vector<std::vector <double>> temp;
	temp = R.getCurrentMutationParameter();
	for (unsigned i = 0u; i < temp.size(); i++)
	{
		std::cout << "Mutation param " << i << ":\n";
		for (unsigned j = 0u; j < temp[i].size(); j++)
		{
			std::cout << temp[i][j] << "\n";
		}
		std::cout << "\n\n";
	}

	std::cout << "looping through all selection params:\n";
	temp = R.getCurrentSelectionParameter();
	for (unsigned i = 0u; i < temp.size(); i++)
	{
		std::cout << "Selection param " << i << ":\n";
		for (unsigned j = 0u; j < temp[i].size(); j++)
		{
			std::cout << temp[i][j] << "\n";
		}
		std::cout << "\n\n";
	}
	std::cout << "------------------ TEST THETAKMATRIX ------------------" << std::endl;

}

void testSimulateGenome(Genome& genome)
{
	std::cout << "------------------ TEST SIMULATEGENOME ------------------" << std::endl;


	std::vector<unsigned> geneAssignment(genome.getGenomeSize());
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		if (i < 448) geneAssignment[i] = 0u;
		else geneAssignment[i] = 1u;
	}



	std::cout << "initialize ROCParameter object" << std::endl;
	double sphi_init = 2;
	unsigned numMixtures = 2;
	std::string mixDef = ROCParameter::mutationShared;
	std::cout << "\tSphi init: " << sphi_init << "\n";
	std::cout << "\t# mixtures: " << numMixtures << "\n";
	std::cout << "\tmixture definition: " << mixDef << "\n";
	std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
	ROCParameter parameter(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);
	std::vector<std::string> files = { "Skluyveri_CSP_ChrA.csv", "Skluyveri_CSP_ChrCleft.csv" };
	parameter.initMutationSelectionCategories(files, parameter.getNumMutationCategories(), ROCParameter::dM);
	parameter.initMutationSelectionCategories(files, parameter.getNumSelectionCategories(), ROCParameter::dEta);



	//parameter.InitializeSynthesisRate(genome, sphi_init);
	std::vector<double> phiVals;
	phiVals = parameter.readPhiValues("Skluyveri_ChrA_ChrCleft_phi_est.csv");
	parameter.InitializeSynthesisRate(phiVals);

	std::cout << "done initialize ROCParameter object" << std::endl;
	ROCModel model;
	model.setParameter(parameter);
	genome.simulateGenome(model);
	std::vector <Gene> simGenes = genome.getSimulatedGenes();
	unsigned aaRange[2];
	std::cout << "FREQUENCIES:\n";


	for (unsigned i = 0u; i < simGenes.size(); i++)
	{
		unsigned mixtureElement = parameter.getMixtureAssignment(i);
		unsigned mutationCategory = parameter.getMutationCategory(mixtureElement);
		unsigned selectionCategory = parameter.getSelectionCategory(mixtureElement);
		unsigned expressionCategory = parameter.getSynthesisRateCategory(mixtureElement);
		double phi = parameter.getSynthesisRate(i, expressionCategory, false);

		std::cout << "phi = " << phi << "\n";
		Gene gene = simGenes[i];
		SequenceSummary simSeqSum = gene.geneData;


		for (int j = 0; j < 22; j++)
		{
			SequenceSummary::AAindexToCodonRange(j, false, aaRange);
			std::string curAA = SequenceSummary::IndexToAA(j);
			unsigned numCodons = simSeqSum.GetNumCodonsForAA(curAA);
			double* codonProb = new double[numCodons];
			double* mutation = new double[numCodons - 1];
			double* selection = new double[numCodons - 1];
			if (curAA == "X") std::cout << numCodons << "\n";
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
			std::cout << "size of counts: " << counts.size() << "\n";
			int aaCount = simSeqSum.getAAcountForAA(j);
			std::cout << "amino acid " << curAA << ": " << aaCount << "\n";
			for (unsigned k = 0u; k < counts.size(); k++)
			{
				std::cout << "\tval: " << counts[k] / sum << " VS " << codonProb[k] << " with count of " << counts[k] << "\n";
			}
			std::cout << "\n";
			delete[] codonProb;
			delete[] mutation;
			delete[] selection;
		}

	}
	std::cout << "Writing new fasta file\n";
	genome.writeFasta("SIMULATED_MUTATION_SHARED_Skluyveri_A_andCleft.fasta", true);
	std::cout << "------------------ TEST SIMULATEGENOME ------------------" << std::endl;
}

void testCovMatrixOverloading()
{
	std::cout << "------------------ TEST COVMATRIXOVERLOADING ------------------" << std::endl;
	CovarianceMatrix cm(3);

	cm.printCovarianceMatrix();
	cm * 4;
	std::cout << "\n";
	cm.printCovarianceMatrix();

	std::cout << "------------------ TEST COVMATRIXOVERLOADING ------------------" << std::endl;
}

void testWriteRestartFile(Genome &genome)
{
	std::cout << "------------------ TEST WRITERESTARTFILE ------------------" << std::endl;
	std::vector<unsigned> geneAssignment(genome.getGenomeSize());
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		if (i < 448) geneAssignment[i] = 0u;
		else geneAssignment[i] = 1u;
	}
	double sphi_init = 2;
	unsigned numMixtures = 2;
	std::string mixDef = ROCParameter::allUnique;
	std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
	ROCParameter parameter = ROCParameter(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);
	std::vector<std::string> files(2);
	files[0] = std::string("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/Skluyveri_CSP_ChrA.csv");
	files[1] = std::string("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/Skluyveri_CSP_ChrCleft.csv");
	parameter.initMutationSelectionCategories(files, parameter.getNumMutationCategories(), ROCParameter::dM);
	parameter.initMutationSelectionCategories(files, parameter.getNumSelectionCategories(), ROCParameter::dEta);
	parameter.InitializeSynthesisRate(genome, sphi_init);
	ROCModel model;
	model.setParameter(parameter);
	model.writeRestartFile("RestartFile1.txt");
	std::cout << "------------------ TEST WRITERESTARTFILE ------------------" << std::endl;
}


void testInitFromRestartFile()
{
	std::cout << "------------------ TEST INITFROMRESTARTFILE ------------------" << std::endl;
	ROCParameter parameter("RestartFile1.txt");
	ROCModel model;
	model.setParameter(parameter);
	model.writeRestartFile("RestartFile2.txt");
	std::cout << "------------------ TEST INITFROMRESTARTFILE ------------------" << std::endl;

}

int main()
{
	enum User { cedric, gabe, jeremy };

	/* Test variables */
	User user = cedric;

	bool read = false;
	unsigned index;

	Genome genome;
	std::cout << "reading fasta file" << std::endl;
	
	switch (user) {
		case cedric:
			genome.readFasta("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/simulatedAllUniqueR.fasta");
			//genome.readFasta("C:/Users/Cedric/Documents/GitHub/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrB_andCleft.fasta");
			break;
		case gabe:
			genome.readFasta("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/simulatedAllUniqueR.fastaa");
			break;
		case jeremy:
			genome.readFasta("C:/Users/Jeremy/Documents/GitHub/RibModelFramework/ribModel/data/simulatedAllUniqueR.fasta");
			break;
	}

	std::cout << "done reading fasta file" << std::endl;
	bool testing = false;

	if (testing)
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
		//testCovMatrixOverloading();
		testWriteRestartFile(genome);
		testInitFromRestartFile();
	}
	else{
		std::cout << "initialize MCMCAlgorithm object" << std::endl;
		int samples = 100;
		int thining = 10;
		int useSamples = 100;
		std::cout << "\t# samples: " << samples << "\n";
		std::cout << "\t thining: " << thining << "\n";
		std::cout << "\t # samples used: " << useSamples << "\n";
		MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thining, 10, true, true, true);
		std::cout << "done initialize MCMCAlgorithm object" << std::endl;
		if (!read)
		{
			std::vector<unsigned> geneAssignment(genome.getGenomeSize());
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 448) geneAssignment[i] = 0u;
				else geneAssignment[i] = 1u;
			}
			std::cout << "initialize ROCParameter object" << std::endl;
			double sphi_init = 2;
			unsigned numMixtures = 2;
			std::string mixDef = ROCParameter::allUnique;
			std::cout << "\tSphi init: " << sphi_init << "\n";
			std::cout << "\t# mixtures: " << numMixtures << "\n";
			std::cout << "\tmixture definition: " << mixDef << "\n";
			std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
			ROCParameter parameter(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			std::vector<std::string> files(2);
			
			switch (user) {
				case cedric:
					files[0] = std::string("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/simulated_CSP0.csv");
					files[1] = std::string("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/simulated_CSP1.csv");
					//files[0] = std::string("C:/Users/Cedric/Documents/GitHub/RibModelFramework/ribModel/data/Skluyveri_CSP_ChrA.csv");
					//files[1] = std::string("C:/Users/Cedric/Documents/GitHub/RibModelFramework/ribModel/data/Skluyveri_CSP_ChrCleft.csv");
					break;
				
				case gabe:
					files[0] = std::string("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/simulated_CSP0.csv");
					files[1] = std::string("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/simulated_CSP1.csv");
					break;
				case jeremy:
					files[0] = std::string("C:/Users/Jeremy/Documents/GitHub/RibModelFramework/ribModel/data/simulated_CSP0.csv");
					files[1] = std::string("C:/Users/Jeremy/Documents/GitHub/RibModelFramework/ribModel/data/simulated_CSP1.csv");
			}
			parameter.initMutationSelectionCategories(files, parameter.getNumMutationCategories(), ROCParameter::dM);
			parameter.initMutationSelectionCategories(files, parameter.getNumSelectionCategories(), ROCParameter::dEta);
			parameter.InitializeSynthesisRate(genome, sphi_init);
			//std::vector<double> phiVals = parameter.readPhiValues("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrCleft_phi_est.csv");
			//parameter.InitializeSynthesisRate(phiVals);
			std::cout << "done initialize ROCParameter object" << std::endl;
			ROCModel model;
			model.setParameter(parameter);
			std::ofstream scuoout("results/scuo.csv");
			for (unsigned n = 0u; n < genome.getGenomeSize(); n++)
			{
				scuoout << genome.getGene(n).getId() << "," << parameter.calculateSCUO(genome.getGene(n), 22) << std::endl;
			}
			scuoout.close();

			mcmc.setRestartFileSettings("RestartFile.txt", 20, true);
			std::cout << "starting MCMC" << std::endl;
			mcmc.run(genome, model);
			std::cout << std::endl << "Finish MCMC" << std::endl;

			std::cout << "Sphi posterior estimate: " << parameter.getSphiPosteriorMean(useSamples) << std::endl;
			std::cout << "Sphi proposal width: " << parameter.getSphiProposalWidth() << std::endl;
			std::cout << "CSP proposal width: \n";
			for (unsigned n = 0; n < model.getGroupListSize(); n++)
			{
				std::string aa = model.getGrouping(n);
				index = SequenceSummary::AAToAAIndex(aa);
				std::cout << SequenceSummary::AminoAcidArray[index] << ": " << parameter.getCurrentCodonSpecificProposalWidth(index) << "\n";
			}
		}
		else
		{
			ROCParameter parameter("RestartFile.txt");
			ROCModel model;
			model.setParameter(parameter);
			std::ofstream scuoout("results/scuo.csv");
			for (unsigned n = 0u; n < genome.getGenomeSize(); n++)
			{
				scuoout << genome.getGene(n).getId() << "," << parameter.calculateSCUO(genome.getGene(n), 22) << std::endl;
			}
			scuoout.close();
			mcmc.setRestartFileSettings("RestartFile.txt", 20, true);
			std::cout << "starting MCMC" << std::endl;
			mcmc.run(genome, model);
			std::cout << std::endl << "Finish MCMC" << std::endl;

			std::cout << "Sphi posterior estimate: " << parameter.getSphiPosteriorMean(useSamples) << std::endl;
			std::cout << "Sphi proposal width: " << parameter.getSphiProposalWidth() << std::endl;
			std::cout << "CSP proposal width: \n";
			for (unsigned n = 0; n < model.getGroupListSize(); n++)
			{
				std::string aa = model.getGrouping(n);
				index = SequenceSummary::AAToAAIndex(aa);
				std::cout << SequenceSummary::AminoAcidArray[index] << ": " << parameter.getCurrentCodonSpecificProposalWidth(index) << "\n";
			}
		}



		//These files used to be written here:
		//mutationPosterior_Cat#.csv
		//selectionPosterior_Cat#.csv
		//expressionPosterior_Cat#.csv
		//mixAssignment.csv
	}
	return 0;
}




