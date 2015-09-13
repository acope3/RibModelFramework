#include <iostream>
#include <time.h>
#include <fstream>
#include <algorithm>
#include <sstream>

#include "include/MCMCAlgorithm.h"


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

/*void testSCUO(Genome& genome)
{
	std::cout << "------------------ SCUO VALUES ------------------" << std::endl;
	for (unsigned n = 0u; n < genome.getGenomeSize(); n++)
	{
		std::cout << genome.getGene(n).getId() << "\t" << Parameter::calculateSCUO(genome.getGene(n), 22) << std::endl;
	}
	std::cout << "------------------ SCUO VALUES ------------------" << std::endl;
}
*/
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
	std::vector<std::string> Mfiles = { "Skluyveri_CSP_ChrA.csv", "Skluyveri_CSP_ChrCleft.csv" };
	std::cout << "files array good to go\n";

	//THIS DOES NOT CURRENTLY WORK!!!!
	R.initMutationCategories(Mfiles, R.getNumMutationCategories());
	R.initSelectionCategories(Mfiles, R.getNumSelectionCategories());


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

void testWriteRestartFile()
{
	Genome genome;
	genome.readFasta("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/fake.fasta", false);
	std::cout << "------------------ TEST WRITERESTARTFILE ------------------" << std::endl;
	std::vector<unsigned> geneAssignment(genome.getGenomeSize());
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		geneAssignment[i] = 0u;
	}
	double sphi_init = 2;
	unsigned numMixtures = 1;
	std::string mixDef = ROCParameter::allUnique;
	std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
	ROCParameter parameter(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);
	std::vector<std::string> files(2);
	files[0] = std::string("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/simulated_mutation0.csv");
	files[1] = std::string("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/simulated_mutation1.csv");
	parameter.initMutationCategories(files, parameter.getNumMutationCategories());

	files[0] = std::string("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/simulated_selection0.csv");
	files[1] = std::string("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/simulated_selection1.csv");
	parameter.initSelectionCategories(files, parameter.getNumSelectionCategories());
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

void simulateRFPData()
{
	Genome genome;
	genome.readRFPFile("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/rfp.counts.by.codon.and.gene.GSE63789.wt.csv");
	std::vector<unsigned> geneAssignment(genome.getGenomeSize());
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		geneAssignment[i] = 0u;
	}
	double sphi_init = 2;
	unsigned numMixtures = 1;
	std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
	RFPParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, "allUnique");

	std::vector<std::string> files;
	files.push_back("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/RFPAlphaValues.csv");
	tmp.initMutationSelectionCategories(files, 1, RFPParameter::alp);
	files[0] = "/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/RFPLambdaPrimeValues.csv";
	tmp.initMutationSelectionCategories(files, 1, RFPParameter::lmPri);
	std::vector<double> phi = tmp.readPhiValues("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/RFPPhiValues.csv");
	tmp.InitializeSynthesisRate(phi);


	RFPModel model;

	model.setParameter(tmp);

	std::cout <<"init done\n";
	model.simulateGenome(genome);
	std::cout <<"writing file\n";
	genome.writeRFPFile("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/SimulatedRFPData.csv", true);


}

void simulateROCData()
{
	Genome genome;
	//genome.readFasta("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/Skluyveri_chromosomeA.fasta");
	genome.readFasta("C:/Users/Jeremy/Documents/GitHub/RibModelFramework/ribModel/data/Skluyveri_chromosomeA.fasta");

	std::vector<unsigned> geneAssignment(genome.getGenomeSize());
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
    {
        if (i < 448) geneAssignment[i] = 0u;
        else geneAssignment[i] = 1u;
    }

	double sphi_init = 2;
	unsigned numMixtures = 2;
	std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;

	ROCParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, Parameter::allUnique);
	std::vector<std::string> files = { "/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/Skluyveri_mutation_ChrA.csv",
									   "/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/Skluyveri_mutation_ChrCleft.csv" };
	tmp.initMutationCategories(files, tmp.getNumMutationCategories());
	files[0] = "/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/Skluyveri_selection_ChrA.csv";
	files[1] = "/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/Skluyveri_selection_ChrCleft.csv";
	tmp.initSelectionCategories(files, tmp.getNumSelectionCategories());
	ROCModel model(false);

	model.setParameter(tmp);

	std::cout <<"init done\n";

	model.simulateGenome(genome);

	//genome.writeFasta("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/SimulatedROCData.fasta", true);
	genome.writeFasta("C:/Users/Jeremy/Documents/GitHub/RibModelFramework/ribModel/data/SimulatedROCData.fasta", true);
}




void testInitMutationSelection()
{
	Genome genome;
	genome.readRFPFile("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/rfp.counts.by.codon.and.gene.GSE63789.wt.csv");
	std::vector<unsigned> geneAssignment(genome.getGenomeSize());
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		geneAssignment[i] = 0u;
	}
	double sphi_init = 2;
	unsigned numMixtures = 1;
	std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
	RFPParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, "allUnique");

	std::vector<std::string> files;
	files.push_back("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/RFPAlphaValues.csv");
	tmp.initMutationSelectionCategories(files, 1, RFPParameter::alp);


}

void testRFPVarianceAndMean()
{
	Genome genome;
	genome.readRFPFile("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/rfp.counts.by.codon.and.gene.GSE63789.wt.csv");
	std::vector<unsigned> geneAssignment(genome.getGenomeSize());
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		geneAssignment[i] = 0u;
	}
	double sphi_init = 2;
	unsigned numMixtures = 1;
	std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
	RFPParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, "allUnique");

	std::cout <<"about to calculate\n";
	tmp.calculateRFPMean(genome);

	Genome genome2;
	genome2.readRFPFile("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/SimulatedRFPData.csv");
	std::vector<unsigned> geneAssignment2(genome2.getGenomeSize());
	for (unsigned i = 0u; i < genome2.getGenomeSize(); i++)
	{
		geneAssignment2[i] = 0u;
	}
	RFPParameter tmp2(sphi_init, numMixtures, geneAssignment2, mixtureDefinitionMatrix, true, "allUnique");

	tmp2.calculateRFPMean(genome2);


}


void testReadMutationValues()
{
	std::cout <<"----------TEST INITMUTATIONCATEGORIES----------\n";

	Genome genome;
	genome.readFasta("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/Skluyveri_chromosomeA.fasta");


	std::vector<unsigned> geneAssignment(genome.getGenomeSize());
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		if (i < 448) geneAssignment[i] = 0u;
		else geneAssignment[i] = 1u;
	}

	double sphi_init = 2;
	unsigned numMixtures = 2;
	std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;

	ROCParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, Parameter::allUnique);
//	std::vector<std::string> files = { "/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/new_simulated_mutation0.csv",
	//								   "/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/new_simulated_mutation1.csv" };
	std::vector<std::string> files = { "C:/Users/Jeremy/Documents/GitHub/RibModelFramework/ribModel/data/new_simulated_mutation0.csv",
									   "C:/Users/Jeremy/Documents/GitHub/RibModelFramework/ribModel/data/new_simulated_mutation1.csv" };


	tmp.initMutationCategories(files, 2);

	std::vector <std::vector <double>> v = tmp.getCurrentMutationParameter();

	for (unsigned i = 0; i < v.size(); i++)
	{
		for(unsigned j = 0; j < v[i].size(); j++)
		{
			std::cout << v[i][j] <<"\n";
		}
		std::cout <<"\n\n\n";
	}


	std::cout <<"----------END TEST INITMUTATIONCATEGORIES----------\n";
}

void testMultiplePhi()
{
	std::cout <<"----------TEST MULTIPLEPHI----------\n";

	Genome genome;
	genome.readFasta("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/test.fasta");
	genome.readObservedPhiValues("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/testReadObservedPhiValuesValid.csv", false);

	for (unsigned i = 0; i < genome.getGenomeSize(); i++)
	{
		Gene gene = genome.getGene(i);
		std::cout << gene.getId() <<": ";
		for (unsigned j = 0; j < gene.observedPhiValues.size(); j++)
		{
			std::cout << gene.observedPhiValues[j] <<", ";
		}
		std::cout <<"\n";
	}


	std::cout <<"----------END TEST MULTIPLEPHI----------\n";
}


void codonTableTest()
{
	CodonTable CT(1, false);
	CT.setupCodonTable();
	CT.AAToCodon("C", false);

	return;
}

void testProcessSequence()
{
	std::cout <<"Testing process sequence\n";
	CodonTable::createCodonTable(1, true);
	SequenceSummary SS("ATGCTCATTCTCAGTGCTGCCTCGTAG");
}
//---------------------------------------------------------------------------------------------//
//-----------------------------------------END OF UNIT TESTING---------------------------------//
//---------------------------------------------------------------------------------------------//

int main()
{
	enum User { cedric, gabe, jeremy }; //For file IO paths
	enum ModelToRun { ROC, RFP, FONSE };


	User user = gabe;
	ModelToRun modelToRun = ROC;

	bool read = false; //Initialize parameter by reading a restart file
	bool testing = false; //Run unit testing
	bool withPhi = false; //Using multiple phi values per gene


	if (testing) //Unit testing
	{
		//testLogNormDensity();
		//testSCUO(genome);
		//testCovarianceMatrix();
		//testRandMultiNom(3);
		//testThetaKMatrix();
		//testCovMatrixOverloading();
		//testWriteRestartFile();
		//testInitFromRestartFile();
		//simulateRFPData();
		//simulateROCData();
		//testInitMutationSelection();
		//testRFPVarianceAndMean();
		//testReadMutationValues();
		//testGeneSequenceSummary();
		//testReadMutationValues();
		//testMultiplePhi();
		//codonTableTest();
		testProcessSequence();
	}
	else //Not doing unit testing, running a model
	{
		unsigned index;
		std::cout << "initialize MCMCAlgorithm object" << std::endl;
		int samples = 100;
		int thining = 10;
		int useSamples = 100;

		std::cout << "Initializing MCMCAlgorithm object" << std::endl;
		std::cout << "\t# samples: " << samples << "\n";
		std::cout << "\t thining: " << thining << "\n";
		std::cout << "\t # samples used: " << useSamples << "\n";


		MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thining, 10, true, true, true);
		mcmc.setRestartFileSettings("RestartFile.txt", 20, true);
		std::cout << "Done initializing MCMCAlgorithm object" <<"\n\n";
		

		std::cout << "Initializing Genome object" << std::endl;
		unsigned tableId = 1;
		bool splitAA = false;
		Genome genome(tableId, splitAA);


		switch (user) {
			case cedric:
				if (modelToRun == ROC || modelToRun == FONSE)
				{
					genome.readFasta("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/simulatedAllUniqueR.fasta");
					//genome.readFasta("C:/Users/Cedric/Documents/GitHub/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrB_andCleft.fasta");
				}
				else if (modelToRun == RFP)
				{
					genome.readRFPFile("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/rfp.counts.by.codon.and.gene.GSE63789.wt.csv");
				}
				else {}
				break;

			case gabe:
				if (modelToRun == ROC || modelToRun == FONSE)
				{
					std::cout <<"\tReading Fasta File\n";
					genome.readFasta("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/simulatedAllUniqueR.fasta");
					if (withPhi)
					{
						genome.readObservedPhiValues("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/simulatedAllUniqueR_phi.csv", false);
					}
				}
				else if (modelToRun == RFP)
				{
					genome.readRFPFile("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/rfp.counts.by.codon.and.gene.GSE63789.wt.csv");
				}
				else {}
				break;

			case jeremy:
				if (modelToRun == ROC || modelToRun == FONSE)
				{
					genome.readFasta("C:/Users/Jeremy/Documents/GitHub/RibModelFramework/ribModel/data/simulatedAllUniqueR.fasta");
					if (withPhi) {
						genome.readObservedPhiValues("C:/Users/Jeremy/Documents/GitHub/RibModelFramework/ribModel/data/simulatedAllUniqueR_phi.csv", false);
					}
				}
				else if (modelToRun == RFP)
				{
					genome.readRFPFile("C:/Users/Jeremy/Documents/GitHub/RibModelFramework/ribModel/data/rfp.counts.by.codon.and.gene.GSE63789.wt.csv");
				}
				else {}
				break;
		}
		std::cout << "Done initializing Genome object" <<"\n\n";


		std::cout << "Initializing shared parameter variables\n";
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());


		/* For 2 mixtures */
		for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
		{
			if (i < 500) geneAssignment[i] = 0u;
			else geneAssignment[i] = 1u;
		}

		/* For 1 mixture */
		/*for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
		{
			geneAssignment[i] = 0u;
		}*/
		double sphi_init = 2;
		unsigned numMixtures = 2;
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		std::cout <<"Done initializing shared parameter variables\n";


		if (modelToRun == ROC)
		{

			ROCParameter parameter;
			if (read)
			{
				std::cout <<"Initializing ROCParameter object from Restart File\n";
				ROCParameter tmp("RestartFile.txt");
				parameter = tmp;
				std::cout <<"Done initializing ROCParameter object from Restart File\n";
			}
			else
			{
			std::cout << "initialize ROCParameter object" << std::endl;
			std::string mixDef = ROCParameter::allUnique;
			std::cout << "\tSphi init: " << sphi_init << "\n";
			std::cout << "\t# mixtures: " << numMixtures << "\n";
			std::cout << "\tmixture definition: " << mixDef << "\n";

			ROCParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);
			std::vector<std::string> files(2);

			switch (user) {
				case cedric:
					files[0] = std::string("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/simulated_mutation0.csv");
					files[1] = std::string("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/simulated_mutation1.csv");
					tmp.initMutationCategories(files, tmp.getNumMutationCategories());

					files[0] = std::string("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/simulated_selection0.csv");
					files[1] = std::string("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/simulated_selection1.csv");
					tmp.initSelectionCategories(files, tmp.getNumSelectionCategories());
					break;

				case gabe:
					files[0] = std::string("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/simulated_mutation0.csv");
					files[1] = std::string("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/simulated_mutation1.csv");
				//	tmp.initMutationCategories(files, tmp.getNumMutationCategories());

					files[0] = std::string("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/simulated_selection0.csv");
					files[1] = std::string("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/simulated_selection1.csv");
				//	tmp.initSelectionCategories(files, tmp.getNumSelectionCategories());
					break;
				case jeremy:
					files[0] = std::string("C:/Users/Jeremy/Documents/GitHub/RibModelFramework/ribModel/data/simulated_mutation0.csv");
					files[1] = std::string("C:/Users/Jeremy/Documents/GitHub/RibModelFramework/ribModel/data/simulated_mutation1.csv");
				//	tmp.initMutationCategories(files, tmp.getNumMutationCategories());

					files[0] = std::string("C:/Users/Jeremy/Documents/GitHub/RibModelFramework/ribModel/data/simulated_selection0.csv");
					files[1] = std::string("C:/Users/Jeremy/Documents/GitHub/RibModelFramework/ribModel/data/simulated_selection1.csv");
				//	tmp.initSelectionCategories(files, tmp.getNumSelectionCategories());
					break;
			}
			tmp.InitializeSynthesisRate(genome, sphi_init);
			//std::vector<double> phiVals = parameter.readPhiValues("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrCleft_phi_est.csv");
			//parameter.InitializeSynthesisRate(phiVals);
			parameter = tmp;
			std::cout << "Done initializing ROCParameter object\n\n";
			}

			std::cout <<"Initializing ROCModel object\n";
			ROCModel model(withPhi);
			model.setParameter(parameter);
			std::ofstream scuoout("results/scuo.csv");
			for (unsigned n = 0u; n < genome.getGenomeSize(); n++)
			{
				scuoout << genome.getGene(n).getId() << "," << parameter.calculateSCUO(genome.getGene(n)) << std::endl;
			}
			scuoout.close();
			std::cout <<"Done initializing ROCModel object\n\n";

			std::cout <<"Temparary exit\n";
			

			std::cout << "starting MCMC for ROC" << std::endl;
			mcmc.run(genome, model, 1, 20);
			std::cout << std::endl << "Finished MCMC for ROC" << std::endl;
			

			std::cout << "Sphi posterior estimate: " << parameter.getSphiPosteriorMean(useSamples) << std::endl;
			std::cout << "Sphi proposal width: " << parameter.getCurrentSphiProposalWidth() << std::endl;
			std::cout << "CSP proposal width: \n";
			CodonTable *codonTable = CodonTable::getInstance();
			for (unsigned n = 0; n < model.getGroupListSize(); n++)
			{
				std::string aa = model.getGrouping(n);
				index = codonTable -> AAToAAIndex(aa);
				std::cout << CodonTable::aminoAcidArray[index] << ": " << parameter.getCurrentCodonSpecificProposalWidth(index) << "\n";
			}
		}
		else if (modelToRun == RFP)
		{
			RFPParameter parameter;
			if (read)
			{
				std::cout <<"Initializing RFPParameter object from Restart File\n";
				RFPParameter tmp("RestartFile1.txt");
				parameter = tmp;
				std::cout <<"Done initializing RFPParameter object from Restart File\n";
			}
			else
			{
				std::cout << "initialize RFPParameter object" << std::endl;
				std::string mixDef = Parameter::allUnique;
				std::cout << "\tSphi init: " << sphi_init << "\n";
				std::cout << "\t# mixtures: " << numMixtures << "\n";
				std::cout << "\tmixture definition: " << mixDef << "\n";
				RFPParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);
				tmp.InitializeSynthesisRate(genome, sphi_init);
				parameter = tmp;
				std::cout << "Done initialize RFPParameter object" << std::endl;
			}
			

			std::cout <<"Initializing RFPModel object\n";
			RFPModel model;
			model.setParameter(parameter);
			//TODO: what does this do in the RFP case? Do we even need it????
		/*	std::ofstream scuoout("results/scuo.csv");
			for (unsigned n = 0u; n < genome.getGenomeSize(); n++)
			{
				scuoout << genome.getGene(n).getId() << "," << parameter.calculateSCUO(genome.getGene(n), 64) << std::endl;
			}
			scuoout.close();
			*/std::cout <<"Done initializing RFPModel object\n";
		
			std::cout << "starting MCMC for RFP" << std::endl;
      mcmc.run(genome, model, 8);
      std::cout << std::endl << "Finished MCMC for RFP" << std::endl;


      std::cout << "Sphi posterior estimate: " << parameter.getSphiPosteriorMean(useSamples) << std::endl;
      std::cout << "Sphi proposal width: " << parameter.getCurrentSphiProposalWidth() << std::endl;

		}
		else if (modelToRun == FONSE)
		{
			FONSEParameter parameter;
			if (read)
			{
				std::cout << "Initializing ROCParameter object from Restart File\n";
				FONSEParameter tmp("RestartFile.txt");
				parameter = tmp;
				std::cout << "Done initializing ROCParameter object from Restart File\n";
			}
			else
			{
				std::cout << "initialize ROCParameter object" << std::endl;
				std::string mixDef = FONSEParameter::allUnique;
				std::cout << "\tSphi init: " << sphi_init << "\n";
				std::cout << "\t# mixtures: " << numMixtures << "\n";
				std::cout << "\tmixture definition: " << mixDef << "\n";

				FONSEParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);
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
					break;
				}

				tmp.initMutationSelectionCategories(files, tmp.getNumMutationCategories(), FONSEParameter::dM);
				tmp.initMutationSelectionCategories(files, tmp.getNumSelectionCategories(), FONSEParameter::dOmega);
				tmp.InitializeSynthesisRate(genome, sphi_init);
				//std::vector<double> phiVals = parameter.readPhiValues("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrCleft_phi_est.csv");
				//parameter.InitializeSynthesisRate(phiVals);
				parameter = tmp;
				std::cout << "done initialize FONSEParameter object" << std::endl;
			}

			std::cout << "Initializing FONSEModel object\n";
			FONSEModel model;
			model.setParameter(parameter);
			std::ofstream scuoout("results/scuo.csv");
			for (unsigned n = 0u; n < genome.getGenomeSize(); n++)
			{
				scuoout << genome.getGene(n).getId() << "," << parameter.calculateSCUO(genome.getGene(n)) << std::endl;
			}
			scuoout.close();
			std::cout << "Done initializing FONSEModel object\n";


			std::cout << "starting MCMC for FONSE" << std::endl;
			mcmc.run(genome, model, 4);
			std::cout << std::endl << "Finished MCMC for FONSE" << std::endl;


			std::cout << "Sphi posterior estimate: " << parameter.getSphiPosteriorMean(useSamples) << std::endl;
			std::cout << "Sphi proposal width: " << parameter.getCurrentSphiProposalWidth() << std::endl;
			std::cout << "CSP proposal width: \n";
			CodonTable *codonTable = CodonTable::getInstance();
			for (unsigned n = 0; n < model.getGroupListSize(); n++)
			{
				std::string aa = model.getGrouping(n);
				index = codonTable -> AAToAAIndex(aa);
				std::cout << CodonTable::aminoAcidArray[index] << ": " << parameter.getCurrentCodonSpecificProposalWidth(index) << "\n";
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
