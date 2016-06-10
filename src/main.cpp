#include "include/MCMCAlgorithm.h"
#include "include/Testing.h"
#include <vector>

#ifdef CEDRIC

int main()
{
	std::cout << "Initializing MCMCAlgorithm object---------------" << std::endl;
	int samples = 500;
	int thining = 10;
	int useSamples = 100;
	std::cout << "\t# Samples: " << samples << "\n";
	std::cout << "\tThining: " << thining << "\n";
	std::cout << "\t# Samples used: " << useSamples << "\n";
	MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thining, 100, true, true, true);
	//mcmc.setRestartFileSettings("RestartFile.txt", 20, true);
	std::cout << "Done!-------------------------------\n\n\n";


	std::cout << "initialize Genome object--------------------------" << std::endl;
	bool withPhi = true;

	Genome genome;
	//genome.readFasta("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/data/twoMixtures/simulatedAllUniqueR_unevenMixtures.fasta");
	genome.readFasta("F:/GitHub/RibModelDev/data/twoMixtures/simulatedAllUniqueR.fasta");
	//genome.readFasta("E:/RibosomeModel/RibModelDev/data/twoMixtures/simulatedAllUniqueR_unevenMixtures.fasta");
	if(withPhi)
	{
		genome.readObservedPhiValues("F:/GitHub/RibModelDev/data/twoMixtures/simulatedAllUniqueR_phi_withPhiSet.csv", false);
		//genome.readObservedPhiValues("E:/RibosomeModel/RibModelDev/data/twoMixtures/simulatedAllUniqueR_phi_unevenMixtures.csv", false);
	}

	std::cout << "Done!-------------------------------\n\n\n";
	std::cout << "Initializing shared parameter variables---------------\n";

	std::cout << "Done!-------------------------------\n\n\n";
	std::cout << "Initializing shared parameter variables---------------\n";
	std::vector<unsigned> geneAssignment(genome.getGenomeSize());

	unsigned numMixtures = 2;
	std::vector<double> sphi_init(numMixtures, 1);

	/* For 2 mixture */
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		geneAssignment[i] = ( ((double)rand() / (double)RAND_MAX) < 0.5 ? 0u : 1u );

	}
	std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
	std::cout << "Done!------------------------\n\n\n";

	std::cout << "initialize ROCParameter object" << std::endl;
	std::string mixDef = ROCParameter::allUnique;
	ROCParameter parameter(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

	for (unsigned i = 0u; i < numMixtures; i++)
	{
		unsigned selectionCategry = parameter.getSelectionCategory(i);
		std::cout << "Sphi_init for selection category " << selectionCategry << ": " << sphi_init[selectionCategry] << std::endl;
	}
	std::cout << "\t# mixtures: " << numMixtures << "\n";
	std::cout << "\tmixture definition: " << mixDef << "\n";

	std::vector<std::string> files(2);
	files[0] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_mutation0.csv");
	files[1] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_mutation1.csv");
	//files[0] = std::string("E:/RibosomeModel/RibModelDev/data/twoMixtures/simulated_mutation0.csv");
	//files[1] = std::string("E:/RibosomeModel/RibModelDev/data/twoMixtures/simulated_mutation1.csv");
	parameter.initMutationCategories(files, parameter.getNumMutationCategories());
	files[0] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_selection0.csv");
	files[1] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_selection1.csv");
	//files[0] = std::string("E:/RibosomeModel/RibModelDev/data/twoMixtures/simulated_selection0.csv");
	//files[1] = std::string("E:/RibosomeModel/RibModelDev/data/twoMixtures/simulated_selection1.csv");
	parameter.initSelectionCategories(files, parameter.getNumSelectionCategories());

	parameter.InitializeSynthesisRate(genome, sphi_init[0]);
	//std::vector<double> phiVals = parameter.readPhiValues("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrCleft_phi_est.csv");
	//parameter.InitializeSynthesisRate(phiVals);
	std::cout << "done initialize ROCParameter object" << std::endl;


	std::cout << "Initializing ROCModel object\n";

	ROCModel model(withPhi);
	model.setParameter(parameter);


	std::cout << "starting MCMC for ROC" << std::endl;
	mcmc.run(genome, model, 1, 0);
	std::cout << std::endl << "Finished MCMC for ROC" << std::endl;
	double temp[] = { 0.025, 0.5, 0.975 };
	std::vector<double> probs(temp, temp + sizeof(temp) / sizeof(double));
	std::vector<double> quants = parameter.getCodonSpecificQuantile(1, 100, std::string("GCA"), 0, probs, true);

	std::cout << std::endl << "Exiting" << std::endl;
}
#endif // CEDRIC

#ifdef GABE
int main()
{


	if (1)
	{
		//testSequenceSummary();
		//testGene();
		testGenome("/Users/roxasoath1/Desktop/RibModelDevScripts/RibModelDev/data/UnitTestingData");

		/*
		Genome genome;
		genome.readRFPFile("/Users/roxasoath1/Desktop/RibModelDevScripts/RibModelDev/data/rfp/rfp.counts.by.codon.and.gene.GSE63789.wt.csv");
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
		{
			geneAssignment[i] = 0u;
		}
		unsigned numMixtures = 1;
		std::vector<double> sphi_init(numMixtures, 2);
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		RFPParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, "allUnique");

		std::vector<std::string> files;
		files.push_back("/Users/roxasoath1/Desktop/TONEWTON/RFPAlphaValues.csv");
		tmp.initMutationSelectionCategories(files, 1, RFPParameter::alp);
		files[0] = "/Users/roxasoath1/Desktop/TONEWTON/RFPLambdaPrimeValues.csv";
		tmp.initMutationSelectionCategories(files, 1, RFPParameter::lmPri);
		std::vector<double> phi = tmp.readPhiValues("/Users/roxasoath1/Desktop/TONEWTON/RFPPsiValues.csv");
		tmp.InitializeSynthesisRate(phi);


		RFPModel model;

		model.setParameter(tmp);

		std::cout <<"init done\n";
		model.simulateGenome(genome);
		std::cout <<"writing file\n";
		genome.writeRFPFile("/Users/roxasoath1/Desktop/RibModelDevScripts/RibModelDev/data/rfp/simulatedRFPData.csv", true);
*/
		exit(1);
	}
	std::string modelToRun = "RFP"; //can also be ROC or FONSE
	bool withPhi = false;
	bool fromRestart = true;
	unsigned numMixtures = 1;


	std::cout << "Initializing MCMCAlgorithm object---------------" << std::endl;
	int samples = 10;
	int thining = 10;
	int useSamples = 100;
	std::cout << "\t# Samples: " << samples << "\n";
	std::cout << "\tThining: " << thining << "\n";
	std::cout << "\t # Samples used: " << useSamples << "\n";
	MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thining, 10, true, true, true);
	mcmc.setRestartFileSettings("RestartFile.txt", 20, true);
	std::cout << "Done!-------------------------------\n\n\n";




	if (modelToRun == "ROC")
	{
		std::cout << "Initializing Genome object--------------------------" << std::endl;
		Genome genome;
		genome.readFasta("/Users/roxasoath1/Desktop/RibModelDevScripts/RibModelDev/data/twoMixtures/simulatedAllUniqueR.fasta");
		if (withPhi)
		{
			genome.readObservedPhiValues("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/simulatedAllUniqueR_phi.csv", false);
		}
		std::cout << "Done!-------------------------------\n\n\n";



		std::cout << "Initializing shared parameter variables---------------\n";
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		std::cout << "Done!------------------------\n\n\n";



		std::cout << "Initializing ROCParameter object--------------------\n" << std::endl;
		ROCParameter parameter;

		if (fromRestart)
		{
			ROCParameter tmp("/Users/roxasoath1/Desktop/RibModelFramework/DevRscripts/10restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = ROCParameter::allUnique;
			ROCParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++)
			{
				unsigned selectionCategry = tmp.getSelectionCategory(i);
				std::cout << "Sphi_init for selection category " << selectionCategry << ": " << sphi_init[selectionCategry] << std::endl;
			}
			std::cout << "\t# mixtures: " << numMixtures << "\n";
			std::cout << "\tmixture definition: " << mixDef << "\n";

			std::vector<std::string> files(2);
			files[0] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_mutation0.csv");
			files[1] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_mutation1.csv");
			tmp.initMutationCategories(files, tmp.getNumMutationCategories());
			files[0] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_selection0.csv");
			files[1] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_selection1.csv");
			tmp.initSelectionCategories(files, tmp.getNumSelectionCategories());

			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			//std::vector<double> phiVals = parameter.readPhiValues("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrCleft_phi_est.csv");
			//parameter.InitializeSynthesisRate(phiVals);
		}
		std::cout << "Done!--------------------------------\n\n\n" << std::endl;



		std::cout << "Initializing ROCModel object--------------------------\n";

		ROCModel model;
		model.setParameter(parameter);
		std::cout << "Done!----------------------------------\n\n\n" << std::endl;


		std::cout << "Running MCMC.............\n" << std::endl;
		mcmc.run(genome, model, 1, 0);
		std::cout << "Done!----------------------------------\n\n\n" << std::endl;
	} //END OF ROC
	else if (modelToRun == "RFP")
	{
		std::cout << "Initializing Genome object--------------------------" << std::endl;
		Genome genome;
		genome.readRFPFile("/Users/roxasoath1/Desktop/RibModelDevScripts/RibModelDev/data/rfp/rfp.counts.by.codon.and.gene.GSE63789.wt.csv");
		std::cout << "Done!-------------------------------\n\n\n";



		std::cout << "Initializing shared parameter variables---------------\n";
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		std::cout << "Done!------------------------\n\n\n";



		std::cout << "Initializing RFPParameter object--------------------\n" << std::endl;
		RFPParameter parameter;

		if (fromRestart)
		{
			RFPParameter tmp("/Users/roxasoath1/Desktop/RibModelFramework/10_restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = Parameter::allUnique;
			RFPParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++) {
				unsigned selectionCategry = tmp.getSelectionCategory(i);
				std::cout << "Sphi_init for selection category " << selectionCategry << ": " <<
				sphi_init[selectionCategry] << std::endl;
			}
			std::cout << "\t# mixtures: " << numMixtures << "\n";
			std::cout << "\tmixture definition: " << mixDef << "\n";

			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			parameter = tmp;
		}
		std::cout << "Done!--------------------------------\n\n\n" << std::endl;



		std::cout << "Initializing RFPModel object--------------------------\n";

		RFPModel model;
		model.setParameter(parameter);
		std::cout << "Done!----------------------------------\n\n\n" << std::endl;


		std::cout << "Running MCMC.............\n" << std::endl;
		mcmc.run(genome, model, 1, 0);
		std::cout << "Done!----------------------------------\n\n\n" << std::endl;

	} //END OF RFP
	else if (modelToRun == "FONSE")
	{
		std::cout << "initialize Genome object--------------------------" << std::endl;
		Genome genome;
		genome.readFasta("/Users/roxasoath1/Desktop/RibModelDevScripts/RibModelDev/data/FONSE/genome_2000.fasta");
		std::cout << "Done!-------------------------------\n\n\n";



		std::cout << "Initializing shared parameter variables---------------\n";
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		std::cout << "Done!------------------------\n\n\n";



		FONSEParameter parameter;
		std::cout << "initialize Parameter object" << std::endl;
		if (fromRestart)
		{
			FONSEParameter tmp("/Users/roxasoath1/Desktop/RibModelDevScripts/RibModelDev/DevRscripts/10restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = ROCParameter::allUnique;
			FONSEParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++) {
				unsigned selectionCategry = tmp.getSelectionCategory(i);
				std::cout << "Sphi_init for selection category " << selectionCategry << ": " <<
				sphi_init[selectionCategry] << std::endl;
			}
			std::cout << "\t# mixtures: " << numMixtures << "\n";
			std::cout << "\tmixture definition: " << mixDef << "\n";

			std::vector<std::string> files(1);
			files[0] = std::string(
					"/Users/roxasoath1/Desktop/RibModelDevScripts/RibModelDev/data/FONSE/genome_2000.mutation.csv");
			tmp.initMutationCategories(files, tmp.getNumMutationCategories());
			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			//std::vector<double> phiVals = parameter.readPhiValues("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrCleft_phi_est.csv");
			//parameter.InitializeSynthesisRate(phiVals);
			parameter = tmp;
			std::cout << "done initialize Parameter object" << std::endl;
		}


		std::cout << "Initializing Model object\n";

		FONSEModel model;
		model.setParameter(parameter);


		std::cout << "starting MCMC for ROC" << std::endl;
		mcmc.run(genome, model, 4, 0);
		std::cout << std::endl << "Finished MCMC for ROC" << std::endl;

	}
}

#endif // GABE



#ifdef JEREMY
int main()
{
	unsigned index;
	bool fromRestart = false;
	std::string modelToRun = "ROC";
	bool withPhi = false;

	
	std::cout << "Initializing MCMCAlgorithm object---------------" << std::endl;
	int samples = 100;
	int thining = 10;
	int useSamples = 100;
	unsigned numMixtures = 1;
	std::cout << "\t# Samples: " << samples << "\n";
	std::cout << "\tThining: " << thining << "\n";
	std::cout << "\t # Samples used: " << useSamples << "\n";
	MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thining, 10, true, true, true);
	mcmc.setRestartFileSettings("RestartFile.txt", 20, true);
	std::cout << "Done!-------------------------------\n\n\n";

	if (modelToRun == "FONSE") {
		std::cout << "initialize Genome object--------------------------" << std::endl;
		Genome genome;
		genome.readFasta("C:/Users/Jeremy/Documents/GitHub/RibModelDev/data/FONSE/nse2000.fasta");
		std::cout << "Done!-------------------------------\n\n\n";


		std::cout << "Initializing shared parameter variables---------------\n";
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());

		std::vector<double> sphi_init(numMixtures, 1);

		/* For 1 mixture */
		for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
		{
			geneAssignment[i] = 0u;
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		std::cout << "Done!------------------------\n\n\n";

		//ROCParameter parameter;
		FONSEParameter parameter;
		std::cout << "initialize Parameter object" << std::endl;
		std::string mixDef = Parameter::allUnique;
		if (fromRestart)
		{
			FONSEParameter tmp("C:/Users/Jeremy/Documents/GitHub/RibModelDev/DevRScripts/10restartfile.rst");
			parameter = tmp;
		}
		else
		{
			//ROCParameter tmp(sphi_init, nu/mMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);
			FONSEParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++)
			{
				unsigned selectionCategry = tmp.getSelectionCategory(i);
				std::cout << "Sphi_init for selection category " << selectionCategry << ": " << sphi_init[selectionCategry] << std::endl;
			}
			std::cout << "\t# mixtures: " << numMixtures << "\n";
			std::cout << "\tmixture definition: " << mixDef << "\n";

			//std::vector<std::string> files(1);
			//files[0] = std::string("C:/Users/Jeremy/Documents/GitHub/RibModelDev/data/FONSE/Scereviciae.mut.csv");
			//tmp.initMutationCategories(files, tmp.getNumMutationCategories());
			//tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			std::vector<double> phiVals = parameter.readPhiValues("C:/Users/Jeremy/Documents/GitHub/RibModelDev/data/FONSE/nse2000.phi.csv");
			tmp.InitializeSynthesisRate(phiVals);
			parameter = tmp;
		}
		std::cout << "done initialize Parameter object" << std::endl;


		std::cout << "Initializing Model object\n";

		bool withPhi = true;
		FONSEModel model;
		//ROCModel model(withPhi);
		model.setParameter(parameter);


		std::cout << "starting MCMC for ROC" << std::endl;
		mcmc.run(genome, model, 4, 0);
		std::cout << std::endl << "Finished MCMC for ROC" << std::endl;
	}
	else if (modelToRun == "RFP")
	{
		std::cout << "Initializing Genome object--------------------------" << std::endl;
		Genome genome;
		genome.readRFPFile("C:/Users/Jeremy/Documents/GitHub/RibModelDev/data/rfp/rfp.counts.by.codon.and.gene.GSE63789.wt.csv");
		std::cout << "Done!-------------------------------\n\n\n";



		std::cout << "Initializing shared parameter variables---------------\n";
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		std::cout << "Done!------------------------\n\n\n";



		std::cout << "Initializing RFPParameter object--------------------\n" << std::endl;
		RFPParameter parameter;

		if (fromRestart)
		{
			RFPParameter tmp("/Users/roxasoath1/Desktop/RibModelFramework/DevRscripts/10restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = Parameter::allUnique;
			RFPParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++) {
				unsigned selectionCategry = tmp.getSelectionCategory(i);
				std::cout << "Sphi_init for selection category " << selectionCategry << ": " <<
					sphi_init[selectionCategry] << std::endl;
			}
			std::cout << "\t# mixtures: " << numMixtures << "\n";
			std::cout << "\tmixture definition: " << mixDef << "\n";

			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			parameter = tmp;
		}
		std::cout << "Done!--------------------------------\n\n\n" << std::endl;



		std::cout << "Initializing RFPModel object--------------------------\n";

		RFPModel model;
		model.setParameter(parameter);
		std::cout << "Done!----------------------------------\n\n\n" << std::endl;


		std::cout << "Running MCMC.............\n" << std::endl;
		mcmc.run(genome, model, 1, 0);
		std::cout << "Done!----------------------------------\n\n\n" << std::endl;
	} //END OF RFP
	if (modelToRun == "ROC")
	{
		std::cout << "Initializing Genome object--------------------------" << std::endl;
		Genome genome;
		genome.readFasta("C:/Users/Jeremy/Documents/GitHub/RibModelDev/data/realGenomes/Skluyveri.fasta");
		if (withPhi)
		{
			genome.readObservedPhiValues("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/simulatedAllUniqueR_phi.csv", false);
		}
		std::cout << "Done!-------------------------------\n\n\n";



		std::cout << "Initializing shared parameter variables---------------\n";
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		std::cout << "Done!------------------------\n\n\n";



		std::cout << "Initializing ROCParameter object--------------------\n" << std::endl;
		ROCParameter parameter;

		if (fromRestart)
		{
			ROCParameter tmp("/Users/roxasoath1/Desktop/RibModelFramework/DevRscripts/10restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = ROCParameter::allUnique;
			ROCParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++)
			{
				unsigned selectionCategry = tmp.getSelectionCategory(i);
				std::cout << "Sphi_init for selection category " << selectionCategry << ": " << sphi_init[selectionCategry] << std::endl;
			}
			std::cout << "\t# mixtures: " << numMixtures << "\n";
			std::cout << "\tmixture definition: " << mixDef << "\n";

			/*std::vector<std::string> files(2);
			files[0] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_mutation0.csv");
			files[1] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_mutation1.csv");
			tmp.initMutationCategories(files, tmp.getNumMutationCategories());
			files[0] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_selection0.csv");
			files[1] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_selection1.csv");
			tmp.initSelectionCategories(files, tmp.getNumSelectionCategories());

			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			*/
			//std::vector<double> phiVals = parameter.readPhiValues("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrCleft_phi_est.csv");
			//parameter.InitializeSynthesisRate(phiVals);
			parameter = tmp;
		}
		std::cout << "Done!--------------------------------\n\n\n" << std::endl;



		std::cout << "Initializing ROCModel object--------------------------\n";

		ROCModel model;
		model.setParameter(parameter);
		std::cout << "Done!----------------------------------\n\n\n" << std::endl;


		std::cout << "Running MCMC.............\n" << std::endl;
		std::cout << numMixtures << std::endl;
		mcmc.run(genome, model, 1, 0);
		std::cout << "Done!----------------------------------\n\n\n" << std::endl;
	} //END OF ROC
}

#endif // JEREMY


#ifdef HOLLIS
int main()
{
	std::string pathBegin = "/Users/hollisbui/";

	// UNIT TESTING
	//testUtility();
	//testSequenceSummary();
	//testGene();
	//testGenome(pathBegin + "RibModelFramework/tests/testthat/UnitTestingData");
	//testCovarianceMatrix();
	//testParameter();
	//testTrace();
	//testRFPParameter();
	testMCMCAlgorithm();
	exit(0);


	std::string modelToRun = "RFP"; //can be RFP, ROC or FONSE
	bool withPhi = false;
	bool fromRestart = false;
	unsigned numMixtures = 1;


	my_print("Initializing MCMCAlgorithm object---------------\n");
	unsigned samples = 10;
	unsigned thining = 10;
	int useSamples = 100;
	my_print("\t# Samples: %\n", samples);
	my_print("\tThining: %\n", thining);
	my_print("\t # Samples used: %\n", useSamples);
	MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thining, 10, true, true, true);
	//mcmc.setRestartFileSettings("RestartFile.txt", 20, true);
	my_print("Done!-------------------------------\n\n\n");


	if (modelToRun == "ROC")
	{
		my_print("Initializing Genome object--------------------------\n");
		Genome genome;
		genome.readFasta(pathBegin + "RibModelDev/data/twoMixtures/simulatedAllUniqueR.fasta");
		if (withPhi)
		{
			genome.readObservedPhiValues(pathBegin + "RibModelFramework/ribModel/data/simulatedAllUniqueR_phi.csv", false);
		}
		my_print("Done!-------------------------------\n\n\n");


		my_print("Initializing shared parameter variables---------------\n");
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		my_print("Done!------------------------\n\n\n");


		my_print("Initializing ROCParameter object--------------------\n\n");
		ROCParameter parameter;

		if (fromRestart)
		{
			ROCParameter tmp(pathBegin + "RibModelFramework/DevRscripts/10restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = ROCParameter::allUnique;
			ROCParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++)
			{
				unsigned selectionCategry = tmp.getSelectionCategory(i);
				my_print("Sphi_init for selection category %: %\n", selectionCategry, sphi_init[selectionCategry]);
			}
			my_print("\t# mixtures: %\n", numMixtures);
			my_print("\tmixture definition: %\n", mixDef);

			std::vector<std::string> files(2);
			files[0] = std::string(pathBegin + "RibModelDev/data/twoMixtures/simulated_mutation0.csv");
			files[1] = std::string(pathBegin + "RibModelDev/data/twoMixtures/simulated_mutation1.csv");
			tmp.initMutationCategories(files, tmp.getNumMutationCategories());
			files[0] = std::string(pathBegin + "RibModelDev/data/twoMixtures/simulated_selection0.csv");
			files[1] = std::string(pathBegin + "RibModelDev/data/twoMixtures/simulated_selection1.csv");
			tmp.initSelectionCategories(files, tmp.getNumSelectionCategories());

			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			//std::vector<double> phiVals =
			// parameter.readPhiValues(pathBegin + "RibModelDev/data/realGenomes/Skluyveri_ChrA_ChrCleft_phi_est.csv");
			//parameter.InitializeSynthesisRate(phiVals);
			parameter = tmp;
		}
		my_print("Done!--------------------------------\n\n\n");


		my_print("Initializing ROCModel object--------------------------\n");
		ROCModel model;
		model.setParameter(parameter);
		my_print("Done!----------------------------------\n\n\n");


		my_print("Running MCMC.............\n\n");
		mcmc.run(genome, model, 1, 0);
		my_print("Done!----------------------------------\n\n\n");
	} //END OF ROC
	else if (modelToRun == "RFP")
	{
		my_print("Initializing Genome object--------------------------\n");
		Genome genome;
		genome.readRFPFile(pathBegin + "RibModelDev/data/rfp/rfp.counts.by.codon.and.gene.GSE63789.wt.csv");
		my_print("Done!-------------------------------\n\n\n");


		my_print("Initializing shared parameter variables---------------\n");
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		my_print("Done!------------------------\n\n\n");


		my_print("Initializing RFPParameter object--------------------\n\n");
		RFPParameter parameter;

		if (fromRestart)
		{
			RFPParameter tmp(pathBegin + "RibModelFramework/10_restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = Parameter::allUnique;
			RFPParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++) {
				unsigned selectionCategry = tmp.getSelectionCategory(i);
				my_print("Sphi_init for selection category %: %\n", selectionCategry, sphi_init[selectionCategry]);
			}
			my_print("\t# mixtures: %\n", numMixtures);
			my_print("\tmixture definition: %\n", mixDef);

			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			parameter = tmp;
		}
		my_print("Done!--------------------------------\n\n\n");


		my_print("Initializing RFPModel object--------------------------\n");
		RFPModel model;
		model.setParameter(parameter);
		my_print("Done!----------------------------------\n\n\n");


		my_print("Running MCMC.............\n\n");
		mcmc.run(genome, model, 1, 0);
		my_print("Done!----------------------------------\n\n\n");

	} //END OF RFP
	else if (modelToRun == "FONSE")
	{
		my_print("initialize Genome object--------------------------\n");
		Genome genome;
		genome.readFasta(pathBegin + "RibModelDev/data/singleMixture/genome_2000.fasta");
		my_print("Done!-------------------------------\n\n\n");


		my_print("Initializing shared parameter variables---------------\n");
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		my_print("Done!------------------------\n\n\n");


		FONSEParameter parameter;
		my_print("initialize Parameter object\n");
		if (fromRestart)
		{
			FONSEParameter tmp(pathBegin + "RibModelDev/DevRscripts/10restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = ROCParameter::allUnique;
			FONSEParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++) {
				unsigned selectionCategory = tmp.getSelectionCategory(i);
				my_print("Sphi_init for selection category %: %\n", selectionCategory, sphi_init[selectionCategory]);
			}
			my_print("\t# mixtures: %\n", numMixtures);
			my_print("\tmixture definition: %\n", mixDef);

			std::vector<std::string> files(1);
			files[0] = std::string(pathBegin + "RibModelDev/data/singleMixture/genome_2000.mutation.csv");
			tmp.initMutationCategories(files, tmp.getNumMutationCategories());
			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			//std::vector<double> phiVals = parameter.readPhiValues(
			// pathBegin + "RibModelDev/data/realGenomes/Skluyveri_ChrA_ChrCleft_phi_est.csv");
			//parameter.InitializeSynthesisRate(phiVals);
			parameter = tmp;
			my_print("done initialize Parameter object\n");
		}


		my_print("Initializing Model object\n");
		FONSEModel model;
		model.setParameter(parameter);
		my_print("Done!------------------------\n\n\n");


		my_print("starting MCMC for ROC\n");
		mcmc.run(genome, model, 4, 0);
		my_print("\nFinished MCMC for ROC\n");

	}
}

#endif // Hollis
