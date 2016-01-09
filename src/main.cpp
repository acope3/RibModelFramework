#include "include/MCMCAlgorithm.h"
#include <iostream>
#include <vector>

#ifdef CEDRIC

int main()
{
	std::cout << "Initializing MCMCAlgorithm object---------------" << std::endl;
	int samples = 1000;
	int thining = 10;
	int useSamples = 100;
	std::cout << "\t# Samples: " << samples << "\n";
	std::cout << "\tThining: " << thining << "\n";
	std::cout << "\t# Samples used: " << useSamples << "\n";
	MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thining, 10, true, true, true);
	//mcmc.setRestartFileSettings("RestartFile.txt", 20, true);
	std::cout << "Done!-------------------------------\n\n\n";


	std::cout << "initialize Genome object--------------------------" << std::endl;
	bool withPhi = false;

	Genome genome;
	//genome.readFasta("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/data/twoMixtures/simulatedAllUniqueR_unevenMixtures.fasta");
	genome.readFasta("F:/GitHub/RibModelDev/data/twoMixtures/simulatedAllUniqueR_unevenMixtures.fasta");
	if(withPhi)
	{
		genome.readObservedPhiValues("F:/GitHub/RibModelDev/data/twoMixtures/simulatedAllUniqueR_phi_unevenMixtures.csv", false);
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
	parameter.initMutationCategories(files, parameter.getNumMutationCategories());
	files[0] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_selection0.csv");
	files[1] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_selection1.csv");
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

	std::cout << std::endl << "Exiting" << std::endl;
}
#endif // CEDRIC

#ifdef GABE
int main()
{
	std::string modelToRun = "RFP"; //can also be ROC or FONSE
	bool withPhi = false;
	bool fromRestart = false;
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

		Trace trace = parameter.getTraceObject();
		std::string codon = "ATG";
		std::vector<double> v;
		v = trace.getCodonSpecificParameterTraceByMixtureElementForCodon(0, codon, 0);
		for (int i = 0; i < v.size(); i++)
		{
			std::cout <<v[i] <<"\n";
		}
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
	std::string modelToRun = "FONSE";


	
	std::cout << "Initializing MCMCAlgorithm object---------------" << std::endl;
	int samples = 10;
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
		genome.readFasta("C:/Users/Jeremy/Documents/GitHub/RibModelDev/data/FONSE/genome_2000.fasta");
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
		std::string mixDef = ROCParameter::selectionShared;
		if (fromRestart)
		{
			FONSEParameter tmp("C:/Users/Jeremy/Documents/GitHub/RibModelDev/DevRScripts/10restartfile.rst");
			parameter = tmp;
		}
		else
		{
			//ROCParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);
			FONSEParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++)
			{
				unsigned selectionCategry = tmp.getSelectionCategory(i);
				std::cout << "Sphi_init for selection category " << selectionCategry << ": " << sphi_init[selectionCategry] << std::endl;
			}
			std::cout << "\t# mixtures: " << numMixtures << "\n";
			std::cout << "\tmixture definition: " << mixDef << "\n";

			std::vector<std::string> files(1);
			files[0] = std::string("C:/Users/Jeremy/Documents/GitHub/RibModelDev/data/FONSE/genome_2000.mutation.csv");
			tmp.initMutationCategories(files, tmp.getNumMutationCategories());
			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			//std::vector<double> phiVals = parameter.readPhiValues("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrCleft_phi_est.csv");
			//parameter.InitializeSynthesisRate(phiVals);
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
}

#endif // JEREMY


