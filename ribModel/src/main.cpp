#include "include/MCMCAlgorithm.h"
#include <iostream>
#include <vector>

#ifdef CEDRIC
int main()
{
	unsigned index;
	std::cout << "Initializing MCMCAlgorithm object---------------" << std::endl;
	int samples = 1000;
	int thining = 10;
	int useSamples = 100;
	std::cout << "\t# Samples: " << samples << "\n";
	std::cout << "\tThining: " << thining << "\n";
	std::cout << "\t # Samples used: " << useSamples << "\n";
	MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thining, 10, true, true, true);
	//mcmc.setRestartFileSettings("RestartFile.txt", 20, true);
	std::cout << "Done!-------------------------------\n\n\n";


	std::cout << "initialize Genome object--------------------------" << std::endl;
	Genome genome;
	genome.readFasta("F:/GitHub/RibModelFramework/data/twoMixtures/simulatedAllUniqueR.fasta");
	genome.readObservedPhiValues("F:/GitHub/RibModelFramework/data/twoMixtures/simulatedAllUniqueR_phi.csv", false);
	std::cout << "Done!-------------------------------\n\n\n";
	std::cout << "Initializing shared parameter variables---------------\n";

	std::cout << "Done!-------------------------------\n\n\n";
	std::cout << "Initializing shared parameter variables---------------\n";
	std::vector<unsigned> geneAssignment(genome.getGenomeSize());

	unsigned numMixtures = 1;
	std::vector<double> sphi_init(numMixtures, 1);

	/* For 1 mixture */
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		geneAssignment[i] = 0u;
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
	files[0] = std::string("F:/GitHub/RibModelFramework/data/twoMixtures/simulated_mutation0.csv");
	files[1] = std::string("F:/GitHub/RibModelFramework/data/twoMixtures/simulated_mutation1.csv");
	parameter.initMutationCategories(files, parameter.getNumMutationCategories());
	files[0] = std::string("F:/GitHub/RibModelFramework/data/twoMixtures/simulated_selection0.csv");
	files[1] = std::string("F:/GitHub/RibModelFramework/data/twoMixtures/simulated_selection1.csv");
	parameter.initSelectionCategories(files, parameter.getNumSelectionCategories());

	parameter.InitializeSynthesisRate(genome, sphi_init[0]);
	//std::vector<double> phiVals = parameter.readPhiValues("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrCleft_phi_est.csv");
	//parameter.InitializeSynthesisRate(phiVals);
	std::cout << "done initialize ROCParameter object" << std::endl;


	std::cout << "Initializing ROCModel object\n";

	bool withPhi = true;
	ROCModel model(withPhi);
	model.setParameter(parameter);


	std::cout << "starting MCMC for ROC" << std::endl;
	mcmc.run(genome, model, 1, 0);
	std::cout << std::endl << "Finished MCMC for ROC" << std::endl;
}
#endif // CEDRIC

#ifdef GABE
int main()
{
	unsigned index;
	std::cout << "Initializing MCMCAlgorithm object---------------" << std::endl;
	int samples = 100;
	int thining = 10;
	int useSamples = 100;
	std::cout << "\t# Samples: " << samples << "\n";
	std::cout << "\tThining: " << thining << "\n";
	std::cout << "\t # Samples used: " << useSamples << "\n";
	MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thining, 10, true, true, true);
	mcmc.setRestartFileSettings("RestartFile.txt", 20, true);
	std::cout << "Done!-------------------------------\n\n\n";


	std::cout << "Initializing Genome object--------------------------" << std::endl;
	Genome genome;
	genome.readRFPFile("/Users/roxasoath1/Desktop/RibModelFramework/data/rfp/rfp.counts.by.codon.and.gene.GSE63789.wt.csv");
	std::cout << "Done!-------------------------------\n\n\n";


	std::cout << "Initializing shared parameter variables---------------\n";
	std::vector<unsigned> geneAssignment(genome.getGenomeSize());
	unsigned numMixtures = 1;
	std::vector<double> sphi_init(numMixtures, 1);

	/* For 1 mixture */
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		geneAssignment[i] = 0u;
	}
	std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
	std::cout << "Done!------------------------\n\n\n";

	
	RFPParameter parameter;
	std::cout << "Initializing RFPParameter object--------------------\n" << std::endl;
	std::string mixDef = Parameter::selectionShared;
	RFPParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

	for (unsigned i = 0u; i < numMixtures; i++)
	{
		unsigned selectionCategry = tmp.getSelectionCategory(i);
		std::cout << "Sphi_init for selection category " << selectionCategry << ": " << sphi_init[selectionCategry] << std::endl;
	}
	std::cout << "\t# mixtures: " << numMixtures << "\n";
	std::cout << "\tmixture definition: " << mixDef << "\n";

	tmp.InitializeSynthesisRate(genome, sphi_init[0]);
	std::cout << "Done!--------------------------------\n\n\n" << std::endl;


	std::cout << "Initializing RFPModel object--------------------------\n";

	bool withPhi = true;
	RFPModel model;
	model.setParameter(parameter);
	std::cout << "Done!----------------------------------\n\n\n" << std::endl;


	std::cout << "Running MCMC.............\n" << std::endl;
	mcmc.run(genome, model, 1, 0);
	std::cout << "Done!----------------------------------\n\n\n" << std::endl;
}

#endif // GABE


#ifdef JEREMY
int main()
{
	unsigned index;
	std::cout << "Initializing MCMCAlgorithm object---------------" << std::endl;
	int samples = 100;
	int thining = 10;
	int useSamples = 100;
	std::cout << "\t# Samples: " << samples << "\n";
	std::cout << "\tThining: " << thining << "\n";
	std::cout << "\t # Samples used: " << useSamples << "\n";
	MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thining, 10, true, true, true);
	mcmc.setRestartFileSettings("RestartFile.txt", 20, true);
	std::cout << "Done!-------------------------------\n\n\n";


	std::cout << "initialize Genome object--------------------------" << std::endl;
	Genome genome;
	genome.readFasta("C:/Users/Cedric/Documents/GitHub/RibModelFramework/data/realGenomes/Skluyveri.fasta");
	std::cout << "Done!-------------------------------\n\n\n";
	std::cout << "Initializing shared parameter variables---------------\n";

	std::cout << "Done!-------------------------------\n\n\n";
	std::cout << "Initializing shared parameter variables---------------\n";
	std::vector<unsigned> geneAssignment(genome.getGenomeSize());

	unsigned numMixtures = 1;
	std::vector<double> sphi_init(numMixtures, 1);

	/* For 1 mixture */
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		geneAssignment[i] = 0u;
	}
	std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
	std::cout << "Done!------------------------\n\n\n";

	ROCParameter parameter;
	std::cout << "initialize ROCParameter object" << std::endl;
	std::string mixDef = ROCParameter::selectionShared;
	ROCParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

	for (unsigned i = 0u; i < numMixtures; i++)
	{
		unsigned selectionCategry = tmp.getSelectionCategory(i);
		std::cout << "Sphi_init for selection category " << selectionCategry << ": " << sphi_init[selectionCategry] << std::endl;
	}
	std::cout << "\t# mixtures: " << numMixtures << "\n";
	std::cout << "\tmixture definition: " << mixDef << "\n";

	std::vector<std::string> files(2);
	files[0] = std::string("C:/Users/Cedric/Documents/GitHub/RibModelFramework/data/realGenomes/Skluyveri_mutation_ChrA.csv");
	files[1] = std::string("C:/Users/Cedric/Documents/GitHub/RibModelFramework/data/realGenomes/Skluyveri_mutation_ChrCleft.csv");
	tmp.initMutationCategories(files, tmp.getNumMutationCategories());
	tmp.InitializeSynthesisRate(genome, sphi_init[0]);
	//std::vector<double> phiVals = parameter.readPhiValues("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrCleft_phi_est.csv");
	//parameter.InitializeSynthesisRate(phiVals);
	parameter = tmp;
	std::cout << "done initialize ROCParameter object" << std::endl;


	std::cout << "Initializing ROCModel object\n";

	bool withPhi = true;
	ROCModel model(withPhi);
	model.setParameter(parameter);


	std::cout << "starting MCMC for ROC" << std::endl;
	mcmc.run(genome, model, 1, 0);
	std::cout << std::endl << "Finished MCMC for ROC" << std::endl;
}

#endif // JEREMY


