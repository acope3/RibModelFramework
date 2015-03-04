#include <iostream>
#include <time.h>
#include <fstream>

#include "../include/Genome.h"
#include "../include/Gene.h"
#include "../include/SequenceSummary.h"
#include "../include/ROCParameter.h"
#include "../include/ROCModel.h"
#include "../include/MCMCAlgorithm.h"

using namespace std;

int main()
{
    std::cout << "Hello world!" << std::endl << std::endl;


    std::cout << "reading fasta file" << std::endl;
    Genome genome;
    //genome.readFasta("../../inst/testGenome.fasta");
    genome.readFasta("/home/clandere/CodonUsageBias/organisms/yeast/data/LKluyveri/Skluyveri_chromosomeA.fasta");
    //genome.readFasta("/home/clandere/CodonUsageBias/organisms/yeast/data/LKluyveri/Skluyveri.fasta");
    //genome.writeFasta("../../inst/resGenome.fasta

/*
    Gene gene = genome.getGene(0);
    for(int i = 0; i < genome.getGenomeSize(); i++)
    {
        //gene.setDeltaEtaCategory(1);
        std::cout << gene.getDeltaEtaCategory() << std::endl;
        //gene.setMutationCategory(1);
        std::cout << gene.getMutationCategory() << std::endl;
    }

    std::cout << ROCParameter::dM << std::endl;
    std::cout << ROCParameter::dEta << std::endl;

    std::cout << gene.toString() << std::endl;
    //SequenceSummary seqsum = gene.getSequenceSummary();

    unsigned* test = SequenceSummary::AAToCodonRange('A', true);
    std::cout << test[0] << " to " << test[1] << std::endl;
    test = SequenceSummary::AAToCodonRange('C', true);
    std::cout << test[0] << " to " << test[1] << std::endl;
    test = SequenceSummary::AAToCodonRange('D', true);
    std::cout << test[0] << " to " << test[1] << std::endl;
*/

    int samples = 2000;
    int thining = 10;
    int useSamples = 1000;

    ROCModel model = ROCModel();
    ROCParameter parameter = ROCParameter(1, 1, genome.getGenomeSize(), 2.0);
    parameter.InitializeExpression(genome, 2.0);
    MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thining, true, false, true);

    std::cout << "starting MCMC" << std::endl;
    mcmc.run(genome, model, parameter);
    std::cout << std::endl << "Finish MCMC" << std::endl;


    std::cout << parameter.getSphiPosteriorMean(useSamples) << std::endl;
    std::ofstream phiout("/home/clandere/CodonUsageBias/organisms/yeast/results/test.phi");
    for(int n = 0; n < genome.getGenomeSize(); n++)
    {
        phiout << genome.getGene(n).getId() << "," << parameter.getExpressionPosteriorMean(useSamples, n) << std::endl;
    }

    return 0;
}
