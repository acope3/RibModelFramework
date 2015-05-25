#ifndef GENOME_H
#define GENOME_H

#include <vector>
#include <string>

//#include <Rcpp.h>

#include "../include/Gene.h"


//IMPORTANT NOTE: forward declarations used. Includes are in genome.cpp. 
//Used to solve circular dependices. See http://www.cplusplus.com/forum/general/125/
//for more information.
class ROCParameter;
class ROCModel;


class Genome
{
    private:
        std::vector<Gene> genes;
				std::vector<Gene> simulatedGenes;
    public:
        //constructor/destructor
        explicit Genome();
        virtual ~Genome();
        Genome(const Genome& other);
        Genome& operator=(const Genome& other);

        void readFasta(char* filename);
        void writeFasta(char* filename, bool simulated = false);
        void addGene(const Gene& gene);
        void getCountsForAA(char aa, unsigned codonCounts[][5]);
				std::vector <Gene> getGenes() {return genes;}
				std::vector <Gene> getSimulatedGenes() {return simulatedGenes;}
        Gene& getGene(int index);
        Gene& getGene(std::string id);
				void simulateGenome(ROCParameter& parameter, ROCModel& model);
        int getGenomeSize() {return genes.size();}

    protected:
};

//using namespace Rcpp;
//RCPP_MODULE(mGenome)
//{
//    class_<Genome>("Genome")
//    .constructor<>()
//    .constructor<std::string, std::string, std::string>()
//
//    .method("readFasta", &Genome::readFasta)
//    .method("writeFasta", &Genome::writeFasta)
//    ;
//}

#endif // GENOME_H
