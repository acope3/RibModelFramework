#ifndef GENOME_H
#define GENOME_H

#include <vector>
#include <string>

//#include <Rcpp.h>

#include "../include/Gene.h"
class Genome
{
    private:
        std::vector<Gene> genes;

    public:
        //constructor/destructor
        explicit Genome();
        virtual ~Genome();
        Genome(const Genome& other);
        Genome& operator=(const Genome& other);

        void readFasta(char* filename);
        void writeFasta(char* filename);
        void addGene(const Gene& gene);

        Gene& getGene(int index);
        Gene& getGene(std::string id);

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
