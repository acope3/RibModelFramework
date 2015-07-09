#ifndef GENOME_H
#define GENOME_H

#include <vector>
#include <string>

//#include <Rcpp.h>

#include "Gene.h"


//IMPORTANT NOTE: forward declarations used. Includes are in genome.cpp. 
//Used to solve circular dependices. See http://www.cplusplus.com/forum/general/125/
//for more information.
class Model;


class Genome
{
	private:
		std::vector<Gene> genes;
		std::vector<Gene> simulatedGenes;
	public:
		//constructor/destructor
		explicit Genome();
		virtual ~Genome();
		Genome& operator=(const Genome& other);

		void readFasta(std::string filename, bool Append = false);
		void writeFasta(std::string filename, bool simulated = false);
		void readRFPFile(std::string filename);
		void addGene(const Gene& gene);
		std::vector<unsigned> getCodonCountsPerGene(std::string codon);
		std::vector <Gene> getGenes() {return genes;}
		std::vector <Gene> getSimulatedGenes() {return simulatedGenes;}
		Gene& getGene(unsigned index);
		Gene& getGene(std::string id);
		void simulateGenome(Model& model);
		unsigned getGenomeSize() {return genes.size();}
		void clear();
		Genome getGenomeForGeneIndicies(std::vector <unsigned> indicies);
		Genome getGenomeForGeneIndiciesR(std::vector <unsigned> indicies);
		bool checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound);
		//R wrapper functions
		Gene& getGeneByIndex(unsigned index);

	protected:
};

#endif // GENOME_H
