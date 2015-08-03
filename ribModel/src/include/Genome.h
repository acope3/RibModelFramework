#ifndef GENOME_H
#define GENOME_H

#include <vector>
#include <string>
#include <map>

#include "Gene.h"


class Model;


class Genome
{
	private:

		std::vector<Gene> genes;
		std::vector<Gene> simulatedGenes;

		unsigned numGenesWithPhi;
	public:

		//Constructors & destructors:
		explicit Genome();
		virtual ~Genome();
		Genome& operator=(const Genome& other);


		//File I/O functions:
		void readFasta(std::string filename, bool Append = false);
		void writeFasta(std::string filename, bool simulated = false);
		void readRFPFile(std::string filename);
		void writeRFPFile(std::string filename, bool simulated = false);
		void readObservedPhiValues(std::string filename, bool byId = true);


		//Gene functions:
		void addGene(const Gene& gene, bool simulated = false);
		std::vector <Gene> getGenes(bool simulated = false);
		Gene& getGene(unsigned index);
		Gene& getGene(std::string id);


		//Other functions:
		bool checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound);
		unsigned getGenomeSize();
		void clear();
		Genome getGenomeForGeneIndicies(std::vector <unsigned> indicies);
		std::vector<unsigned> getCodonCountsPerGene(std::string codon);



		//R wrapper functions:
		Gene& getGeneByIndex(unsigned index);
		Gene& getGeneById(std::string ID);
		Genome getGenomeForGeneIndiciesR(std::vector <unsigned> indicies);

	protected:
};

#endif // GENOME_H
