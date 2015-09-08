#ifndef GENOME_H
#define GENOME_H

#include <vector>
#include <string>
#include <map>

#include "Gene.h"
#include "CodonTable.h"


class Model;


class Genome
{
	private:

		std::vector<Gene> genes;
		std::vector<Gene> simulatedGenes;
		std::vector <unsigned> numGenesWithPhi;

	public:

		//Constructors & destructors:
		Genome();
		explicit Genome(unsigned tableId, bool splitAA = false);
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
		unsigned getNumGenesWithPhi(unsigned index);
		Gene& getGene(unsigned index, bool simulated = false);
		Gene& getGene(std::string id, bool simulated = false);


		//Other functions:
		bool checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound);
		unsigned getGenomeSize();
		void clear();
		Genome getGenomeForGeneIndicies(std::vector <unsigned> indicies, bool simulated = false); //NOTE: If simulated is true, it will return a genome with the simulated genes, but the returned genome's genes vector will contain the simulated genes.
		std::vector<unsigned> getCodonCountsPerGene(std::string codon);



		//R wrapper functions:
		Gene& getGeneByIndex(unsigned index, bool simulated = false);
		Gene& getGeneById(std::string ID, bool simulated = false);
		Genome getGenomeForGeneIndiciesR(std::vector <unsigned> indicies, bool simulated = false);

	protected:
};

#endif // GENOME_H
