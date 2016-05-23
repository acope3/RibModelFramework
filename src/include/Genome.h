#ifndef GENOME_H
#define GENOME_H

#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <cmath>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

#include "Gene.h"

class Model;
class Genome
{
	private:

		std::vector<Gene> genes;
		std::vector<Gene> simulatedGenes;
		std::vector <unsigned> numGenesWithPhi; //Number of phi sets is vector size, value is number of genes
												//with a phi value for that set. Values should currently be equal.

	public:

		//Constructors & Destructors:
		explicit Genome();
		Genome& operator=(const Genome& other);
		bool operator==(const Genome& other) const;
		virtual ~Genome();


		//File I/O Functions:
		void readFasta(std::string filename, bool Append = false);
		void writeFasta(std::string filename, bool simulated = false);
		void readRFPFile(std::string filename);
		void writeRFPFile(std::string filename, bool simulated = false);
		void readPANSEFile(std::string filename, bool Append = false);
		//TODO: Add writePANSEFile function
		void readObservedPhiValues(std::string filename, bool byId = true);


		//Gene Functions:
		void addGene(const Gene& gene, bool simulated = false);
		std::vector <Gene> getGenes(bool simulated = false);
		unsigned getNumGenesWithPhiForIndex(unsigned index);
		Gene& getGene(unsigned index, bool simulated = false);
		Gene& getGene(std::string id, bool simulated = false);


		//Other Functions:
		unsigned getGenomeSize(bool simulated = false);
		void clear();
		Genome getGenomeForGeneIndices(std::vector <unsigned> indices, bool simulated = false); //NOTE: If simulated is true, it will return a genome with the simulated genes, but the returned genome's genes vector will contain the simulated genes.
		std::vector<unsigned> getCodonCountsPerGene(std::string codon);


		//Testing Functions:
		std::vector <unsigned> getNumGenesWithPhi();
		void setNumGenesWithPhi(std::vector <unsigned> newVector);


		//R Section:

#ifndef STANDALONE

		bool checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound);
		Gene& getGeneByIndex(unsigned index, bool simulated = false);
		Gene& getGeneById(std::string ID, bool simulated = false);
		Genome getGenomeForGeneIndicesR(std::vector <unsigned> indices, bool simulated = false);

#endif //STANDALONE

	protected:
};

#endif // GENOME_H
