#include "include/Genome.h"

#include <iostream>     // std::cout, std::cerr
#include <cstring>
#include <fstream>
#include <ctime>
#include <sstream> // ostringstream

Genome::Genome()
{
	//ctor
}

Genome::~Genome()
{
	//dtor
}

Genome& Genome::operator=(const Genome& rhs)
{
	if (this == &rhs) return *this; // handle self assignment
	genes = rhs.genes;
	simulatedGenes = rhs.simulatedGenes;
	//assignment operator
	return *this;
}

void Genome::addGene(const Gene& gene, bool simulated)
{
	if (!simulated)
		genes.push_back(gene);
	else
		simulatedGenes.push_back(gene);
}


std::vector<unsigned> Genome::getCodonCountsPerGene(std::string codon)
{
	std::vector<unsigned> codonCounts;
	codonCounts.resize(genes.size());
	unsigned codonIndex = SequenceSummary::CodonToIndex(codon);
	for(unsigned i = 0u; i < genes.size(); i++)
	{
		Gene gene = genes[i];
		SequenceSummary seqsum = gene.getSequenceSummary();
		codonCounts[i] = seqsum.getCodonCountForCodon(codonIndex);
	}
	return codonCounts;
}

void Genome::writeFasta (std::string filename, bool simulated)
{
	try {
		std::ofstream Fout;
		Fout.open(filename.c_str());
		if (Fout.fail())
		{
			std::cerr <<"Error in Genome::writeFasta: Cannot open output Fasta file " << filename <<"\n";
		}
		else
		{
			if (simulated)
			{
				for(unsigned i = 0u; i < simulatedGenes.size(); i++)
				{
					Fout << ">" << simulatedGenes[i].getId() << " " << simulatedGenes[i].getDescription() <<"\n";
					for(int j = 0; j < simulatedGenes[i].length(); j++)
					{
						Fout << simulatedGenes[i].getNucleotideAt(j);
						if((j + 1) % 60 == 0) Fout << std::endl;
					}
					Fout << std::endl;
				}
			}
			else
			{
				for(unsigned i = 0u; i < genes.size(); i++)
				{
					Fout << ">" << genes[i].getId() << " " << genes[i].getDescription() << std::endl;
					for(int j = 0; j < genes[i].length(); j++)
					{
						Fout << genes[i].getNucleotideAt(j);
						if((j + 1) % 60 == 0) Fout << std::endl;
					}
					Fout << std::endl;
				}
			}
		} // end else
		Fout.close();
	} // end try
	catch(char* pMsg) { std::cerr << std::endl << "Exception:" << pMsg << std::endl; }
}

void Genome::readFasta(std::string filename, bool Append) // read Fasta format sequences
{
	try
	{
		if (!Append)
		{
			clear();
		}
		std::ifstream Fin;
		Fin.open(filename.c_str());
		if (Fin.fail())
		{
			std::cerr << "Error in Genome::readFasta: Cannot open input Fasta file " << filename << "\n";
		}
		else
		{

			bool fastaFormat = false;
			std::string buf;
			int newLine;


			Gene tmpGene;
			std::string tempSeq = "";
			for (;;)
			{
				// read a new line in every cycle
				std::getline(Fin, buf);
				/* The new line is marked with an integer, which corresponds
					 to one of the following cases:

					 1. New sequence line, which starts with "> ".
					 2. Sequence line, which starts with any charactor except for ">"
					 3. End of file, detected by function eof().
				 */

				if(buf[0] == '>' ) newLine=1;
				else if(Fin.eof()) newLine=3;
				else newLine=2;

				if( newLine == 1 )
				{ // this is a start of a new chain.
					if( !fastaFormat )
					{
						// if it is the first chain, just build a new chain object
						tmpGene.clear();
						fastaFormat = true;
					} else
					{
						// otherwise, need to store the old chain first, then build a new chain
						tmpGene.setSequence(tempSeq);
						//tmpGene.cleanSeq();
						addGene(tmpGene);

						tmpGene.clear();
						tempSeq = "";
					}
					tmpGene.setDescription( buf.substr(1,buf.size()-1) );
					int pos = buf.find(" ") - 1;
					tmpGene.setId( buf.substr(1,pos) );
				}

				if( newLine == 2 && fastaFormat )
				{ // sequence line
					tempSeq.append(buf);
					//tmpGene.seq.append(buf);
				}

				if( newLine == 3 )
				{ // end of file
					if( !fastaFormat ) throw std::string("Genome::readFasta throws: ") + std::string(filename) + std::string(" is not in Fasta format.");
					else
					{
						// otherwise, need to store the old chain first, then to
						// build a new chain
						tmpGene.setSequence(tempSeq);
						//tmpGene.cleanSeq();
						addGene(tmpGene);
						break;
					}
				}
			} // end while
		} // end else
	} // end try
	catch(char* pMsg) { std::cerr << std::endl << "Exception:" << pMsg << std::endl; }
}

void Genome::readRFPFile(std::string filename)
{
  std::ifstream Fin;
  Fin.open(filename.c_str());
  if (Fin.fail())
  {
    std::cerr << "Error in Genome::readRFPFile: Cannot open input RFP file " << filename << "\n";
  }

  std::string tmp;
  std::getline(Fin, tmp); //trash the first line
  std::string prevID = "";
  Gene tmpGene;
  SequenceSummary SS;
  bool first = true;
  std::string seq = "";

  while (getline(Fin,tmp))
  {
    unsigned pos = tmp.find(",");
    std::string ID = tmp.substr(0, pos);
    unsigned pos2 = tmp.find(",", pos + 1);
    std::string value = tmp.substr(pos + 1, pos2 - (pos + 1));
    unsigned tmpRFP = std::atoi(value.c_str());

    pos = tmp.find(",", pos2 + 1);
    value = tmp.substr(pos2 + 1, pos - (pos2 + 1));
    unsigned counts = std::atoi(value.c_str());

    std::string codon = tmp.substr(pos + 1, 3);
    for (unsigned i = 0; i < counts; i++)
      seq += codon;

    if (first)
    {
      prevID = ID;
      first = false;
    }
    if (ID != prevID)
    {
      tmpGene.setId(prevID);
      tmpGene.setDescription("No description for RFP Model");
      tmpGene.setSequence(seq);
      addGene(tmpGene); //add to genome
      tmpGene.clear();
      seq = "";
	}

    prevID = ID;
    unsigned index = SequenceSummary::CodonToIndex(codon);
    tmpGene.geneData.setRFPObserved(index, tmpRFP);
	}


  tmpGene.setId(prevID);
  tmpGene.setDescription("No description for RFP Model");
  tmpGene.setSequence(seq);

  addGene(tmpGene); //add to genome

  Fin.close();
}


void Genome::writeRFPFile(std::string filename, bool simulated)
{
	std::ofstream Fout;
	Fout.open(filename.c_str());
	if (Fout.fail())
	{
		std::cerr <<"Error in Genome::writeRFPFile: Cannot open output RFP file " << filename <<"\n";
	}

	Fout <<"ORF,RFP_Counts,Codon_Counts,Codon\n";
	for (unsigned geneIndex = 0; geneIndex < genes.size(); geneIndex++)
	{
		Gene *currentGene;
		if (simulated)
			currentGene = &simulatedGenes[geneIndex];
		else
			currentGene = &genes[geneIndex];


		for (unsigned codonIndex = 0; codonIndex < 64; codonIndex++)
		{
			std::string codon = SequenceSummary::codonArray[codonIndex];

			Fout << currentGene->getId() <<",";
			Fout << currentGene->geneData.getRFPObserved(codonIndex) <<",";
			Fout << currentGene->geneData.getCodonCountForCodon(codonIndex) <<"," << codon <<"\n";
		}
	}
	Fout.close();
}


Gene& Genome::getGene(unsigned index)
{
	Gene gene;

	if (index >= genes.size()) 
	{
		std::cerr << "Error in Genome::getGene: Index " << index << " is out of bounds.\n";
	}

	return index >= genes.size() ? gene : genes[index];
}
Gene& Genome::getGene(std::string id)
{
	unsigned i = 0;
	bool geneFound = false;
	while(!geneFound)
	{
		Gene tempGene = genes[i];
		geneFound = (tempGene.getId().compare(id) == 0);
		i++;
	}
	// i is increase after a potential finding, therefore i-1 is correct
	return genes[i-1];
}


void Genome::readObservedPhiValues(std::string filename, bool byId)
{
	std::ifstream input;
	std::string tmp;

	std::cout <<"Opening " << filename <<"\n";
	input.open(filename);
	if (input.fail())
	{
		std::cerr <<"Error opening file in readObservedPhiValues\n";
		std::exit(1);
	}



	if (genes.size() == 0)
	{
		std::cerr << "Genome is empty, function will not execute!\n";
	}
	else
	{
		if (byId)
		{
			//Mapping is done so the genes can be found
			std::map <std::string, Gene* > genomeMapping;
			for (unsigned i = 0; i < genes.size(); i++)
			{
				genomeMapping.insert(make_pair(genes[i].getId(), &genes[i]));
			}


			while (input >> tmp)
			{
				std::size_t pos = tmp.find(",");
				std::string geneID = tmp.substr(0, pos);
				std::map <std::string, Gene * >::iterator it;
				it = genomeMapping.find(geneID);
				if (it == genomeMapping.end()) {
					std::cerr << "Gene " << geneID << " not found!\n";
				}
				else //gene is found
				{

					bool notDone = true;
					while (notDone)
					{
						std::size_t pos2 = tmp.find(",", pos + 1);
						if (pos2 == std::string::npos)
						{
							notDone = false;
						}
						std::string val = tmp.substr(pos + 1, pos2 - (pos + 1));
						it->second -> observedPhiValues.push_back(std::atof(val.c_str())); //make vector private again
						pos = pos2;
					}
				}
			}


		} //end of putting in by ID
		else //doing this by index
		{
			unsigned geneIndex = 0;
			while (input >> tmp)
			{
				if (geneIndex >= genes.size())
				{
					std::cerr <<"GeneIndex exceeds the number of genes in the genome. Exiting function\n";
					break;
				}
				std::size_t pos = tmp.find(",");
				bool notDone = true;
				while (notDone)
				{
					std::size_t pos2 = tmp.find(",", pos + 1);
					if (pos2 == std::string::npos)
					{
						notDone = false;
					}
					std::string val = tmp.substr(pos + 1, pos2 - (pos + 1));
					genes[geneIndex].observedPhiValues.push_back(std::atof(val.c_str()));
					pos = pos2;
				}
				geneIndex++;
			}
		}//end of reading by index
	}
}




void Genome::clear()
{
	genes.clear();
	simulatedGenes.clear();
}

Genome Genome::getGenomeForGeneIndicies(std::vector <unsigned> indicies)
{
	Genome genome;

	for (unsigned i = 0; i < indicies.size(); i++)
	{
		genome.addGene(genes[indicies[i]]);
	}

	return genome;
}

Genome Genome::getGenomeForGeneIndiciesR(std::vector <unsigned> indicies)
{
	Genome genome;
	bool bad = false;
	for (unsigned i = 0; i < indicies.size(); i++)
	{
		if (indicies[i] == 0 || indicies[i] > genes.size())
		{
			std::cerr << "Problem with index " << indicies[i] << ".\n";
			std::cerr << "returning an empty genome\n";
			bad = true;
			break;
		}
		else
		{
			indicies[i]--;
		}
	}
	if (!bad)
	{
		genome = getGenomeForGeneIndicies(indicies);
	}

	return genome;

}
Gene& Genome::getGeneByIndex(unsigned index)
{
	Gene gene;
	bool checker = checkIndex(index, 1, genes.size());
	return checker ? genes[index - 1] : gene;
}

bool Genome::checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound)
{
	bool check = false;
	if (lowerbound <= index && index <= upperbound)
	{
		check = true;
	}
	else
	{
		std::cerr <<"Error with the index\nGIVEN: " << index <<"\n";
		std::cerr <<"MUST BE BETWEEN:	" << lowerbound << " & " << upperbound <<"\n";
	}
	return check;
}

// ---------------------------------------------------------------------------
// ----------------------------- RCPP STUFF ----------------------------------
// ---------------------------------------------------------------------------
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;

	RCPP_EXPOSED_CLASS(Gene)
RCPP_EXPOSED_CLASS(Genome)

RCPP_MODULE(Genome_mod)
{
	class_<Genome>("Genome")
		.constructor("empty constructor")

		.method("readFasta", &Genome::readFasta, "reads a genome into the object")
		.method("writeFasta", &Genome::writeFasta, "writes the genome to a fasta file")
		.method("readRFPFile", &Genome::readRFPFile, "reads RFP data in for the RFP model")
		.method("getGenomeSize", &Genome::getGenomeSize, "returns how many genes are in the genome")
		.method("clear", &Genome::clear, "clears the genome")
		.method("getCodonCountsPerGene", &Genome::getCodonCountsPerGene, "returns a vector of codon counts for a given gene")

		//R Wrapper function
		.method("getGeneByIndex", &Genome::getGeneByIndex, "returns a gene for a given index")
		.method("getGenomeForGeneIndicies", &Genome::getGenomeForGeneIndiciesR, "returns a new genome based on the ones requested in the given vector")
		;
}
#endif
