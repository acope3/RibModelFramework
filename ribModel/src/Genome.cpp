#include "include/Genome.h"


#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>


Genome::Genome()
{
	//ctor
	numGenesWithPhi = 0;
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
	numGenesWithPhi = rhs.numGenesWithPhi;
	//assignment operator
	return *this;
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
					if(!fastaFormat)
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
					std::size_t pos = buf.find(" ") - 1;
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
					Fout << ">" << simulatedGenes[i].getDescription() <<"\n";
					for(unsigned j = 0u; j < simulatedGenes[i].length(); j++)
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
					Fout << ">" << genes[i].getDescription() << std::endl;
					for(unsigned j = 0u; j < genes[i].length(); j++)
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
	bool first = true;
	std::string seq = "";

	while (getline(Fin,tmp))
	{
		std::size_t pos = tmp.find(",");
		std::string ID = tmp.substr(0, pos);


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
		std::size_t pos2 = tmp.find(",", pos + 1);
		std::string value = tmp.substr(pos + 1, pos2 - (pos + 1));
		unsigned tmpRFP = (unsigned)std::atoi(value.c_str());
		pos = tmp.find(",", pos2 + 1);
		value = tmp.substr(pos2 + 1, pos - (pos2 + 1));
		unsigned counts = (unsigned)std::atoi(value.c_str());

		std::string codon = tmp.substr(pos + 1, 3);
		for (unsigned i = 0; i < counts; i++)
			seq += codon;

		prevID = ID;
		unsigned index = SequenceSummary::codonToIndex(codon);
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


void Genome::readObservedPhiValues(std::string filename, bool byId)
{
	std::cout << numGenesWithPhi <<"\n";
	std::ifstream input;
	std::string tmp;
	unsigned numPhi = 0;
	bool exitfunction = false;

	input.open(filename);
	if (input.fail())
	{
		std::cerr <<"Error opening file in readObservedPhiValues\n";
		std::exit(1);
	}

	std::getline(input, tmp); //Trash the header line

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

			bool first = true;
			while (std::getline(input, tmp))
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
					unsigned count = 0;
					std::string val ="";
					bool notDone = true;
					while (notDone)
					{
						std::size_t pos2 = tmp.find(",", pos + 1);
						if (pos2 == std::string::npos)
						{
							val = tmp.substr(pos + 1, tmp.size() + 1);
							notDone = false;
							if (first)
							{
								first = false;
								numPhi = (unsigned) it -> second -> observedPhiValues.size() + 1;
							}
							else
							{
								if (count + 1 != numPhi)
								{
									std::cerr << geneID <<": has a differnt number of phi values given than other genes: ";
									std::cerr << it-> second -> getId() <<"\n";
									std::cerr << it-> second -> observedPhiValues.size() + 1 <<". Exiting function.\n";
									exitfunction = true;
									for (unsigned a = 0; a < getGenomeSize(); a++)
									{
										genes[a].observedPhiValues.clear();
									}
									break;
								}
							}
						}
						else {
							val = tmp.substr(pos + 1, pos2 - (pos + 1));
						}
						double value = std::atof(val.c_str());
						if (value <=  0 || std::isnan(value))
						{
							if (value == 0 || std::isnan(value)) value = -1;
							else
							{
								std::cerr <<"WARNING! Negative phi value given - values should not be on the log scale. Negative Value stored.";
							}
						}

						it->second->observedPhiValues.push_back(value); //make vector private again
						if (it->second->observedPhiValues.size() == 1) numGenesWithPhi++;
						pos = pos2;
						count++;
					}
				}
				if (exitfunction) break;
			}
			if (!exitfunction)
			{
				for (unsigned geneIndex = 0; geneIndex < getGenomeSize(); geneIndex++)
				{
					if (getGene(geneIndex).observedPhiValues.size() != numPhi)
					{
						Gene *gene = &(getGene(geneIndex));
						std::cerr <<"Gene # " << geneIndex <<" (" << gene->getId() <<") does not have any phi values.";
						std::cerr <<" Filling with -1's\n";
						gene->observedPhiValues.resize(numPhi, -1);
					}
				}
			}

		} //end of putting in by ID
		else //doing this by index
		{
			unsigned geneIndex = 0;
			bool first = true;
			while (std::getline(input, tmp))
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
						if (first)
						{
							first = false;
							numPhi = (unsigned) genes[geneIndex].observedPhiValues.size() + 1;
						}
						notDone = false;
						if (numPhi != genes[geneIndex].observedPhiValues.size() + 1)
						{
							std::cerr <<"gene " << geneIndex <<": has a differnt number of phi values given than other genes, exiting.\n";
							for (unsigned a = 0; a < getGenomeSize(); a++)
							{
								genes[a].observedPhiValues.clear();
							}
							exitfunction = true;
							break;
						}
					}

					std::string val = tmp.substr(pos + 1, pos2 - (pos + 1));

					double value = std::atof(val.c_str());
					if (value <=  0 || std::isnan(value))
					{
						if (value == 0 || std::isnan(value)) value = -1;
						else
						{
							std::cerr <<"WARNING! Negative phi value given - values should not be on the log scale. Negative Value stored.";
						}
					}
					genes[geneIndex].observedPhiValues.push_back(value);
					if (genes[geneIndex].observedPhiValues.size() == 1) numGenesWithPhi++;
					pos = pos2;
				}
				geneIndex++;
				if (exitfunction)
				{
					break;
				}
			}
			if (!exitfunction)
			{
				if (numGenesWithPhi != getGenomeSize())
				{
					std::cerr <<"The last " << getGenomeSize() - numGenesWithPhi <<" genes do not have phi values.";
					std::cerr <<"Please check your file to make sure every gene has a phi value. Filling empty genes";
					std::cerr <<"with -1's for calculations.\n";

					for (unsigned a = numGenesWithPhi; a < getGenomeSize(); a++)
					{
						Gene *gene = &(getGene(a));
						gene->observedPhiValues.resize(numPhi, -1);
					}
				}
			}
		}//end of reading by index
		input.close();
	}
}


void Genome::addGene(const Gene& gene, bool simulated)
{
	if (!simulated)
		genes.push_back(gene);
	else
		simulatedGenes.push_back(gene);
}


std::vector <Gene> Genome::getGenes(bool simulated)
{
	return  !simulated ? genes : simulatedGenes;
}


unsigned Genome::getNumGenesWithPhi()
{
	return numGenesWithPhi;
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


unsigned Genome::getGenomeSize()
{
	return (unsigned)genes.size();
}


void Genome::clear()
{
	genes.clear();
	simulatedGenes.clear();
	numGenesWithPhi = 0;
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


std::vector<unsigned> Genome::getCodonCountsPerGene(std::string codon)
{
	std::vector<unsigned> codonCounts(genes.size());
	unsigned codonIndex = SequenceSummary::codonToIndex(codon);
	for(unsigned i = 0u; i < genes.size(); i++)
	{
		Gene gene = genes[i];
		SequenceSummary seqsum = gene.getSequenceSummary();
		codonCounts[i] = seqsum.getCodonCountForCodon(codonIndex);
	}
	return codonCounts;
}


//---------------------R WRAPPER FUNCTIONS---------------------//

Gene& Genome::getGeneByIndex(unsigned index)
{
	Gene gene;
	bool checker = checkIndex(index, 1, (unsigned)genes.size());
	return checker ? genes[index - 1] : gene;
}


Gene& Genome::getGeneById(std::string ID)
{
	return getGene(ID);
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
		.method("writeRFPFile", &Genome::writeRFPFile)
		.method("readObservedPhiValues", &Genome::readObservedPhiValues)

		.method("addGene", &Genome::addGene) //TEST THAT ONLY!
		.method("getGenes", &Genome::getGenes) //TEST THAT ONLY!
		.method("getNumGenesWithPhi", &Genome::getNumGenesWithPhi) //TEST THAT ONLY!

		.method("checkIndex", &Genome::checkIndex) //TEST THAT ONLY!
		.method("getGenomeSize", &Genome::getGenomeSize, "returns how many genes are in the genome")
		.method("clear", &Genome::clear, "clears the genome")
		.method("getCodonCountsPerGene", &Genome::getCodonCountsPerGene, "returns a vector of codon counts for a given gene")

		//R Wrapper function
		.method("getGeneByIndex", &Genome::getGeneByIndex, "returns a gene for a given index")
		.method("getGeneById", &Genome::getGeneById) //TEST THAT ONLY!
		.method("getGenomeForGeneIndicies", &Genome::getGenomeForGeneIndiciesR, "returns a new genome based on the ones requested in the given vector")
		;
}
#endif
