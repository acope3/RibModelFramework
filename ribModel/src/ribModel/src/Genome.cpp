#include "../include/Genome.h"
#include "../include/ROCParameter.h"//these two files must be included here to get at the implimentation
#include "../include/ROCModel.h" 		//for simulateGenome. They cannot be included in the header file. See genome.h for
																		//more information on circular dependices/forward declarations
#include <iostream>     // std::cout
#include <cstring>
#include <fstream>
#include <ctime>

Genome::Genome()
{
	//ctor
}

Genome::~Genome()
{
	//dtor
}

Genome::Genome(const Genome& other)
{
	//copy ctor
	genes = other.genes;
}

Genome& Genome::operator=(const Genome& rhs)
{
	if (this == &rhs) return *this; // handle self assignment
	genes = rhs.genes;
	//assignment operator
	return *this;
}

void Genome::addGene(const Gene& gene)
{
	genes.push_back(gene);
}


void Genome::getCountsForAA(char aa, unsigned codonCounts[][5])
{
	for(int i = 0; i < genes.size(); i++)
	{
		Gene gene = genes[i];
		SequenceSummary seqsum = gene.getSequenceSummary();
		unsigned codonRange[2];
		SequenceSummary::AAToCodonRange(aa, false, codonRange);
		// get codon counts for AA
		unsigned j = 0u;
		for(unsigned k = codonRange[0]; k < codonRange[1]; k++, j++)
		{
			codonCounts[i][j] = seqsum.getCodonCountForCodon(k);
		}
	}
}

void Genome::writeFasta (char* filename, bool simulated)
{
	try {
		std::ofstream Fout;
		Fout.open(filename);
		if (Fout.fail())
		{
			std::cerr <<"Cannot open output Fasta file " << filename <<"\n";
			exit(1);
		}
		if (simulated)
		{
			for(int i = 0; i < simulatedGenes.size(); i++)
			{
				Fout << ">" << simulatedGenes[i].getId() <<" " << simulatedGenes[i].getDescription() <<"\n";
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
			for(int i = 0; i < genes.size(); i++)
			{
				Fout << ">" << genes[i].getId() <<" " << genes[i].getDescription() << std::endl;
				for(int j = 0; j < genes[i].length(); j++)
				{
					Fout << genes[i].getNucleotideAt(j);
					if((j + 1) % 60 == 0) Fout << std::endl;
				}
				Fout << std::endl;
			}
		}
		Fout.close();
	}
	catch(char* pMsg) { std::cerr << std::endl << "Exception:" << pMsg << std::endl; }
}

void Genome::readFasta(char* filename) // read Fasta format sequences
{
	try
	{
		std::ifstream Fin;
		Fin.open(filename);
		if (Fin.fail())
		{
			std::cerr <<"Genome::readFasta throws: Cannot open input Fasta file " << filename <<"\n";
			std::exit(1);
		}

		bool fastaFormat = false;
		std::string buf;
		int newLine;


		Gene tmpGene;
		std::string tempSeq = "";
		while (1)
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

					/*
						 std::cout << tmpGene.getId() << std::endl;
						 for(int i = 0; i < 22; i++)
						 {
						 std::cout << SequenceSummary::IndexToAA(i) << ":"<< tmpGene.geneData.getAAcountForAA(i) << "\t";
						 }
						 std::cout << std::endl << std::endl;
					 */
					tmpGene.clear();
					tempSeq = "";
				}
				tmpGene.setDescription( buf.substr(1,buf.size()-1) );
				int pos = buf.find(" ");
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
					/*
						 std::cout << tmpGene.getId() << std::endl;
						 for(int i = 0; i < 22; i++)
						 {
						 std::cout << SequenceSummary::IndexToAA(i) << ":"<< tmpGene.geneData.getAAcountForAA(i) << "\t";
						 }
						 std::cout << std::endl << std::endl;
					 */
					break;
				}
			}
		}
	}
	catch(char* pMsg) { std::cerr << std::endl << "Exception:" << pMsg << std::endl; }
}



Gene& Genome::getGene(int index)
{
	return genes[index];
}
Gene& Genome::getGene(std::string id)
{
	int i = 0;
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

void Genome::simulateGenome(ROCParameter& parameter, ROCModel& model)
{
	int i, j, k;
	int aaCount;
	int numCodons;
	unsigned codonIndex;
	unsigned aaRange[2];
	std::string tmpSeq;
	std::string codon;
	char curAA;
	std::string tmpDesc;
	int numParam = parameter.getNumParam();

	std::srand(std::time(0));
	simulatedGenes.resize(genes.size());
	tmpDesc = "Simulated Gene";


	for (i = 0; i < genes.size(); i++) //loop over all genes in the genome
	{
		Gene gene = genes[i];
		SequenceSummary seqSum = gene.geneData;
		tmpSeq = ""; //reset the sequence to blank
		tmpSeq += "ATG"; //Always will have the start amino acid


		unsigned mixtureElement = parameter.getMixtureAssignment(i);
		unsigned mutationCategory = parameter.getMutationCategory(mixtureElement);
		unsigned selectionCategory = parameter.getSelectionCategory(mixtureElement);
		unsigned expressionCategory = parameter.getExpressionCategory(mixtureElement);
		double phi = parameter.getExpression(i, expressionCategory, false);
		
		std::string tmpID = gene.getId() + "_MixtureElement" + std::to_string(mixtureElement);
		for (j = 0; j < 22; j++) //loop over each amino acid, naa[]
		{
			aaCount = seqSum.getAAcountForAA(j);
			curAA = seqSum.AminoAcidArray[j];
			if (curAA == 'X') continue;
			numCodons = seqSum.GetNumCodonsForAA(curAA);
			if (curAA == 'M') aaCount -= 1;
			double codonProb[numCodons]; //size the arrays to the proper size based on # of codons.
			double mutation[numCodons - 1];
			double selection[numCodons - 1];

			//get the probability vector for each amino acid
			if (curAA == 'M' || curAA == 'W')
			{
				codonProb[0] = 1;
			}
			else
			{
				parameter.getParameterForCategory(mutationCategory, ROCParameter::dM, curAA, false, mutation);
				parameter.getParameterForCategory(selectionCategory, ROCParameter::dEta, curAA, false, selection);
				model.calculateCodonProbabilityVector(numCodons, mutation, selection, phi, codonProb, numParam);
			}
			for (k = 0; k < aaCount; k++)
			{
				//std::cout <<"\t\taaCount is " << aaCount <<"\n";
				codonIndex = ROCParameter::randMultinom(codonProb, numCodons);
				seqSum.AAToCodonRange(curAA, false, aaRange); //need the first spot in the array where the codons for curAA are
				codon = seqSum.IndexToCodon(aaRange[0] + codonIndex);//get the correct codon based off codonIndex
				tmpSeq += codon;
			}
			//std::cout <<"\tDone with amino acid " << curAA <<"\n";
		}

		codon =	seqSum.IndexToCodon((rand() % 3) + 61); //randomly choose a stop codon, from range 61-63
		tmpSeq += codon;
		Gene tmpGene(tmpSeq, tmpID, tmpDesc);
		simulatedGenes[i] = tmpGene;
		//std::cout <<"done with gene " << i <<"\n";
	}
}

