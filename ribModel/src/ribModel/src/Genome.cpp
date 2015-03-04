#include "../include/Genome.h"


#include <iostream>     // std::cout
#include <cstring>
#include <fstream>

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


void Genome::writeFasta ( char* filename )
{
    try {
        std::ofstream Fout(filename);
        if(!Fout) throw strcat("Cannot open output Fasta file ", filename);

        for(int i = 0; i < genes.size(); i++) {
            Fout << ">" << genes[i].getDescription() << std::endl;
            for(int j = 0; j < genes[i].length(); j++) {
                Fout << genes[i].getNucleotideAt(j);
                if((j + 1) % 60 == 0) Fout << std::endl;
            }
            Fout << std::endl;
        }
    }

    catch(char* pMsg) { std::cerr << std::endl << "Exception:" << pMsg << std::endl; }
}

void Genome::readFasta(char* filename) // read Fasta format sequences
{
    try
    {
        std::ifstream Fin(filename);
        if(!Fin) throw strcat("Genome::readFasta throws: Cannot open input Fasta file ", filename);

        bool fastaFormat = false;
        std::string buf;
        int newLine;


        Gene tmpGene;
        std::string tempSeq = "";
        for(;;)
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



