#include "include/Gene.h"


#include <iostream>


Gene::Gene() : seq(""), id(""), description("")
{
    //ctor
}


Gene::Gene(std::string _seq, std::string _id, std::string _desc) : seq(_seq), id(_id), description(_desc)
{
    cleanSeq();
	if (seq.length() % 3 == 0)
	{
		geneData.processSequence(seq);
	}
	else 
	{
		std::cerr << "Gene: " << id << " has sequence length NOT multiple of 3 after cleaning of the sequence!" <<
				"\nGene data is NOT processed! \nValid characters are A,C,T,G, and N \n";
	}
}

Gene::~Gene()
{
    //dtor
}


Gene::Gene(const Gene& other)
{
    seq = other.seq;
    id = other.id;
    description = other.description;
    geneData = other.geneData;
	observedPhiValues = other.observedPhiValues;
}


Gene& Gene::operator=(const Gene& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    seq = rhs.seq;
    id = rhs.id;
    description = rhs.description;
    geneData = rhs.geneData;
	observedPhiValues = rhs.observedPhiValues;
    //assignment operator
    return *this;
}


void Gene::cleanSeq()
{
    std::string valid = "ACGTN";
    for (unsigned i = 0; i < seq.length(); i++) {
        if (valid.find(seq[i]) == std::string::npos) {
            seq.erase(i);
        }
    }
}


std::string Gene::getId()
{
    return id;
}


void Gene::setId(std::string _id)
{
    id = _id;
}


std::string Gene::getDescription()
{
    return description;
}


void Gene::setDescription(std::string _desc)
{
    description = _desc;
}


std::string Gene::getSequence()
{
    return seq;
}


void Gene::setSequence(std::string _seq)
{
    std::transform(_seq.begin(), _seq.end(), _seq.begin(), ::toupper);
    seq = _seq;
    cleanSeq();
	if (seq.length() % 3 == 0)
	{
		bool check = geneData.processSequence(seq);
		if (!check)
		{
			std::cerr <<"Error with gene " << id <<"\nBad codons found!\n";
		}
	}
	else
	{
		std::cerr << "Gene: " << id << " has sequence length NOT multiple of 3 after cleaning of the sequence!" <<
				"\nGene data is NOT processed! \nValid characters are A,C,T,G, and N \n";
	}
}


SequenceSummary& Gene::getSequenceSummary()
{
    return geneData;
}


std::vector<double> Gene::getObservedPhiValues()
{
    return observedPhiValues;
}


char Gene::getNucleotideAt(unsigned i)
{
    return seq[i];
}


void Gene::clear()
{
  seq = "";
  id = "";
  description = "";
  geneData.clear();
}

unsigned Gene::length()
{
    return (unsigned)seq.size();
}


Gene Gene::reverseComplement()
{
  Gene tmpGene;
  tmpGene.id = id;
  tmpGene.description = description;
  tmpGene.seq = seq;

  std::reverse_copy(seq.begin(), seq.end(), tmpGene.seq.begin());
  std::transform(tmpGene.seq.begin(), tmpGene.seq.end(), tmpGene.seq.begin(),
            SequenceSummary::complimentNucleotide);
  return tmpGene;
}


std::string Gene::toAASequence()
{

    std::string aaseq = "";
    for(unsigned i = 0; i < seq.length(); i+=3)
    {
        std::string codon = seq.substr(i, 3);
        aaseq += SequenceSummary::codonToAA(codon);
    }
    return aaseq;
}


//---------------------R WRAPPER FUNCTIONS---------------------//

void Gene::cleanSeqR()
{
    cleanSeq();
}


unsigned Gene::getAACount(std::string aa)
{
    return geneData.getAACountForAAR(aa);
}


unsigned Gene::getCodonCount(std::string& codon)
{
    return geneData.getCodonCountForCodonR(codon);
}

void Gene::setRFPObserved(unsigned index, unsigned value)
{
	geneData.setRFPObserved(index, value);
}

unsigned Gene::getRFPObserved(std::string codon)
{
    return geneData.getRFPObservedForCodonR(codon);
}


std::vector <unsigned> Gene::getCodonPositions(std::string codon)
{
    return geneData.getCodonPositionsForCodonR(codon);
}



// ---------------------------------------------------------------------------
// ----------------------------- RCPP STUFF ----------------------------------
// ---------------------------------------------------------------------------
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;

RCPP_EXPOSED_CLASS(Gene) //Exposed because of functions that return a gene.
RCPP_EXPOSED_CLASS(SequenceSummary)

RCPP_MODULE(Gene_mod)
{
  class_<Gene>( "Gene" )
    
	.constructor("empty constructor")
    .constructor<std::string, std::string, std::string >("Initialize a gene by giving the id, description, and sequence string")

    //Private functions:
    .method("cleanSeq", &Gene::cleanSeqR) //TEST THAT ONLY!


	//Public functions & variables:
	.field("geneData", &Gene::geneData)
	.property("id", &Gene::getId, &Gene::setId)
    .property("description", &Gene::getDescription, &Gene::setDescription)
    .property("seq", &Gene::getSequence, &Gene::setSequence)

    .method("getSequenceSummary", &Gene::getSequenceSummary) //TEST THAT ONLY!
    .method("getObservedPhiValues", &Gene::getObservedPhiValues)
    .method("getNucleotideAt", &Gene::getNucleotideAt) //TEST THAT ONLY!
	.method("clear", &Gene::clear, "clears the id, sequence, and description in the object")
    .method("length", &Gene::length, "returns the length of sequence")
    .method("reverseComplement", &Gene::reverseComplement) //TEST THAT ONLY!
    .method("toAASequence", &Gene::toAASequence)


	.method("getAACount", &Gene::getAACount, "returns the number of amino acids that are in the sequence for a given amino acid")
	.method("getCodonCount", &Gene::getCodonCount, "returns the number of codons that are in the sequence for a given codon")
	.method("getRFPObserved", &Gene::getRFPObserved)
	.method("setRFPObserved", &Gene::setRFPObserved)
	.method("getCodonPositions", &Gene::getCodonPositions)
  ;
}
#endif


