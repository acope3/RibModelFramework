#include "include/Gene.h"
#include "include/SequenceSummary.h"

#include <algorithm>
#include <iostream>
Gene::Gene() : seq(""), id(""), description("")
{

}

Gene::Gene(std::string _seq, std::string _id, std::string _desc) : seq(_seq), id(_id), description(_desc)
{
    cleanSeq();
    geneData.processSequence(seq);
}

Gene::Gene(const Gene& other)
{
    seq = other.seq;
    id = other.id;
    description = other.description;
    geneData = other.geneData;
}

Gene& Gene::operator=(const Gene& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    seq = rhs.seq;
    id = rhs.id;
    description = rhs.description;
    geneData = rhs.geneData;
    //assignment operator
    return *this;
}

Gene::~Gene()
{
    //dtor
}


void Gene::setSequence(std::string _seq)
{
    std::transform(_seq.begin(), _seq.end(), _seq.begin(), ::toupper);
    seq = _seq;
    cleanSeq();
    geneData.processSequence(seq);
}

void Gene::clear()
{
  seq = "";
  id = "";
  description = "";
  geneData.clear();
}

bool notNucleotide(char ch)
{
    std::string str = "ACGT";
    std::string ts(1,ch);
    return ( str.find(ts) == std::string::npos );
}

void Gene::cleanSeq()
{
  std::string::iterator endp = std::remove_if(seq.begin(), seq.end(), notNucleotide);
  seq.erase(endp, seq.end());
}

char complimentNucleotide(char ch)
{
  //std::string str = "ACGT";
  std::string ts(1,ch);
  if( ch == 'A' ) return 'T';
  else if( ch == 'T' ) return 'A';
  else if( ch == 'C' ) return 'G';
  else return 'C';
}

Gene Gene::reverseCompliment()
{
  Gene tmpGene;
  tmpGene.id = id;
  tmpGene.description = description;
  tmpGene.seq = seq;

  std::reverse_copy(seq.begin(), seq.end(), tmpGene.seq.begin());
  std::transform(tmpGene.seq.begin(), tmpGene.seq.end(), tmpGene.seq.begin(),
            complimentNucleotide);
  return tmpGene;
}


std::string Gene::toAAsequence()
{

    std::string aaseq = "";
    for(unsigned i = 0; i < seq.length(); i+=3)
    {
        std::string codon = seq.substr(i, 3);
        //std::cout << codon << ": " << Gene::CodonToAA(codon) << "\n" << std::endl;
        aaseq.append(1, SequenceSummary::CodonToAA(codon));
    }
    return aaseq;
}



// ---------------------------------------------------------------------------
// ----------------------------- RCPP STUFF ----------------------------------
// ---------------------------------------------------------------------------
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;

RCPP_EXPOSED_CLASS(Gene)

RCPP_MODULE(Gene_mod)
{
  class_<Gene>( "Gene" )
    
		.constructor("empty constructor")
    .constructor<std::string, std::string, std::string >("Initialize a gene by giving the id, description, and sequence string")

    .property("id", &Gene::getId, &Gene::setId)
    .property("description", &Gene::getDescription, &Gene::setDescription)
    .property("seq", &Gene::getSequence, &Gene::setSequence)

		.method("getAACount", &Gene::getAACount)
		.method("getCodonCount", &Gene::getCodonCount)
		.method("clear", &Gene::clear)
    .method("length", &Gene::length)
  ;
}
#endif

