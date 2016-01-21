#include "include/Gene.h"

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif



//--------------------------------------------------//
// ---------- Constructors & Destructors ---------- //
//--------------------------------------------------//


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
#ifndef STANDALONE
		Rf_warning("Gene: %s has sequence length NOT multiple of 3 after cleaning of the sequence!\nGene data is NOT processed! \nValid characters are A,C,T,G, and N \n", id.c_str());
#else
		std::cerr << "Gene: " << id << " has sequence length NOT multiple of 3 after cleaning of the sequence!" <<
				"\nGene data is NOT processed! \nValid characters are A,C,T,G, and N \n";
#endif
	}
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


Gene::~Gene()
{
    //dtor
}





//-------------------------------------------------//
//---------- Data Manipulation Functions ----------//
//-------------------------------------------------//


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
#ifndef STANDALONE
			Rf_warning("Error with gene %s\nBad codons found!\n", id.c_str());
#else
			std::cerr << "Error with gene " << id << "\nBad codons found!\n";
#endif
		}
	}
	else
	{
#ifndef STANDALONE
		Rf_warning("Gene: %s has sequence length NOT multiple of 3 after cleaning of the sequence!\nGene data is NOT processed! \nValid characters are A,C,T,G, and N \n", id.c_str());
#else
		std::cerr << "Gene: " << id << " has sequence length NOT multiple of 3 after cleaning of the sequence!" <<
				"\nGene data is NOT processed! \nValid characters are A,C,T,G, and N \n";
#endif
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


void Gene::setObservedPhiValues(std::vector <double> values)
{
    observedPhiValues = values;
}

double Gene::getObservedSynthesisRate(unsigned index)
{
	return observedPhiValues[index];
}

unsigned Gene::getNumObservedSynthesisSets()
{
	return observedPhiValues.size();
}

char Gene::getNucleotideAt(unsigned i)
{
    return seq[i];
}





//-------------------------------------//
//---------- Other Functions ----------//
//-------------------------------------//


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







// -----------------------------------------------------------------------------------------------------//
// ---------------------------------------- R SECTION --------------------------------------------------//
// -----------------------------------------------------------------------------------------------------//



#ifndef STANDALONE

unsigned Gene::getAACount(std::string aa)
{
    bool error = false;
    unsigned rv = 0;

    //TODO: more extranious testing on input (capital letters, valid letters) to prevent R crashing (on all fcts below as well).
    if (aa.size() != 1)
    {
        error = true;
        std::cerr <<"Invalid string given. Returning 0.\n";
    }
    if (!error)
    {
        rv = geneData.getAACountForAA(aa);
    }
    return rv;
}


unsigned Gene::getCodonCount(std::string& codon)
{
    bool error = false;
    unsigned rv = 0;
    if (codon.size() != 3)
    {
        error = true;
        std::cerr <<"Invalid codon given. Returning 0.\n";
    }

    if(!error)
    {
        rv = geneData.getCodonCountForCodon(codon);
    }
    return rv;
}

unsigned Gene::getRFPObserved(std::string codon)
{
    return geneData.getRFPObserved(codon);
}


std::vector <unsigned> Gene::getCodonPositions(std::string codon)
{
    std::vector <unsigned> rv;
    std::vector <unsigned> *tmp;
    tmp = geneData.getCodonPositions(codon);
    for (unsigned i = 0; i < tmp -> size(); i++)
    {
        rv.push_back(tmp->at(i));
    }
    return rv;
}




//---------------------------------//
//---------- RCPP Module ----------//
//---------------------------------//


RCPP_EXPOSED_CLASS(SequenceSummary) //TODO: try removing this at the end now.

RCPP_MODULE(Gene_mod)
{
  class_<Gene>( "Gene" )
    
	.constructor("empty constructor")
    .constructor<std::string, std::string, std::string >("Initialize a gene by giving the id, description, and sequence string")


	//Public functions & variables:
	.property("id", &Gene::getId, &Gene::setId)
    .property("description", &Gene::getDescription, &Gene::setDescription)
    .property("seq", &Gene::getSequence, &Gene::setSequence)

    .method("getObservedPhiValues", &Gene::getObservedPhiValues)
    .method("length", &Gene::length, "returns the length of sequence")


	.method("getAACount", &Gene::getAACount, "returns the number of amino acids that are in the sequence for a given amino acid")
	.method("getCodonCount", &Gene::getCodonCount, "returns the number of codons that are in the sequence for a given codon")
	.method("getRFPObserved", &Gene::getRFPObserved)
	.method("getCodonPositions", &Gene::getCodonPositions)
  ;
}
#endif


