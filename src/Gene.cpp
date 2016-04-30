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
		geneData.processSequence(_seq);
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
	observedSynthesisRateValues = other.observedSynthesisRateValues;
}


Gene& Gene::operator=(const Gene& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    seq = rhs.seq;
    id = rhs.id;
    description = rhs.description;
    geneData = rhs.geneData;
    observedSynthesisRateValues = rhs.observedSynthesisRateValues;
    //assignment operator
    return *this;
}


bool Gene::operator==(const Gene& other) const
{
    bool match = true;


    if(this->seq != other.seq) { match = false; std::cerr<<"faila\n";}
    if(this->id != other.id) { match = false; std::cerr<<"failb\n";}
    if(this->description != other.description) { match = false; std::cerr<<"failc\n";}
    if(this->observedSynthesisRateValues != other.observedSynthesisRateValues) { match = false; std::cerr<<"faild\n";}
    if(!(this->geneData == other.geneData)) { match = false; std::cerr<<"faile\n";} //if structures aren't equal, genes aren't equal.

    return match;
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


SequenceSummary *Gene::getSequenceSummary()
{
    SequenceSummary *rv = &geneData;
    return rv;
}


std::vector<double> Gene::getObservedSynthesisRateValues()
{
    return observedSynthesisRateValues;
}


void Gene::setObservedSynthesisRateValues(std::vector <double> values)
{
    observedSynthesisRateValues = values;
}

double Gene::getObservedSynthesisRate(unsigned index)
{
	return observedSynthesisRateValues[index];
}

unsigned Gene::getNumObservedSynthesisSets()
{
	return observedSynthesisRateValues.size();
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
    unsigned rv = 0;

    if (SequenceSummary::aaToIndex.end() != SequenceSummary::aaToIndex.find(aa))
    {
        rv = geneData.getAACountForAA(aa);
    }
    else
    {
        Rprintf("Invalid string given. Returning 0.\n");
    }
    return rv;
}


unsigned Gene::getCodonCount(std::string& codon)
{
    unsigned rv = 0;

    if (SequenceSummary::codonToIndexWithReference.end() != SequenceSummary::codonToIndexWithReference.find(codon))
    {
        rv = geneData.getCodonCountForCodon(codon);
    }
    else
    {
        Rprintf("Invalid codon given. Returning 0.\n");
    }
    return rv;
}

unsigned Gene::getRFPObserved(std::string codon)
{
    unsigned rv = 0;

    if (SequenceSummary::codonToIndexWithReference.end() != SequenceSummary::codonToIndexWithReference.find(codon))
    {
        rv = geneData.getRFPObserved(codon);
    }
    else
    {
        Rprintf("Invalid codon given. Returning 0.\n");
    }
    return rv;
}


std::vector <unsigned> Gene::getCodonPositions(std::string codon)
{
    std::vector <unsigned> rv;
    std::vector <unsigned> *tmp;
    tmp = &rv; //So if an invalid codon is given, tmp will point to an empty vector.


    if (SequenceSummary::codonToIndexWithReference.end() != SequenceSummary::codonToIndexWithReference.find(codon))
    {
        tmp = geneData.getCodonPositions(codon);
    }
    else
    {
        Rprintf("Invalid codon given. Returning empty vector.\n");
    }


    for (unsigned i = 0; i < tmp -> size(); i++)
    {
        rv.push_back(tmp->at(i));
    }
    return rv;
}




//---------------------------------//
//---------- RCPP Module ----------//
//---------------------------------//


RCPP_MODULE(Gene_mod)
{
  class_<Gene>( "Gene" )
    
	.constructor("empty constructor")
    .constructor<std::string, std::string, std::string >("Initialize a gene by giving the id, description, and sequence string")


	//Public functions & variables:
	.property("id", &Gene::getId, &Gene::setId)
    .property("description", &Gene::getDescription, &Gene::setDescription)
    .property("seq", &Gene::getSequence, &Gene::setSequence)

    .method("getObservedSynthesisRateValues", &Gene::getObservedSynthesisRateValues)
    .method("length", &Gene::length, "returns the length of sequence")


	.method("getAACount", &Gene::getAACount, "returns the number of amino acids that are in the sequence for a given amino acid")
	.method("getCodonCount", &Gene::getCodonCount, "returns the number of codons that are in the sequence for a given codon")
	.method("getRFPObserved", &Gene::getRFPObserved)
	.method("getCodonPositions", &Gene::getCodonPositions)
  ;
}
#endif


