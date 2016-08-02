#include "include/Gene.h"

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif



//--------------------------------------------------//
// ---------- Constructors & Destructors ---------- //
//--------------------------------------------------//


/* Gene constructor (RCPP EXPOSED)
 * Arguments: None
 * Blank constructor for Gene class. Sets the sequence, id, and
 * description fields to empty strings.
*/
Gene::Gene() : seq(""), id(""), description("")
{
    //ctor
}


/* Gene constructor (RCPP EXPOSED)
 * Arguments: sequence, id, description
 * Gene constructor that will set the sequence, id, and description strings in the Gene class. It will
 * then generate the sequence summary object for the gene based off the set sequence as long as the sequence
 * is a multiple of three (three is needed because of the size of codons).
*/
Gene::Gene(std::string _seq, std::string _id, std::string _desc) : seq(_seq), id(_id), description(_desc)
{
	if (seq.length() % 3 == 0)
		geneData.processSequence(_seq);
	else
    {
        my_printError("WARNING: Gene: % has sequence length NOT multiple of 3 after cleaning of the sequence!\n", id);
        my_printError("Gene data is NOT processed!\nValid characters are A,C,T,G, and N \n");
    }
}


/* Gene copy constructor (NOT EXPOSED)
 * Arguments: Gene object
 * Copy constructor for gene. All fields, public and private, will be set for the
 * given gene.
*/
Gene::Gene(const Gene& other)
{
    seq = other.seq;
    id = other.id;
    description = other.description;
    geneData = other.geneData;
	observedSynthesisRateValues = other.observedSynthesisRateValues;
}


/* Gene = operator (NOT EXPOSED)
 * Arguments: Gene object
 * Overloaded definition of the assignment operator ( = ). Function is
 * similar to that of the copy constructor.
*/
Gene& Gene::operator=(const Gene& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    seq = rhs.seq;
    id = rhs.id;
    description = rhs.description;
    geneData = rhs.geneData;
    observedSynthesisRateValues = rhs.observedSynthesisRateValues;
    return *this;
}


/* Gene == operator (NOT EXPOSED)
 * Arguments: Gene object
 * Overloaded definition of the equality operator ( == ). Will compare all
 * fields, private and public, contained in gene. Returns false if one of the
 * comparisons fail.
*/
bool Gene::operator==(const Gene& other) const
{
    bool match = true;

    if (this->seq != other.seq) { match = false; }
    if (this->id != other.id) { match = false; }
    if (this->description != other.description) { match = false; }
    if (this->observedSynthesisRateValues != other.observedSynthesisRateValues) { match = false; }
    if (!(this->geneData == other.geneData)) { match = false; } //if structures aren't equal, genes aren't equal.

    return match;
}


/* Gene deconstructor (NOT EXPOSED)
 * Arguments: None
 * Standard deconstructor for the gene object.
*/
Gene::~Gene()
{
    //dtor
}





//-------------------------------------------------//
//---------- Data Manipulation Functions ----------//
//-------------------------------------------------//


/* getId (RCPP EXPOSED)
 * Arguments: None
 * Returns the gene's id.
*/
std::string Gene::getId()
{
    return id;
}


/* setId (RCPP EXPOSED)
 * Arguments: id
 * Takes the given id string and sets the gene's id field.
*/
void Gene::setId(std::string _id)
{
    id = _id;
}


/* getDescription (RCPP EXPOSED)
 * Arguments: None
 * Returns the gene's description.
*/
std::string Gene::getDescription()
{
    return description;
}


/* setDescription (RCPP EXPOSED)
 * Arguments: description
 * Takes the given description string and sets the
 * gene's description field.
*/
void Gene::setDescription(std::string _desc)
{
    description = _desc;
}


/* getSequence (RCPP EXPOSED)
 * Arguments: None
 * Returns the gene's sequence string.
*/
std::string Gene::getSequence()
{
    return seq;
}


/* setSequence (RCPP EXPOSED)
 * Arguments: sequence
 * Takes the specified sequence string, clearing the current sequence summary. Provided
 * that the new string length is a multiple of 3, the string
 * is processed and set. If not, it is not processed.
 * NOTE: The seq string will still be set, even if it is invalid.
 * NOTE: As part of changing the sequence, the sequence summary is also cleared.
*/
void Gene::setSequence(std::string _seq)
{
    geneData.clear();
    std::transform(_seq.begin(), _seq.end(), _seq.begin(), ::toupper);
    seq = _seq;
	if (seq.length() % 3 == 0)
	{
		bool check = geneData.processSequence(seq);
		if (!check)
			my_printError("WARNING: Error with gene %\nBad codons found!\n", id);
	}
	else
    {
		my_printError("WARNING: Gene: % has sequence length NOT multiple of 3 after cleaning of the sequence!\n");
        my_printError("Gene data is NOT processed!\nValid characters are A,C,T,G, and N \n", id);
    }
}


std::vector <unsigned> Gene::getRFP_count()
{
    return geneData.getRFP_count();
}


void Gene::addRFP_count(std::vector <unsigned> RFP_counts)
{
    geneData.setRFP_count(RFP_counts);
}


/* getSequenceSummary (NOT EXPOSED)
 * Arguments: None
 * Returns a pointer to the stored sequence summary object.
*/
SequenceSummary *Gene::getSequenceSummary()
{
    SequenceSummary *rv = &geneData;
    return rv;
}


/* getObservedSynthesisRateValues (RCPP EXPOSED)
 * Arguments: None
 * Returns the vector containing the set synthesis rate values
 * for the gene.
*/
std::vector<double> Gene::getObservedSynthesisRateValues()
{
    return observedSynthesisRateValues;
}


/* setObservedSynthesisRateValues (NOT EXPOSED)
 * Arguments: vector of synthesis rate values
 * Sets the argument vector for the gene's synthesis rate values.
*/
void Gene::setObservedSynthesisRateValues(std::vector <double> values)
{
    observedSynthesisRateValues = values;
}


/* getObservedSynthesisRate (NOT EXPOSED)
 * Arguments: index to the synthesis rate vector
 * Returns the value in the synthesis rate vector for the given index. NOTE: Could crash if the
 * given index is out of bounds.
*/
double Gene::getObservedSynthesisRate(unsigned index)
{
	return observedSynthesisRateValues[index];
}


/* getNumObservedSynthesisSets (NOT EXPOSED)
 * Arguments: None
 * Returns the number of values stored in the synthesis rate vector.
*/
unsigned Gene::getNumObservedSynthesisSets()
{
	return observedSynthesisRateValues.size();
}


/* getNucleotideAt (NOT EXPOSED)
 * Arguments: index of the sequence string
 * Returns the nucleotide at the given index in the seq string
 * in the gene. NOTE: could crash if the index is out of bounds.
*/
char Gene::getNucleotideAt(unsigned i)
{
    return seq[i];
}





//-------------------------------------//
//---------- Other Functions ----------//
//-------------------------------------//


/* clear (NOT EXPOSED)
 * Arguments: None
 * Clears all the gene's variables.
*/
void Gene::clear()
{
  seq = "";
  id = "";
  description = "";
  geneData.clear();
  observedSynthesisRateValues.clear();
}


/* length (RCPP EXPOSED)
 * Arguments: None
 * Returns the length of the sequence string (ie, the number of nucleotides).
*/
unsigned Gene::length()
{
    return (unsigned)seq.size();
}


/* reverseComplement (NOT EXPOSED)
 * Arguments: None
 * Returns a new gene that will differ only in the sequence string.
 * The sequence string will be reversed and will contain each nucleotide's
 * complement.
*/
Gene Gene::reverseComplement()
{
  Gene tmpGene;
  tmpGene.id = id;
  tmpGene.description = description;
  tmpGene.seq = seq;
  tmpGene.observedSynthesisRateValues = observedSynthesisRateValues;

  std::reverse_copy(seq.begin(), seq.end(), tmpGene.seq.begin());
  std::transform(tmpGene.seq.begin(), tmpGene.seq.end(), tmpGene.seq.begin(),
            SequenceSummary::complimentNucleotide);
  return tmpGene;
}


/* toAASequence (NOT EXPOSED)
 * Arguments: None
 * Returns a string of amino acids corresponding to the gene's
 * sequence string. The string is looked at as codons and the codons
 * are then mapped to their respective amino acid. NOTE: This could crash if
 * the stored seq string is not of length three.
*/
std::string Gene::toAASequence()
{
    std::string aaseq = "";
    for (unsigned i = 0; i < seq.length(); i+=3)
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
        my_print("Invalid string given. Returning 0.\n");
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
        my_print("Invalid codon given. Returning 0.\n");
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
        my_print("Invalid codon given. Returning 0.\n");
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
        my_print("Invalid codon given. Returning empty vector.\n");
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


