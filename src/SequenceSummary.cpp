#include "include/SequenceSummary.h"

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

const std::string SequenceSummary::Ser2 = "Z";

const std::vector<std::string> SequenceSummary::AminoAcidArray = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
	"M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", SequenceSummary::Ser2, "X"};

const std::string SequenceSummary::codonArray[] =
		{"GCA", "GCC", "GCG", "GCT", "TGC", "TGT", "GAC", "GAT", "GAA", "GAG",
		 "TTC", "TTT", "GGA", "GGC", "GGG", "GGT", "CAC", "CAT", "ATA", "ATC",
		 "ATT", "AAA", "AAG", "CTA", "CTC", "CTG", "CTT", "TTA", "TTG", "ATG",
		 "AAC", "AAT", "CCA", "CCC", "CCG", "CCT", "CAA", "CAG", "AGA", "AGG",
		 "CGA", "CGC", "CGG", "CGT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC",
		 "ACG", "ACT", "GTA", "GTC", "GTG", "GTT", "TGG", "TAC", "TAT", "AGC",
		 "AGT", "TAA", "TAG", "TGA"};

const std::string SequenceSummary::codonArrayParameter[] =
		{"GCA", "GCC", "GCG", "TGC", "GAC",
		 "GAA", "TTC", "GGA", "GGC", "GGG",
		 "CAC", "ATA", "ATC", "AAA", "CTA",
		 "CTC", "CTG", "CTT", "TTA", "AAC",
		 "CCA", "CCC", "CCG", "CAA", "AGA",
		 "AGG", "CGA", "CGC", "CGG", "TCA",
		 "TCC", "TCG", "ACA", "ACC", "ACG",
		 "GTA", "GTC", "GTG", "TAC", "AGC"};



const std::map<std::string, unsigned> SequenceSummary::aaToIndex = {{"A", 0}, {"C", 1}, {"D", 2}, {"E", 3}, {"F", 4},
	{"G", 5}, {"H", 6}, {"I", 7}, {"K", 8}, {"L", 9}, {"M", 10}, {"N", 11}, {"P", 12}, {"Q", 13}, {"R", 14}, {"S", 15},
	{"T", 16}, {"V", 17}, {"W", 18}, {"Y", 19}, {SequenceSummary::Ser2, 20}, {"X", 21}};

const std::map<std::string, unsigned> SequenceSummary::codonToIndexWithReference = {{"GCA", 0}, {"GCC", 1}, {"GCG", 2},
	{"GCT", 3}, {"TGC", 4}, {"TGT", 5}, {"GAC", 6}, {"GAT", 7}, {"GAA", 8}, {"GAG", 9}, {"TTC", 10}, {"TTT", 11},
	{"GGA", 12}, {"GGC", 13}, {"GGG", 14}, {"GGT", 15}, {"CAC", 16}, {"CAT", 17}, {"ATA", 18}, {"ATC", 19}, {"ATT", 20},
	{"AAA", 21}, {"AAG", 22}, {"CTA", 23}, {"CTC", 24}, {"CTG", 25}, {"CTT", 26}, {"TTA", 27}, {"TTG", 28}, {"ATG", 29},
	{"AAC", 30}, {"AAT", 31}, {"CCA", 32}, {"CCC", 33}, {"CCG", 34}, {"CCT", 35}, {"CAA", 36}, {"CAG", 37}, {"AGA", 38},
	{"AGG", 39}, {"CGA", 40}, {"CGC", 41}, {"CGG", 42}, {"CGT", 43}, {"TCA", 44}, {"TCC", 45}, {"TCG", 46}, {"TCT", 47},
	{"ACA", 48}, {"ACC", 49}, {"ACG", 50}, {"ACT", 51}, {"GTA", 52}, {"GTC", 53}, {"GTG", 54}, {"GTT", 55}, {"TGG", 56},
	{"TAC", 57}, {"TAT", 58}, {"AGC", 59}, {"AGT", 60}, {"TAA", 61}, {"TAG", 62}, {"TGA", 63}};

const std::map<std::string, unsigned> SequenceSummary::codonToIndexWithoutReference = {{"GCA", 0}, {"GCC", 1},
	{"GCG", 2}, {"TGC", 3}, {"GAC", 4}, {"GAA", 5}, {"TTC", 6}, {"GGA", 7}, {"GGC", 8}, {"GGG", 9}, {"CAC", 10},
	{"ATA", 11}, {"ATC", 12}, {"AAA", 13}, {"CTA", 14}, {"CTC", 15}, {"CTG", 16}, {"CTT", 17}, {"TTA", 18}, {"AAC", 19},
	{"CCA", 20}, {"CCC", 21}, {"CCG", 22}, {"CAA", 23}, {"AGA", 24}, {"AGG", 25}, {"CGA", 26}, {"CGC", 27}, {"CGG", 28},
	{"TCA", 29}, {"TCC", 30}, {"TCG", 31}, {"ACA", 32}, {"ACC", 33}, {"ACG", 34}, {"GTA", 35}, {"GTC", 36}, {"GTG", 37},
	{"TAC", 38}, {"AGC", 39}}; //Need to add case for AGT


std::map<std::string, std::vector<std::string>> SequenceSummary::AAToCodonsMapWithReference = {{"A",{"GCA","GCC","GCG","GCT"}}, 
			{"C",{"TGC", "TGT"}}, {"D",{"GAC", "GAT"}}, {"E",{"GAA", "GAG"}}, {"F",{"TTC", "TTT"}},
			{"G",{"GGA", "GGC", "GGG", "GGT"}}, {"H",{"CAC", "CAT"}}, {"I",{"ATA", "ATC", "ATT"}},
			{"K",{"AAA", "AAG"}}, {"L",{"CTA", "CTC", "CTG","CTT", "TTA", "TTG"}}, {"M", {"ATG"}},
			{"N",{"AAC", "AAT"}}, {"P",{"CCA", "CCC", "CCG", "CCT"}}, {"Q",{"CAA", "CAG"}}, {"R",{"AGA","AGG", "CGA", "CGC", "CGG", "CGT"}}, 
			{"S",{"TCA", "TCC", "TCG", "TCT"}}, {"T",{"ACA", "ACC", "ACG", "ACT"}}, {"V",{"GTA", "GTC", "GTG", "GTT"}}, {"W",{"TGG"}},
			{"Y",{"TAC", "TAT"}}, {"Z", {"AGC", "AGT"}}, {"X",{"TAA", "TAG", "TGA"}}};

std::map<std::string, std::vector<std::string>> SequenceSummary::AAToCodonsMapWithoutReference = {{"A",{"GCA","GCC","GCG"}}, 
	{"C",{"TGC"}}, {"D",{"GAC"}}, {"E",{"GAA"}}, {"F",{"TTC"}},
	{"G",{"GGA", "GGC", "GGG"}}, {"H",{"CAC"}}, {"I",{"ATA", "ATC"}},
	{"K",{"AAA"}}, {"L",{"CTA", "CTC", "CTG","CTT", "TTA"}}, {"M", {}},
	{"N",{"AAC"}}, {"P",{"CCA", "CCC", "CCG"}}, {"Q",{"CAA"}}, {"R",{"AGA","AGG", "CGA", "CGC", "CGG"}}, 
	{"S",{"TCA", "TCC", "TCG"}}, {"T",{"ACA", "ACC", "ACG"}}, {"V",{"GTA", "GTC", "GTG"}}, {"W", {}},
	{"Y",{"TAC"}}, {"Z", {"AGC"}},{"X",{"TAA", "TAG"}}};

std::map<std::string, std::string> SequenceSummary::codonToAAMap = {{"GCA", "A"}, {"GCC", "A"}, {"GCG", "A"},
	{"GCT", "A"}, {"TGC", "C"}, {"TGT", "C"}, {"GAC", "D"}, {"GAT",  "D"}, {"GAA",  "E"}, {"GAG", "E"}, {"TTC", "F"}, {"TTT", "F"},
	{"GGA", "G"}, {"GGC", "G"}, {"GGG", "G"}, {"GGT", "G"}, {"CAC", "H"}, {"CAT", "H"}, {"ATA", "I"}, {"ATC", "I"}, {"ATT", "I"},
	{"AAA", "K"}, {"AAG", "K"}, {"CTA", "L"}, {"CTC", "L"}, {"CTG", "L"}, {"CTT", "L"}, {"TTA", "L"}, {"TTG", "L"}, {"ATG", "M"},
	{"AAC", "N"}, {"AAT", "N"}, {"CCA", "P"}, {"CCC", "P"}, {"CCG", "P"}, {"CCT", "P"}, {"CAA", "Q"}, {"CAG", "Q"}, {"AGA", "R"},
	{"AGG", "R"}, {"CGA", "R"}, {"CGC", "R"}, {"CGG", "R"}, {"CGT", "R"}, {"TCA", "S"}, {"TCC", "S"}, {"TCG", "S"}, {"TCT", "S"},
	{"ACA", "T"}, {"ACC", "T"}, {"ACG", "T"}, {"ACT", "T"}, {"GTA", "V"}, {"GTC", "V"}, {"GTG", "V"}, {"GTT", "V"}, {"TGG", "W"},
	{"TAC", "Y"}, {"TAT", "Y"}, {"AGC", "Z"}, {"AGT", "Z"}, {"TAA", "X"}, {"TAG", "X"}, {"TGA", "X"}};




// const std::map<std::string, std::vector<std::string>> SequenceSummary::AAToCodonsMapWithReference_12 = {{"A",{"GCA","GCC","GCG","GCT"}}, 
// 			{"C",{"TGC", "TGT"}}, {"D",{"GAC", "GAT"}}, {"E",{"GAA", "GAG"}}, {"F",{"TTC", "TTT"}},
// 			{"G",{"GGA", "GGC", "GGG", "GGT"}}, {"H",{"CAC", "CAT"}}, {"I",{"ATA", "ATC", "ATT"}},
// 			{"K",{"AAA", "AAG"}}, {"L",{"CTA", "CTC", "CTT", "TTA", "TTG"}}, {"M", {"ATG"}},
// 			{"N",{"AAC", "AAT"}}, {"P",{"CCA", "CCC", "CCG", "CCT"}}, {"Q",{"CAA", "CAG"}}, {"R",{"AGA","AGG","CGA", "CGC", "CGG", "CGT"}}, 
// 			{"S",{"TCA", "TCC", "TCG", "TCT"}}, {"T",{"ACA", "ACC", "ACG", "ACT"}}, {"V",{"GTA", "GTC", "GTG", "GTT"}}, {"W",{"TGG"}},
// 			{"Y",{"TAC", "TAT"}}, {"Z", {"AGC", "AGT", "CTG"}}, {"X",{"TAA", "TAG", "TGA"}}};

// const std::map<std::string, std::vector<std::string>> SequenceSummary::AAToCodonsMapWithoutReference_12 = {{"A",{"GCA","GCC","GCG"}}, 
// 	{"C",{"TGC"}}, {"D",{"GAC"}}, {"E",{"GAA"}}, {"F",{"TTC"}},
// 	{"G",{"GGA", "GGC", "GGG"}}, {"H",{"CAC"}}, {"I",{"ATA", "ATC"}},
// 	{"K",{"AAA"}}, {"L",{"CTA", "CTC","CTT", "TTA"}},{"M", {}},
// 	{"N",{"AAC"}}, {"P",{"CCA", "CCC", "CCG"}}, {"Q",{"CAA"}}, {"R",{"AGA","AGG", "CGA", "CGC", "CGG"}}, 
// 	{"S",{"TCA", "TCC", "TCG"}}, {"T",{"ACA", "ACC", "ACG"}}, {"V",{"GTA", "GTC", "GTG"}}, {"W", {}},
// 	{"Y",{"TAC"}}, {"Z", {"AGC", "AGT"}},{"X",{"TAA", "TAG"}}};

// const std::map<std::string, std::string> SequenceSummary::codonToAAMap_12 = {{"GCA", "A"}, {"GCC", "A"}, {"GCG", "A"},
// 	{"GCT", "A"}, {"TGC", "C"}, {"TGT", "C"}, {"GAC", "D"}, {"GAT",  "D"}, {"GAA",  "E"}, {"GAG", "E"}, {"TTC", "F"}, {"TTT", "F"},
// 	{"GGA", "G"}, {"GGC", "G"}, {"GGG", "G"}, {"GGT", "G"}, {"CAC", "H"}, {"CAT", "H"}, {"ATA", "I"}, {"ATC", "I"}, {"ATT", "I"},
// 	{"AAA", "K"}, {"AAG", "K"}, {"CTA", "L"}, {"CTC", "L"}, {"CTG", "Z"}, {"CTT", "L"}, {"TTA", "L"}, {"TTG", "L"}, {"ATG", "M"},
// 	{"AAC", "N"}, {"AAT", "N"}, {"CCA", "P"}, {"CCC", "P"}, {"CCG", "P"}, {"CCT", "P"}, {"CAA", "Q"}, {"CAG", "Q"}, {"AGA", "R"},
// 	{"AGG", "R"}, {"CGA", "R"}, {"CGC", "R"}, {"CGG", "R"}, {"CGT", "R"}, {"TCA", "S"}, {"TCC", "S"}, {"TCG", "S"}, {"TCT", "S"},
// 	{"ACA", "T"}, {"ACC", "T"}, {"ACG", "T"}, {"ACT", "T"}, {"GTA", "V"}, {"GTC", "V"}, {"GTG", "V"}, {"GTT", "V"}, {"TGG", "W"},
// 	{"TAC", "Y"}, {"TAT", "Y"}, {"AGC", "Z"}, {"AGT", "Z"}, {"TAA", "X"}, {"TAG", "X"}, {"TGA", "X"}};


// std::map<std::string, std::vector<std::string>> SequenceSummary::AAToCodonsMapWithReference;
// std::map<std::string, std::vector<std::string>> SequenceSummary::AAToCodonsMapWithoutReference;
// std::map<std::string, std::string> SequenceSummary::codonToAAMap;





//------------------------------------------------//
//---------- Constructors & Destructors ----------//
//------------------------------------------------//


SequenceSummary::SequenceSummary()
{
	clear();
}


SequenceSummary::SequenceSummary(const std::string& sequence)
{
	clear();
	processSequence(sequence);
}


SequenceSummary::SequenceSummary(const SequenceSummary& other)
{
	codonPositions = other.codonPositions;
	ncodons = other.ncodons;
	naa = other.naa;
	RFPCount = other.RFPCount;
	sumRFPCount = other.sumRFPCount;
	positionCodonID = other.positionCodonID;
}


SequenceSummary& SequenceSummary::operator=(const SequenceSummary& rhs)
{
	if (this == &rhs) return *this; // handle self assignment

	codonPositions = rhs.codonPositions;
    ncodons = rhs.ncodons;
    naa = rhs.naa;
	RFPCount = rhs.RFPCount;
	sumRFPCount = rhs.sumRFPCount;
	positionCodonID = rhs.positionCodonID;

	return *this;
}


bool SequenceSummary::operator==(const SequenceSummary& other) const
{
	bool match = true;

	if (this->codonPositions != other.codonPositions) { match = false; }
	if (this->ncodons != other.ncodons) { match = false; }
	if (this->naa != other.naa) { match = false; }
	if (this->RFPCount != other.RFPCount) {match = false; }
	if (this->sumRFPCount != other.sumRFPCount) {match = false; }
	if (this->positionCodonID != other.positionCodonID) { match = false; }

	return match;
}


SequenceSummary::~SequenceSummary()
{
	//dtor
}





//-------------------------------------------------//
//---------- Data Manipulation Functions ----------//
//-------------------------------------------------//

void SequenceSummary::setCodonTable(unsigned codon_table_number)
{

	if (codon_table_number == 1)
	{
		SequenceSummary::AAToCodonsMapWithReference = {{"A",{"GCA","GCC","GCG","GCT"}}, 
			{"C",{"TGC", "TGT"}}, {"D",{"GAC", "GAT"}}, {"E",{"GAA", "GAG"}}, {"F",{"TTC", "TTT"}},
			{"G",{"GGA", "GGC", "GGG", "GGT"}}, {"H",{"CAC", "CAT"}}, {"I",{"ATA", "ATC", "ATT"}},
			{"K",{"AAA", "AAG"}}, {"L",{"CTA", "CTC", "CTG","CTT", "TTA", "TTG"}}, {"M", {"ATG"}},
			{"N",{"AAC", "AAT"}}, {"P",{"CCA", "CCC", "CCG", "CCT"}}, {"Q",{"CAA", "CAG"}}, {"R",{"AGA","AGG", "CGA", "CGC", "CGG", "CGT"}}, 
			{"S",{"TCA", "TCC", "TCG", "TCT"}}, {"T",{"ACA", "ACC", "ACG", "ACT"}}, {"V",{"GTA", "GTC", "GTG", "GTT"}}, {"W",{"TGG"}},
			{"Y",{"TAC", "TAT"}}, {"Z", {"AGC", "AGT"}}, {"X",{"TAA", "TAG", "TGA"}}};

		SequenceSummary::AAToCodonsMapWithoutReference = {{"A",{"GCA","GCC","GCG"}}, 
			{"C",{"TGC"}}, {"D",{"GAC"}}, {"E",{"GAA"}}, {"F",{"TTC"}},
			{"G",{"GGA", "GGC", "GGG"}}, {"H",{"CAC"}}, {"I",{"ATA", "ATC"}},
			{"K",{"AAA"}}, {"L",{"CTA", "CTC", "CTG","CTT", "TTA"}}, {"M", {}},
			{"N",{"AAC"}}, {"P",{"CCA", "CCC", "CCG"}}, {"Q",{"CAA"}}, {"R",{"AGA","AGG", "CGA", "CGC", "CGG"}}, 
			{"S",{"TCA", "TCC", "TCG"}}, {"T",{"ACA", "ACC", "ACG"}}, {"V",{"GTA", "GTC", "GTG"}}, {"W", {}},
			{"Y",{"TAC"}}, {"Z", {"AGC"}},{"X",{"TAA", "TAG"}}};

		SequenceSummary::codonToAAMap = {{"GCA", "A"}, {"GCC", "A"}, {"GCG", "A"},
			{"GCT", "A"}, {"TGC", "C"}, {"TGT", "C"}, {"GAC", "D"}, {"GAT",  "D"}, {"GAA",  "E"}, {"GAG", "E"}, {"TTC", "F"}, {"TTT", "F"},
			{"GGA", "G"}, {"GGC", "G"}, {"GGG", "G"}, {"GGT", "G"}, {"CAC", "H"}, {"CAT", "H"}, {"ATA", "I"}, {"ATC", "I"}, {"ATT", "I"},
			{"AAA", "K"}, {"AAG", "K"}, {"CTA", "L"}, {"CTC", "L"}, {"CTG", "L"}, {"CTT", "L"}, {"TTA", "L"}, {"TTG", "L"}, {"ATG", "M"},
			{"AAC", "N"}, {"AAT", "N"}, {"CCA", "P"}, {"CCC", "P"}, {"CCG", "P"}, {"CCT", "P"}, {"CAA", "Q"}, {"CAG", "Q"}, {"AGA", "R"},
			{"AGG", "R"}, {"CGA", "R"}, {"CGC", "R"}, {"CGG", "R"}, {"CGT", "R"}, {"TCA", "S"}, {"TCC", "S"}, {"TCG", "S"}, {"TCT", "S"},
			{"ACA", "T"}, {"ACC", "T"}, {"ACG", "T"}, {"ACT", "T"}, {"GTA", "V"}, {"GTC", "V"}, {"GTG", "V"}, {"GTT", "V"}, {"TGG", "W"},
			{"TAC", "Y"}, {"TAT", "Y"}, {"AGC", "Z"}, {"AGT", "Z"}, {"TAA", "X"}, {"TAG", "X"}, {"TGA", "X"}};
	}
	if (codon_table_number == 12)
	{
		std::map<std::string, std::vector<std::string>>::iterator it_1;
		std::map<std::string, std::string>::iterator it_2;

		std::string L("L");
		std::string CTG("CTG");

		it_1 = AAToCodonsMapWithReference.find(L);
		it_1 -> second = {"CTA", "CTC", "CTT", "TTA", "TTG"};

		
		it_1 = AAToCodonsMapWithoutReference.find(L);
		it_1 -> second = {"CTA", "CTC", "CTT", "TTA"};

		

		it_2 = codonToAAMap.find(CTG);
		it_2 -> second = "J";
	}
	
}


unsigned SequenceSummary::getAACountForAA(std::string aa)
{
	return naa[aaToIndex.find(aa)->second];
}


unsigned SequenceSummary::getAACountForAA(unsigned aaIndex)
{
	return naa[aaIndex];
}


unsigned SequenceSummary::getCodonCountForCodon(std::string& codon)
{
	return ncodons[codonToIndex(codon)];
}


unsigned SequenceSummary::getCodonCountForCodon(unsigned codonIndex)
{
	return ncodons[codonIndex];
}


std::vector <unsigned> *SequenceSummary::getCodonPositions(std::string codon)
{
	unsigned codonIndex = codonToIndex(codon);
	return getCodonPositions(codonIndex);
}


std::vector <unsigned> *SequenceSummary::getCodonPositions(unsigned index)
{
	return &codonPositions[index];
}





//-----------------------------------//
//---------- RFP Functions ----------//
//-----------------------------------//


/* initRFPCount (NOT EXPOSED)
 * Arguments: A number representing the number of RFP categories.
 * Initializes the vector of vectors that describes, for each category, the RFPCount vector.
 * Note: When programming in C, this function should be called before manipulating RFPCount.
 * Reminder: It must be called after setting the sequence, since that function
 * resets the sequence summary, including the RFPCount.
 */
void SequenceSummary::initRFPCount(unsigned numCategories)
{
	RFPCount.resize(numCategories);
}


/* getRFPCount (NOT EXPOSED)
 * Arguments: A number representing the RFP category to return (default 0)
 * Returns the RFPCount vector for the category index specified.
 * Note: If initRFPCount is not called beforehand, it is called now to return a vector of 0s.
 */
std::vector<unsigned> SequenceSummary::getRFPCount(unsigned RFPCountColumn)
{
	// Note: If the user forgets to initRFPCount manually, this statement is executed but returns an empty vector.
	if (RFPCount.size() < RFPCountColumn + 1) initRFPCount(RFPCountColumn + 1);
	return RFPCount[RFPCountColumn];
}


/* getSingleRFPCount (NOT EXPOSED)
 * Arguments: The position of a single RFP value to return for the given RFP category (default 0)
 * Returns the integer RFPCount value for the category index at the position specified.
 * Note: If initRFPCount is not called beforehand, it is called now to return a value of 0.
 */
unsigned SequenceSummary::getSingleRFPCount(unsigned position, unsigned RFPCountColumn)
{
	if (RFPCount.size() < RFPCountColumn + 1) initRFPCount(RFPCountColumn + 1);
	return RFPCount[RFPCountColumn][position];
}


/* setRFPCount (NOT EXPOSED)
 * Arguments: A vector argument to set the RFP count to for the given RFP category (default 0)
 * Sets the RFPCount vector for the category index specified to the vector argument given.
 * Note: If initRFPCount is not called beforehand, it is called now.
 */
void SequenceSummary::setRFPCount(std::vector <unsigned> arg, unsigned RFPCountColumn)
{
	if (RFPCount.size() < RFPCountColumn + 1) initRFPCount(RFPCountColumn + 1);
	RFPCount[RFPCountColumn] = arg;
}


/* initSumRFPCount (NOT EXPOSED)
 * Arguments: A number representing the number of RFP categories.
 * Initializes the vector of arrays that describes, for each category, the sumRFPCount vector.
 * Note: When programming in C, this function should be called before manipulating sumRFPCount.
 * Reminder: It must be called after setting the sequence, since that function
 * resets the sequence summary, including the sumRFPCount.
 */
void SequenceSummary::initSumRFPCount(unsigned numCategories)
{
	sumRFPCount.resize(numCategories);
	for (unsigned i = 0; i < numCategories; i++)
		sumRFPCount[i].fill(0);
}


/* getSumRFPCount (NOT EXPOSED)
 * Arguments: A number representing the RFP category to return (default 0)
 * Returns the sumRFPCount array of size 64 for the category index specified.
 * Note: If initSumRFPCount is not called beforehand, it is called now to return an array of 0s.
 */
std::array <unsigned, 64> SequenceSummary::getSumRFPCount(unsigned RFPCountColumn)
{
	if (sumRFPCount.size() < RFPCountColumn + 1) initSumRFPCount(RFPCountColumn + 1);
	return sumRFPCount[RFPCountColumn];
}


/* setSumRFPCount (NOT EXPOSED)
 * Arguments: an array argument to set the sumRFPCount (aka RFPValue) to for the given RFP category (default 0)
 * Sets the sumRFPCount vector for the category index specified to the array argument given.
 * Note: If initSumRFPCount is not called beforehand, it is called now.
 */
void SequenceSummary::setSumRFPCount(std::array <unsigned, 64> arg, unsigned RFPCountColumn)
{
	if (sumRFPCount.size() < RFPCountColumn + 1) initSumRFPCount(RFPCountColumn + 1);
	sumRFPCount[RFPCountColumn] = arg;
}


unsigned SequenceSummary::getSumTotalRFPCount(unsigned RFPCountColumn)
{
    if (sumRFPCount.size() < RFPCountColumn + 1) initSumRFPCount(RFPCountColumn + 1);
    unsigned sum = 0;
    for (unsigned i = 0u; i < sumRFPCount[RFPCountColumn].size(); i++)
    {
        sum += sumRFPCount[RFPCountColumn][i];
    }

    return sum;
}

/* getCodonSpecificSumRFPCount (by codon string) (RCPP EXPOSED VIA WRAPPER)
 * Arguments: A three-character codon string to get the RFP value of, a number representing the RFP category to return (default 0)
 * Returns the RFP value of the codon string for the category index specified.
 * Note: If initSumRFPCount is not called beforehand, it is called now to return a value of 0.
 * Wrapped by Gene::getSumRFPCountForCodon on the R-side.
 */
unsigned SequenceSummary::getCodonSpecificSumRFPCount(std::string codon, unsigned RFPCountColumn)
{
	if (sumRFPCount.size() < RFPCountColumn + 1){
	    initSumRFPCount(RFPCountColumn + 1);
	}
	return sumRFPCount[RFPCountColumn][codonToIndex(codon)];
}


/* getCodonSpecificSumRFPCount (by codon index) (NOT EXPOSED)
 * Arguments: A codon index to get the RFP value of, a number representing the RFP category to return (default 0)
 * Returns the RFP value at the codon index for the category index specified.
 * Note: If initSumRFPCount is not called beforehand, it is called now to return a value of 0.
 */
unsigned SequenceSummary::getCodonSpecificSumRFPCount(unsigned codonIndex, unsigned RFPCountColumn)
{
	if (sumRFPCount.size() < RFPCountColumn + 1){
		initSumRFPCount(RFPCountColumn + 1);
	}
    return sumRFPCount[RFPCountColumn][codonIndex];
}


/* setCodonSpecificSumRFPCount (NOT EXPOSED)
 * Arguments: A codon index, the value to set the RFP value to, and a number representing the RFP category (default 0)
 * Sets the RFP value at the codon index for the category index specified.
 * Note: If initSumRFPCount is not called beforehand, it is called now to return initialize the vector of vectors.
 */
void SequenceSummary::setCodonSpecificSumRFPCount(unsigned codonIndex, unsigned value, unsigned RFPCountColumn)
{
    if (sumRFPCount.size() < RFPCountColumn + 1) initSumRFPCount(RFPCountColumn + 1);
    sumRFPCount[RFPCountColumn][codonIndex] = value;
}


/* getPositionCodonID (NOT EXPOSED)
 * Arguments: None.
 * Returns the vector of codon IDs for each position.
 */
std::vector <unsigned> SequenceSummary::getPositionCodonID()
{
	return positionCodonID;
}


/* setPositionCodonID (NOT EXPOSED)
 * Arguments: An vector to be set as the vector of codonIDs for each position.
 * Sets the positionCodonID vector specified to the vector argument given.
 */
void SequenceSummary::setPositionCodonID(std::vector <unsigned> arg)
{
    positionCodonID = arg;
}


//------------------------------------//
//---------- Other Functions ---------//
//------------------------------------//


void SequenceSummary::clear()
{
	codonPositions.clear();
	RFPCount.clear();
	sumRFPCount.clear();
	ncodons.fill(0);
	naa.fill(0);
}


// Returns a bool for error checking purposes related to setSequence in Gene.cpp
bool SequenceSummary::processSequence(const std::string& sequence)
{
	bool check = true;
	codonPositions.clear();
	codonPositions.resize(64);
	ncodons.fill(0);
	naa.fill(0);
	for (unsigned i = 0u; i < sequence.length(); i += 3)
	{
		std::string codon = sequence.substr(i, 3);
		codon[0] = (char)std::toupper(codon[0]);
		codon[1] = (char)std::toupper(codon[1]);
		codon[2] = (char)std::toupper(codon[2]);

		unsigned codonID = codonToIndex(codon);
		if (codonID != 64) // if codon id == 64 => codon not found. Ignore, probably N
		{
			int aaID = codonToAAIndex(codon);
			ncodons[codonID]++;
			naa[aaID]++;
			codonPositions[codonID].push_back(i / 3);
			
		}
		else
		{
			my_printError("WARNING: Codon % not recognized!\n Codon will be ignored!\n", codon);
			check = false;
		}
	}
	return check;
}


bool SequenceSummary::processPA(std::vector<std::vector<int>> table)
{
    // Table format: Each line of input from a .csv (.pa) file, ordered:
    // unknown size table (nRows, aka table.size()), each row a vector:
    // position, codon, category1, ... (may be more than one category)

	bool check = true;
	codonPositions.resize(64);
	unsigned nRows = (unsigned)table.size();
	positionCodonID.resize(nRows);

	// There should be at least 1 table entry to get to this point, so this should be a valid operation
    unsigned numCats = (unsigned)table[0].size() - 2; // numCats = after position, codon.
	initRFPCount(numCats);
	sumRFPCount.resize(numCats);

	for (unsigned j = 0; j < numCats; j++)
	{
		RFPCount[j].resize(nRows);
		sumRFPCount[j].fill(0);
	}

	for (unsigned i = 0; i < nRows; i++)
	{
		std::vector <int> row = table[i];

		unsigned codonID = (unsigned)row[1];
		std::string codon = indexToCodon(codonID);
		// Note: Don't bother writing a function to convert codonIndex to aaIndex
        // Would just perform the exact same steps anyway; redundant code.
		if (codonID != 64) // if codon id == 64 => codon not found. Ignore, probably N
		{
			int aaID = codonToAAIndex(codon);
			ncodons[codonID]++;
			naa[aaID]++;
			codonPositions[codonID].push_back((unsigned) row[0]);
			positionCodonID[row[0]] = codonID;

			for (unsigned j = 0; j < numCats; j++)
			{
				// Category j has an RFPCount at the position equal to the 2-indexed (after position, codon) value of j.
				RFPCount[j][row[0]] = row[j + 2];
				if (row[j+2] > 0) sumRFPCount[j][codonID] += row[j + 2];
                // Recall: We store RFP counts < 0, but do not need to process this information in calculations
                // So we only add to the sumRFPCount if the value is "valid" (> 0).
			}
		}
		else
		{
			my_printError("WARNING: Codon % not recognized!\n Codon will be ignored!\n", codon);
			check = false;
		}
	}

	return check;
}

//TODO: Turn into equivalent PANSE
bool SequenceSummary::processPANSE(std::vector<std::vector<int>> table)
{
    // Table format: Each line of input from a .csv (.panse) file, ordered:
    // unknown size table (nRows, aka table.size()), each row a vector:
    // position, codon, category1, ... (may be more than one category)

	bool check = true;
	codonPositions.resize(64);
	unsigned nRows = (unsigned)table.size();
	positionCodonID.resize(nRows);

	// There should be at least 1 table entry to get to this point, so this should be a valid operation
    unsigned numCats = (unsigned)table[0].size() - 2; // numCats = after position, codon.
	initRFPCount(numCats);
	sumRFPCount.resize(numCats);

	for (unsigned j = 0; j < numCats; j++)
	{
		RFPCount[j].resize(nRows);
		sumRFPCount[j].fill(0);
	}

	for (unsigned i = 0; i < nRows; i++)
	{
		std::vector <int> row = table[i];

		unsigned codonID = (unsigned)row[1];
		std::string codon = indexToCodon(codonID);
		// Note: Don't bother writing a function to convert codonIndex to aaIndex
        // Would just perform the exact same steps anyway; redundant code.
		if (codonID != 64) // if codon id == 64 => codon not found. Ignore, probably N
		{
			int aaID = codonToAAIndex(codon);
			ncodons[codonID]++;
			naa[aaID]++;
			codonPositions[codonID].push_back((unsigned) row[0]);
			positionCodonID[row[0]] = codonID;

			for (unsigned j = 0; j < numCats; j++)
			{
				// Category j has an RFPCount at the position equal to the 2-indexed (after position, codon) value of j.
				RFPCount[j][row[0]] = row[j + 2];
				if (row[j+2] > 0) sumRFPCount[j][codonID] += row[j + 2];
                // Recall: We store RFP counts < 0, but do not need to process this information in calculations
                // So we only add to the sumRFPCount if the value is "valid" (> 0).
			}
		}
		else
		{
			my_printError("WARNING: Codon % not recognized!\n Codon will be ignored!\n", codon);
			check = false;
		}
	}
	return check;
}





//--------------------------------------//
//---------- Static Functions ----------//
//--------------------------------------//


unsigned SequenceSummary::AAToAAIndex(std::string aa)
{
	return SequenceSummary::aaToIndex.find(aa) -> second;
}


//TODO: test this function. See note in testSequenceSummary.R.
// Note: From function definition in header, default forParamVector is false.
void SequenceSummary::AAIndexToCodonRange(unsigned aaIndex, unsigned& startAAIndex, unsigned& endAAIndex, bool forParamVector)
{
	std::string aa = indexToAA(aaIndex);
	AAToCodonRange(aa, startAAIndex, endAAIndex, forParamVector);
}

<<<<<<< HEAD

=======
>>>>>>> 7e358f8caac68038c33c89bba3922e7d9f4dd6a3
//std::array<unsigned, 2>
// Note: From function definition in header, default forParamVector is false.
// Returns the range of index values in the CodonTable for codons corresponding to a given amino acid. The function is overloaded and uses a wrapper function to map from an amino acid index value rather than a string. (which I believe is only a single char).
//
// param aa corresponds to an entry in groupList which are separately defined for each parametrer type, e.g. in ROCParameter.cpp
// param v2 Second value
// return Product of v1 and v2
// [[Rcpp__export]] Change __ to :: to export function description. Seems to be breaking things as of now. 
void SequenceSummary::AAToCodonRange(std::string aa, unsigned& startAAIndex, unsigned& endAAIndex, bool forParamVector)
{
	//aa = (char)std::toupper(aa[0]); CEDRIC: commented out for performance. Put back in if necessary!
	// switch statement is a lot faster than a chain of if else!
	//unsigned startAAIndex = 0u;
	//unsigned endAAIndex = 0u;
	char AA = aa[0];

	switch (AA)
	{
	case 'A':
		if (!forParamVector) { startAAIndex = 0; endAAIndex = 4; }
		else { startAAIndex = 0; endAAIndex = 3; }
		break;
	case 'C':
		if (!forParamVector) { startAAIndex = 4; endAAIndex = 6; }
		else { startAAIndex = 3; endAAIndex = 4; }
		break;
	case 'D':
		if (!forParamVector) { startAAIndex = 6; endAAIndex = 8; }
		else { startAAIndex = 4; endAAIndex = 5; }
		break;
	case 'E':
		if (!forParamVector) { startAAIndex = 8; endAAIndex = 10; }
		else { startAAIndex = 5; endAAIndex = 6; }
		break;
	case 'F':
		if (!forParamVector) { startAAIndex = 10; endAAIndex = 12; }
		else { startAAIndex = 6; endAAIndex = 7; }
		break;
	case 'G':
		if (!forParamVector) { startAAIndex = 12; endAAIndex = 16; }
		else { startAAIndex = 7; endAAIndex = 10; }
		break;
	case 'H':
		if (!forParamVector) { startAAIndex = 16; endAAIndex = 18; }
		else { startAAIndex = 10; endAAIndex = 11; }
		break;
	case 'I':
		if (!forParamVector) { startAAIndex = 18; endAAIndex = 21; }
		else { startAAIndex = 11; endAAIndex = 13; }
		break;
	case 'K':
		if (!forParamVector) { startAAIndex = 21; endAAIndex = 23; }
		else { startAAIndex = 13; endAAIndex = 14; }
		break;
	case 'L':
		if (!forParamVector) { startAAIndex = 23; endAAIndex = 29; }
		else { startAAIndex = 14; endAAIndex = 19; }
		break;
	case 'M':
		if (!forParamVector) { startAAIndex = 29; endAAIndex = 30; }
		else { startAAIndex = 19; endAAIndex = 19; }
		break;
	case 'N':
		if (!forParamVector) { startAAIndex = 30; endAAIndex = 32; }
		else { startAAIndex = 19; endAAIndex = 20; }
		break;
	case 'P':
		if (!forParamVector) { startAAIndex = 32; endAAIndex = 36; }
		else { startAAIndex = 20; endAAIndex = 23; }
		break;
	case 'Q':
		if (!forParamVector) { startAAIndex = 36; endAAIndex = 38; }
		else { startAAIndex = 23; endAAIndex = 24; }
		break;
	case 'R':
		if (!forParamVector) { startAAIndex = 38; endAAIndex = 44; }
		else { startAAIndex = 24; endAAIndex = 29; }
		break;
	case 'S':
		if (!forParamVector) { startAAIndex = 44; endAAIndex = 48; }
		else { startAAIndex = 29; endAAIndex = 32; }
		break;
	case 'T':
		if (!forParamVector) { startAAIndex = 48; endAAIndex = 52; }
		else { startAAIndex = 32; endAAIndex = 35; }
		break;
	case 'V':
		if (!forParamVector) { startAAIndex = 52; endAAIndex = 56; }
		else { startAAIndex = 35; endAAIndex = 38; }
		break;
	case 'W':
		if (!forParamVector) { startAAIndex = 56; endAAIndex = 57; }
		else { startAAIndex = 38; endAAIndex = 38; }
		break;
	case 'Y':
		if (!forParamVector) { startAAIndex = 57; endAAIndex = 59; }
		else { startAAIndex = 38; endAAIndex = 39; }
		break;
	case 'Z':
		if (!forParamVector) { startAAIndex = 59; endAAIndex = 61; }
		else { startAAIndex = 39; endAAIndex = 40; }
		break;
	case 'X':
		if (!forParamVector) { startAAIndex = 61; endAAIndex = 64; }
		else { startAAIndex = 40; endAAIndex = 40; }
		break;
	default: // INVALID AA
		startAAIndex = 0;
		endAAIndex = 0;
		my_print("%\n", AA);
		my_printError("Invalid AA given, returning 0,0\n");
		break;
	}
}


// Note: From function definition in header, default category is 0.
std::vector<std::string> SequenceSummary::AAToCodon(std::string aa, bool forParamVector)
{
	std::vector<std::string> RV;
	aa = (char) std::toupper(aa[0]);

	if (forParamVector)
	{
		RV = SequenceSummary::AAToCodonsMapWithoutReference.find(aa) -> second;
	}
	else
	{
		RV = SequenceSummary::AAToCodonsMapWithReference.find(aa) -> second;
	}
	return RV;
}

// Note: From function definition in header, default category is 0.
std::vector<unsigned> SequenceSummary::AAToCodonIndex(std::string aa, bool forParamVector)
{
	std::vector<std::string> codon_list;
	std::vector<unsigned> RV;
	aa = (char) std::toupper(aa[0]);

	if (forParamVector)
	{
		codon_list = SequenceSummary::AAToCodonsMapWithoutReference.find(aa) -> second;
		for (unsigned i = 0; i < codon_list.size();i++)
		{
			RV.push_back(SequenceSummary::codonToIndexWithoutReference.find(codon_list[i])->second);
		}
	}
	else
	{
		codon_list = SequenceSummary::AAToCodonsMapWithReference.find(aa) -> second;
		for (unsigned i = 0; i < codon_list.size();i++)
		{
			RV.push_back(SequenceSummary::codonToIndexWithReference.find(codon_list[i])->second);
		}
	}
	return RV;
}


std::string SequenceSummary::codonToAA(std::string& codon)
{
	codon[0] = (char) std::toupper(codon[0]);
	codon[1] = (char) std::toupper(codon[1]);
	codon[2] = (char) std::toupper(codon[2]);

	//std::transform(codon.begin(), codon.end(), codon.begin(), ::toupper);
	std::map<std::string, std::string>::iterator it;
	std::string aa = "#";
	it = SequenceSummary::codonToAAMap.find(codon);
	if (it != SequenceSummary::codonToAAMap.end())
	{
		aa = it -> second;
	}
	else
	{
		my_print("WARNING: Codon % has no corresponding amino acid. Returning A by default",codon);
		aa = "A";
	}
	return aa;
}

// Note: From function definition in header, default forParamVector is false.
unsigned SequenceSummary::codonToIndex(std::string& codon, bool forParamVector)
{
	unsigned i = 0;
	codon[0] = (char) std::toupper(codon[0]);
	codon[1] = (char) std::toupper(codon[1]);
	codon[2] = (char) std::toupper(codon[2]);
	if (((codon[0] != 'A') && (codon[0] != 'C') && (codon[0] != 'G') && (codon[0] != 'T')) ||
		((codon[1] != 'A') && (codon[1] != 'C') && (codon[1] != 'G') && (codon[1] != 'T')) ||
		((codon[2] != 'A') && (codon[2] != 'C') && (codon[2] != 'G') && (codon[2] != 'T')))
	{
		i = 64;
	}
	else
	{
		if (forParamVector)
			i = SequenceSummary::codonToIndexWithoutReference.find(codon) -> second;
		else
			i = SequenceSummary::codonToIndexWithReference.find(codon) -> second;
	}
	return i;
}


unsigned SequenceSummary::codonToAAIndex(std::string& codon)
{
	std::string aa = codonToAA(codon);
	return aaToIndex.find(aa) -> second;
}


std::string SequenceSummary::indexToAA(unsigned aaIndex)
{
	return AminoAcidArray[aaIndex];
}


// Note: From function definition in header, default forParamVector is false.
std::string SequenceSummary::indexToCodon(unsigned index, bool forParamVector)
{
	return forParamVector ? codonArrayParameter[index] : codonArray[index];
}


// Note: From function definition in header, default forParamVector is false.
unsigned SequenceSummary::GetNumCodonsForAA(std::string& aa, bool forParamVector)
{
	std::vector<std::string> RV;
	aa = (char) std::toupper(aa[0]);

	if (forParamVector)
	{
		RV = SequenceSummary::AAToCodonsMapWithoutReference.find(aa) -> second;
	
	}
	else
	{
		RV = SequenceSummary::AAToCodonsMapWithReference.find(aa) -> second;
	}
	return (RV.size());

}


char SequenceSummary::complimentNucleotide(char ch)
{
	if ( ch == 'A' ) return 'T';
	else if ( ch == 'T' ) return 'A';
	else if ( ch == 'C' ) return 'G';
	else return 'C';
}


std::vector<std::string> SequenceSummary::aminoAcids()
{
	return AminoAcidArray;
}


std::vector<std::string> SequenceSummary::codons()
{
	std::vector<std::string> RV;
	for (unsigned i = 0; i < 64; i++) RV.push_back(codonArray[i]);
	return RV;
}





// -----------------------------------------------------------------------------------------------------//
// ---------------------------------------- R SECTION --------------------------------------------------//
// -----------------------------------------------------------------------------------------------------//



#ifndef STANDALONE


//---------------------------------//
//---------- RCPP Module ----------//
//---------------------------------//


RCPP_MODULE(SequenceSummary_mod)
{
	class_<SequenceSummary>( "SequenceSummary" );

		//Static Functions:
		Rcpp::function("AAToCodon", &SequenceSummary::AAToCodon, List::create(_["aa"], _["focal"] = false),
				"returns a vector of codons (either containing the reference codon or excluding it) for a given amino acid"); //Used, but will move into Codon Table

		Rcpp::function("codonToAA", &SequenceSummary::codonToAA, List::create(_["codon"]),
		"returns an amino acid string for a given codon string");
		// Note: Unlike AAToCodon, this function works with individual components rather than an entire vector

		Rcpp::function("aminoAcids", &SequenceSummary::aminoAcids, "returns all Amino Acids as one letter code");
		Rcpp::function("codons", &SequenceSummary::codons, "returns all codons or all reference codons");

		Rcpp::function("setCodonTable", &SequenceSummary::setCodonTable, List::create(_["codon_table_number"]),
				"Sets the codon table."); //Used, but will move into Codon Table


}
#endif
