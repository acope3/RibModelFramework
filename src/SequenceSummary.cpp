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
	{"TAC", 38}, {"AGC", 39}};



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
	codonPositions.resize(other.codonPositions.size());
	for (unsigned i = 0u; i < codonPositions.size(); i++) {
		codonPositions[i] = other.codonPositions[i];
	}

	for (unsigned i = 0u; i < 64; i++) {
		ncodons[i] = other.ncodons[i];
	}

	for (unsigned i = 0u; i < 22; i++) {
		naa[i] = other.naa[i];
	}

	for (unsigned i = 0u; i < 64; i++) {
		RFPObserved[i] = other.RFPObserved[i];
	}
}


SequenceSummary& SequenceSummary::operator=(const SequenceSummary& rhs)
{
	if (this == &rhs) return *this; // handle self assignment

	// TODO(CEDRIC): shouldn't a simple = do the job? see http://www.cplusplus.com/reference/vector/vector/operator=/
	codonPositions.resize(rhs.codonPositions.size());
	for (unsigned i = 0u; i < codonPositions.size(); i++) {
		codonPositions[i] = rhs.codonPositions[i];
	}

	for (unsigned i = 0u; i < 64; i++) {
		ncodons[i] = rhs.ncodons[i];
		RFPObserved[i] = rhs.RFPObserved[i];
	}

	for (unsigned i = 0u; i < 22; i++) {
		naa[i] = rhs.naa[i];
	}

	RFP_count = rhs.RFP_count;

	return *this;
}


bool SequenceSummary::operator==(const SequenceSummary& other) const
{
	bool match = true;

	if (this->naa != other.naa) { match = false; }
	if (this->ncodons != other.ncodons) { match = false; }
	if (this->codonPositions != other.codonPositions) { match = false; }
	if (this->RFPObserved != other.RFPObserved) { match = false; }
	if (this->RFP_count != other.RFP_count) {match = false; }

	return match;
}


SequenceSummary::~SequenceSummary()
{
	//dtor
}





//-------------------------------------------------//
//---------- Data Manipulation Functions ----------//
//-------------------------------------------------//


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


unsigned SequenceSummary::getRFPObserved(std::string codon)
{
	return RFPObserved[codonToIndex(codon)];
}


unsigned SequenceSummary::getRFPObserved(unsigned codonIndex)
{
	return RFPObserved[codonIndex];
}


void SequenceSummary::setRFPObserved(unsigned codonIndex, unsigned value)
{
	RFPObserved[codonIndex] = value;
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

std::vector <unsigned> SequenceSummary::getRFP_count()
{
	return RFP_count;
}

void SequenceSummary::setRFP_count(std::vector <unsigned> arg)
{
	RFP_count = arg;
}


//------------------------------------//
//---------- Other Functions ---------//
//------------------------------------//


void SequenceSummary::clear()
{
	codonPositions.clear();
	RFP_count.clear();
	for(unsigned k = 0; k < 64; k++)
	{
		ncodons[k] = 0;
		RFPObserved[k] = 0;
	}
	for(unsigned k = 0; k < 22; k++) { naa[k] = 0; }
}

bool SequenceSummary::processSequence(const std::string& sequence)
{
	//NOTE! Clear() cannot be called in this function because of the RFP model.
	//RFP sets RFPObserved by codon, and not by setting the sequence. This causes
	//the values to be zero during the MCMC.

	bool check = true;
	int codonID;
	int aaID;
	std::string codon;

	codonPositions.resize(64);

	for (unsigned i = 0u; i < sequence.length(); i += 3)
	{
		codon = sequence.substr(i, 3);
		codon[0] = (char)std::toupper(codon[0]);
		codon[1] = (char)std::toupper(codon[1]);
		codon[2] = (char)std::toupper(codon[2]);

		codonID = codonToIndex(codon);
		if (codonID != 64) // if codon id == 64 => codon not found. Ignore, probably N 
		{
			aaID = codonToAAIndex(codon);
			ncodons[codonID]++;
			naa[aaID]++;
			codonPositions[codonID].push_back(i / 3);
		}
		else
		{
#ifndef STANDALONE
			Rf_warning("Codon %s not recognized!\n Codon will be ignored!\n", codon.c_str());
#else
			std::cerr << "WARNING: Codon " << codon << " not recognized!\n Codon will be ignored!\n";
#endif
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


void SequenceSummary::AAIndexToCodonRange(unsigned aaIndex, unsigned& startAAIndex, unsigned& endAAIndex, bool forParamVector)
{
	std::string aa = indexToAA(aaIndex);
	AAToCodonRange(aa, startAAIndex, endAAIndex, forParamVector);
}

//std::array<unsigned, 2>
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
		std::cout << AA << std::endl;
//#ifndef STANDALONE
//		Rf_warning("Invalid Amino Acid given (%s), returning 0,0\n", aa.c_str());
//#else
		std::cerr << "Invalid AA given, returning 0,0\n";
//#endif
		break;
	}
	//std::array<unsigned, 2> aaRange;
	//aaRange[0] = startAAIndex;
	//aaRange[1] = endAAIndex;

	//return aaRange;
}


std::vector<std::string> SequenceSummary::AAToCodon(std::string aa, bool forParamVector)
{
	std::vector <std::string> RV;
	aa = (char) std::toupper(aa[0]);

	unsigned aaStart;
	unsigned aaEnd;
	SequenceSummary::AAToCodonRange(aa, aaStart, aaEnd, forParamVector);
	if(forParamVector){
		for (unsigned i = aaStart; i < aaEnd; i++)
		{
			RV.push_back(codonArrayParameter[i]);
		}
	}else{
		for (unsigned i = aaStart; i < aaEnd; i++)
		{
			RV.push_back(codonArray[i]);
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
	std::string aa = "#";
	//Phenylalanine
	if(!codon.compare("TTT") || !codon.compare("UUU") || !codon.compare("TTC") || !codon.compare("UUC")) aa = "F";
		//Leucine
	else if(!codon.compare("TTA") || !codon.compare("UUA") || !codon.compare("TTG") || !codon.compare("UUG") ||
			!codon.compare("CTT") || !codon.compare("CUU") || !codon.compare("CTC") || !codon.compare("CUC") ||
			!codon.compare("CTA") || !codon.compare("CUA") || !codon.compare("CTG") || !codon.compare("CUG")) aa = "L";
		//Isoleucine
	else if(!codon.compare("ATT") || !codon.compare("AUU") || !codon.compare("ATC") || !codon.compare("AUC") ||
			!codon.compare("ATA") || !codon.compare("AUA")) aa = "I";
		//Methionine
	else if(!codon.compare("ATG") || !codon.compare("AUG")) aa = "M";
		//Valine
	else if(!codon.compare("GTT") || !codon.compare("GUU") || !codon.compare("GTC") || !codon.compare("GUC") ||
			!codon.compare("GTA") || !codon.compare("GUA") || !codon.compare("GTG") || !codon.compare("GUG")) aa = "V";
		//Serine4
	else if(!codon.compare("TCT") || !codon.compare("UCU") || !codon.compare("TCC") || !codon.compare("UCC") ||
			!codon.compare("TCA") || !codon.compare("UCA") || !codon.compare("TCG") || !codon.compare("UCG")) aa = "S";
		//Proline
	else if(!codon.compare("CCT") || !codon.compare("CCU") || !codon.compare("CCC") ||
			!codon.compare("CCA") || !codon.compare("CCG")) aa = "P";
		//Threonine
	else if(!codon.compare("ACT") || !codon.compare("ACU") || !codon.compare("ACC") ||
			!codon.compare("ACA") || !codon.compare("ACG")) aa = "T";
		//Alanine
	else if(!codon.compare("GCT") || !codon.compare("GCU") || !codon.compare("GCC") ||
			!codon.compare("GCA") || !codon.compare("GCG")) aa = "A";
		//Tyrosine
	else if(!codon.compare("TAT") || !codon.compare("UAU") || !codon.compare("TAC") || !codon.compare("UAC")) aa = "Y";
		//Histidine
	else if(!codon.compare("CAT") || !codon.compare("CAU") || !codon.compare("CAC")) aa = "H";
		//Glutamine
	else if(!codon.compare("CAA") || !codon.compare("CAG")) aa = "Q";
		//Asparagine
	else if(!codon.compare("AAT") || !codon.compare("AAU") || !codon.compare("AAC")) aa = "N";
		//Lysine
	else if(!codon.compare("AAA") || !codon.compare("AAG")) aa = "K";
		//Aspartic Acid
	else if(!codon.compare("GAT") || !codon.compare("GAU") || !codon.compare("GAC")) aa = "D";
		//Glutamic Acid
	else if(!codon.compare("GAA") || !codon.compare("GAG")) aa = "E";
		//Cysteine
	else if(!codon.compare("TGT") || !codon.compare("TAT") || !codon.compare("TGC") || !codon.compare("UGC")) aa = "C";
		//Tryptophan
	else if(!codon.compare("TGG") || !codon.compare("UGG")) aa = "W";
		//Arginine
	else if(!codon.compare("CGT") || !codon.compare("CGU") || !codon.compare("CGC") || !codon.compare("CGA") ||
			!codon.compare("CGG") || !codon.compare("AGA") || !codon.compare("AGG")) aa = "R";
		//Serine2
	else if(!codon.compare("AGT") || !codon.compare("AGU") || !codon.compare("AGC")) aa = SequenceSummary::Ser2;
		//Glycine
	else if(!codon.compare("GGT") || !codon.compare("GGU") || !codon.compare("GGC")  || !codon.compare("GGC") ||
			!codon.compare("GGA")  || !codon.compare("GGG")) aa = "G";
		//Stop
	else if (!codon.compare("TAA") || !codon.compare("UAA") || !codon.compare("TAG") || !codon.compare("UAG") ||
			 !codon.compare("TGA") || !codon.compare("UGA")) aa = "X";

	return aa;
}


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
		if(forParamVector)
		{
			i = SequenceSummary::codonToIndexWithoutReference.find(codon) -> second;
		}else{
			i = SequenceSummary::codonToIndexWithReference.find(codon) -> second;
		}
	}
	return i;
}


unsigned SequenceSummary::codonToAAIndex(std::string& codon)
{
	std::string aa = codonToAA(codon);
	return SequenceSummary::aaToIndex.find(aa) -> second;
}


std::string SequenceSummary::indexToAA(unsigned aaIndex)
{
	return AminoAcidArray[aaIndex];
}


std::string SequenceSummary::indexToCodon(unsigned index, bool forParamVector)
{
	return forParamVector ? codonArrayParameter[index] : codonArray[index];
}


unsigned SequenceSummary::GetNumCodonsForAA(std::string& aa, bool forParamVector)
{
	unsigned ncodon = 0;
	char AA = aa[0];
	switch (AA)
	{
	case 'A':
		ncodon = 4;
		break;
	case 'C':
		ncodon = 2;
		break;
	case 'D':
		ncodon = 2;
		break;
	case 'E':
		ncodon = 2;
		break;
	case 'F':
		ncodon = 2;
		break;
	case 'G':
		ncodon = 4;
		break;
	case 'H':
		ncodon = 2;
		break;
	case 'I':
		ncodon = 3;
		break;
	case 'K':
		ncodon = 2;
		break;
	case 'L':
		ncodon = 6;
		break;
	case 'M':
		ncodon = 1;
		break;
	case 'N':
		ncodon = 2;
		break;
	case 'P':
		ncodon = 4;
		break;
	case 'Q':
		ncodon = 2;
		break;
	case 'R':
		ncodon = 6;
		break;
	case 'S':
		ncodon = 4;
		break;
	case 'T':
		ncodon = 4;
		break;
	case 'V':
		ncodon = 4;
		break;
	case 'W':
		ncodon = 1;
		break;
	case 'Y':
		ncodon = 2;
		break;
	case 'Z':
		ncodon = 2;
		break;
	case 'X':
		ncodon = 3;
		break;
	default: // INVALID AA

#ifndef STANDALONE
		Rf_warning("Invalid Amino Acid given (%s), returning 0,0\n", aa.c_str());
#else
		std::cerr << "Invalid aa given, returning 0\n";
#endif
		break;
	}
	return (forParamVector ? (ncodon - 1) : ncodon);
}


char SequenceSummary::complimentNucleotide(char ch)
{
	if( ch == 'A' ) return 'T';
	else if( ch == 'T' ) return 'A';
	else if( ch == 'C' ) return 'G';
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
		Rcpp::function("AAToCodon", &SequenceSummary::AAToCodon, List::create(_["aa"], _["forParamVector"] = false),
				"returns a vector of codons for a given amino acid"); //Used, but will move into Codon Table
		Rcpp::function("aminoAcids", &SequenceSummary::aminoAcids, "returns all Amino Acids as one letter code");
		Rcpp::function("codons", &SequenceSummary::codons, "returns all codons or all reference codons");

}
#endif

