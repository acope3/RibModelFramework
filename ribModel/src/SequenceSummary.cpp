#include "include/SequenceSummary.h"


#include <iostream>


SequenceSummary::SequenceSummary()
{
	clear();
	//ctor
}


SequenceSummary::SequenceSummary(const std::string& sequence)
{
	clear();
	processSequence(sequence);
}


SequenceSummary::~SequenceSummary()
{
	//dtor
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

	/*std::copy(std::begin(other.naa), std::end(other.naa), std::begin(naa));
	std::copy(std::begin(other.RFPObserved), std::end(other.RFPObserved), std::begin(RFPObserved));
	*/
}


SequenceSummary& SequenceSummary::operator=(const SequenceSummary& rhs)
{
	if (this == &rhs) return *this; // handle self assignment
	codonPositions.resize(rhs.codonPositions.size());
	for (unsigned i = 0u; i < codonPositions.size(); i++) {
		codonPositions[i] = rhs.codonPositions[i];
	}

	for (unsigned i = 0u; i < 64; i++) {
		ncodons[i] = rhs.ncodons[i];
	}

	for (unsigned i = 0u; i < 22; i++) {
		naa[i] = rhs.naa[i];
	}

	for (unsigned i = 0u; i < 64; i++) {
		RFPObserved[i] = rhs.RFPObserved[i];
	}

	/*std::copy(std::begin(rhs.ncodons), std::end(rhs.ncodons), std::begin(ncodons));
	std::copy(std::begin(rhs.naa), std::end(rhs.naa), std::begin(naa));
	std::copy(std::begin(rhs.RFPObserved), std::end(rhs.RFPObserved), std::begin(RFPObserved));
	*/
	//assignment operator
	return *this;
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
	return ncodons[SequenceSummary::CodonToIndex(codon)];
}


unsigned SequenceSummary::getCodonCountForCodon(unsigned codonIndex)
{
	return ncodons[codonIndex];
}


unsigned SequenceSummary::getRFPObserved(unsigned codonIndex)
{
	return RFPObserved[codonIndex];
}


void SequenceSummary::setRFPObserved(unsigned codonIndex, unsigned value)
{
	RFPObserved[codonIndex] = value;
}


std::vector <unsigned> SequenceSummary::getCodonPositions(unsigned index)
{
	std::vector <unsigned> rv;

	for (unsigned i = 0; i < codonPositions[index].size(); i++) {
		rv.push_back(codonPositions[index][i]);
	}

	return rv;
}


void SequenceSummary::clear()
{
	codonPositions.clear();
	for(unsigned k = 0; k < 64; k++)
	{
		ncodons[k] = 0;
		RFPObserved[k] = 0;
	}
	for(unsigned k = 0; k < 22; k++) { naa[k] = 0; }
}


bool SequenceSummary::processSequence(const std::string& sequence)
{
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

		codonID = SequenceSummary::CodonToIndex(codon);
		if (codonID != 64) // if codon id == 64 => codon not found. Ignore, probably N 
		{
			aaID = SequenceSummary::CodonToAAIndex(codon);
			ncodons[codonID]++;
			naa[aaID]++;
			codonPositions[codonID].push_back(i / 3);
		}
		else
		{
			std::cerr << "WARNING: Codon " << codon << " not recognized!\n Codon will be ignored!\n";
			check = false;
		}
	}
	return check;
}





/*
 * STATIC FUNCTIONS
 */
const std::string SequenceSummary::Ser2 = "Z";
const std::vector<std::string> SequenceSummary::AminoAcidArray = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", SequenceSummary::Ser2, "X"};
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

const std::map<std::string, unsigned> SequenceSummary::aaToIndex = {{"A", 0}, {"C", 1}, {"D", 2}, {"E", 3}, {"F", 4}, {"G", 5}, {"H", 6}, {"I", 7},
	{"K", 8}, {"L", 9}, {"M", 10}, {"N", 11}, {"P", 12}, {"Q", 13}, {"R", 14},
	{"S", 15}, {"T", 16}, {"V", 17}, {"W", 18}, {"Y", 19}, {SequenceSummary::Ser2, 20}, {"X", 21}};

const std::map<std::string, unsigned> SequenceSummary::codonToIndexWithReference = {{"GCA", 0}, {"GCC", 1}, {"GCG", 2}, {"GCT", 3},
	{"TGC", 4}, {"TGT", 5}, {"GAC", 6}, {"GAT", 7}, {"GAA", 8}, {"GAG", 9}, {"TTC", 10}, {"TTT", 11}, {"GGA", 12},
	{"GGC", 13}, {"GGG", 14}, {"GGT", 15}, {"CAC", 16}, {"CAT", 17}, {"ATA", 18}, {"ATC", 19}, {"ATT", 20}, {"AAA", 21},
	{"AAG", 22}, {"CTA", 23}, {"CTC", 24}, {"CTG", 25}, {"CTT", 26}, {"TTA", 27}, {"TTG", 28}, {"ATG", 29}, {"AAC", 30},
	{"AAT", 31}, {"CCA", 32}, {"CCC", 33}, {"CCG", 34}, {"CCT", 35}, {"CAA", 36}, {"CAG", 37}, {"AGA", 38}, {"AGG", 39},
	{"CGA", 40}, {"CGC", 41}, {"CGG", 42}, {"CGT", 43}, {"TCA", 44}, {"TCC", 45}, {"TCG", 46}, {"TCT", 47}, {"ACA", 48},
	{"ACC", 49}, {"ACG", 50}, {"ACT", 51}, {"GTA", 52}, {"GTC", 53}, {"GTG", 54}, {"GTT", 55}, {"TGG", 56}, {"TAC", 57},
	{"TAT", 58}, {"AGC", 59}, {"AGT", 60}, {"TAA", 61}, {"TAG", 62}, {"TGA", 63}};

const std::map<std::string, unsigned> SequenceSummary::codonToIndexWithoutReference = {{"GCA", 0}, {"GCC", 1}, {"GCG", 2},
	{"TGC", 3}, {"GAC", 4}, {"GAA", 5}, {"TTC", 6}, {"GGA", 7},
 	{"GGC", 8}, {"GGG", 9}, {"CAC", 10}, {"ATA", 11}, {"ATC", 12}, {"AAA", 13},
	{"CTA", 14}, {"CTC", 15}, {"CTG", 16}, {"CTT", 17}, {"TTA", 18}, {"AAC", 19},
	{"CCA", 20}, {"CCC", 21}, {"CCG", 22}, {"CAA", 23}, {"AGA", 24}, {"AGG", 25},
	{"CGA", 26}, {"CGC", 27}, {"CGG", 28}, {"TCA", 29}, {"TCC", 30}, {"TCG", 31}, {"ACA", 32},
	{"ACC", 33}, {"ACG", 34}, {"GTA", 35}, {"GTC", 36}, {"GTG", 37}, {"TAC", 38},
	{"AGC", 39}};

void SequenceSummary::AAIndexToCodonRange(unsigned aaIndex, bool forParamVector, unsigned aaRange[])
{
	unsigned startAAIndex = 0;
	unsigned endAAIndex = 0;
	std::string aa = IndexToAA(aaIndex);

	if (aa == "A") 
	{
		if (!forParamVector) {startAAIndex = 0; endAAIndex = 4;}
		else { startAAIndex = 0; endAAIndex = 3;} 
	}
	else if (aa == "C") 
	{
		if (!forParamVector) {startAAIndex = 4; endAAIndex = 6;}
		else { startAAIndex = 3; endAAIndex = 4;} 
	}
	else if (aa == "D") 
	{
		if (!forParamVector) {startAAIndex = 6; endAAIndex = 8;}
		else { startAAIndex = 4; endAAIndex = 5;} 
	}
	else if (aa == "E") 
	{
		if (!forParamVector) {startAAIndex = 8; endAAIndex = 10;}
		else { startAAIndex = 5; endAAIndex = 6;} 
	}
	else if (aa == "F") 
	{
		if (!forParamVector) {startAAIndex = 10; endAAIndex = 12;}
		else { startAAIndex = 6; endAAIndex = 7;} 
	}
	else if (aa == "G") 
	{
		if (!forParamVector) {startAAIndex = 12; endAAIndex = 16;}
		else { startAAIndex = 7; endAAIndex = 10;} 
	}
	else if (aa == "H") 
	{
		if (!forParamVector) {startAAIndex = 16; endAAIndex = 18;}
		else { startAAIndex = 10; endAAIndex = 11;} 
	}
	else if (aa == "I") 
	{
		if (!forParamVector) {startAAIndex = 18; endAAIndex = 21;}
		else { startAAIndex = 11; endAAIndex = 13;} 
	}
	else if (aa == "K") 
	{
		if (!forParamVector) {startAAIndex = 21; endAAIndex = 23;}
		else { startAAIndex = 13; endAAIndex = 14;} 
	}
	else if (aa == "L") 
	{
		if (!forParamVector) {startAAIndex = 23; endAAIndex = 29;}
		else { startAAIndex = 14; endAAIndex = 19;} 
	}
	else if (aa == "M") 
	{
		if (!forParamVector) {startAAIndex = 29; endAAIndex = 30;}
		else { startAAIndex = 19; endAAIndex = 19;} 
	}
	else if (aa == "N") 
	{
		if (!forParamVector) {startAAIndex = 30; endAAIndex = 32;}
		else { startAAIndex = 19; endAAIndex = 20;} 
	}
	else if (aa == "P") 
	{
		if (!forParamVector) {startAAIndex = 32; endAAIndex = 36;}
		else { startAAIndex = 20; endAAIndex = 23;} 
	}
	else if (aa == "Q") 
	{
		if (!forParamVector) {startAAIndex = 36; endAAIndex = 38;}
		else { startAAIndex = 23; endAAIndex = 24;} 
	}
	else if (aa == "R") 
	{
		if (!forParamVector) {startAAIndex = 38; endAAIndex = 44;}
		else { startAAIndex = 24; endAAIndex = 29;} 
	}
	else if (aa == "S") 
	{
		if (!forParamVector) {startAAIndex = 44; endAAIndex = 48;}
		else { startAAIndex = 29; endAAIndex = 32;} 
	}
	else if (aa == "T") 
	{
		if (!forParamVector) {startAAIndex = 48; endAAIndex = 52;}
		else { startAAIndex = 32; endAAIndex = 35;} 
	}
	else if (aa == "V") 
	{
		if (!forParamVector) {startAAIndex = 52; endAAIndex = 56;}
		else { startAAIndex = 35; endAAIndex = 38;} 
	}
	else if (aa == "W") 
	{
		if (!forParamVector) {startAAIndex = 56; endAAIndex = 57;}
		else { startAAIndex = 38; endAAIndex = 38;} 
	}
	else if (aa == "Y") 
	{
		if (!forParamVector) {startAAIndex = 57; endAAIndex = 59;}
		else { startAAIndex = 38; endAAIndex = 39;} 
	}
	else if (aa == "Z") 
	{
		if (!forParamVector) {startAAIndex = 59; endAAIndex = 61;}
		else { startAAIndex = 39; endAAIndex = 40;} 
	}
	else if (aa == "X") 
	{
		if (!forParamVector) {startAAIndex = 61; endAAIndex = 64;}
		else { startAAIndex = 40; endAAIndex = 40;} 
	}
	else //Invalid AA
	{
		startAAIndex = 0;
		endAAIndex = 0;
	}

	aaRange[0] = startAAIndex;
	aaRange[1] = endAAIndex;
}

void SequenceSummary::AAToCodonRange(std::string aa, bool forParamVector, unsigned aaRange[])
{
	unsigned aaIndex = aaToIndex.find(aa) -> second; 
	AAIndexToCodonRange(aaIndex, forParamVector, aaRange);
}

std::string SequenceSummary::IndexToAA(int aa)
{
	return AminoAcidArray[aa];
}


unsigned SequenceSummary::CodonToAAIndex(std::string& codon)
{
	std::string aa = SequenceSummary::CodonToAA(codon);
	return SequenceSummary::aaToIndex.find(aa) -> second;
}

unsigned SequenceSummary::CodonToIndex(std::string& codon, bool forParamVector)
{
	unsigned i = 0;
	codon[0] = (char) std::toupper(codon[0]);
	codon[1] = (char) std::toupper(codon[1]);
	codon[2] = (char) std::toupper(codon[2]);
	if(forParamVector)
	{
		i = SequenceSummary::codonToIndexWithoutReference.find(codon) -> second;
	}else{
		i = SequenceSummary::codonToIndexWithReference.find(codon) -> second;
	}
	return i;
}

std::string SequenceSummary::IndexToCodon(unsigned i, bool forParamVector)
{
	return forParamVector ? codonArrayParameter[i] : codonArray[i];
}

unsigned SequenceSummary::AAToAAIndex(std::string aa)
{
	return SequenceSummary::aaToIndex.find(aa) -> second;
}

unsigned SequenceSummary::GetNumCodonsForAA(std::string& aa, bool forParamVector)
{
	aa[0] = (char) std::toupper(aa[0]);
	unsigned ncodon = 0;
	if(aa == "M" || aa == "W") ncodon = 1;
	else if(aa == "C" || aa == "D" || aa == "E" || aa == "F" || aa == "H" || aa == "K" || aa == "N" || aa == "Q" || aa == Ser2 || aa == "Y") ncodon = 2;
	else if(aa == "I" || aa == "X") ncodon = 3;
	else if(aa == "A" || aa == "G" || aa == "P" || aa == "S" || aa == "T" || aa == "V") ncodon = 4;
	else if(aa == "L" || aa == "R") ncodon = 6;

	return (forParamVector ? (ncodon - 1) : ncodon);
}


std::string SequenceSummary::CodonToAA(std::string& codon)
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

std::vector<std::string> SequenceSummary::AAToCodon(std::string aa, bool forParamVector)
{
	std::vector <std::string> RV;
	unsigned aaRange[2];
	aa = (char) std::toupper(aa[0]);

	AAToCodonRange(aa, forParamVector, aaRange);
	if(forParamVector){
		for (unsigned i = aaRange[0]; i < aaRange[1]; i++)
		{
			RV.push_back(codonArrayParameter[i]);
		}
	}else{
		for (unsigned i = aaRange[0]; i < aaRange[1]; i++)
		{
			RV.push_back(codonArray[i]);
		}
	}
	return RV;
}





// ---------------------------------------------------------------------------
// ----------------------------- RCPP STUFF ----------------------------------
// ---------------------------------------------------------------------------
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
RCPP_MODULE(SequenceSummary_mod)
{
	class_<SequenceSummary>( "SequenceSummary" )
		.constructor("empty constructor")
		.constructor<std::string>("Initialize with a DNA Sequence. Sequence must be a multiple of 3")

		.method("getAACountForAA", &SequenceSummary::getAACountForAAR, "returns occurrence of a given amino acid in a sequence")
		.method("getAACountForAAIndex", &SequenceSummary::getAACountForAAIndexR)
		.method("getCodonCountForCodon", &SequenceSummary::getCodonCountForCodonR, "returns occurrence of given codon in sequence")
		.method("getCodonCountForCodonIndex", &SequenceSummary::getCodonCountForCodonIndexR, "returns occurrence of given codon in sequence")
		.method("getRFPObserved", &SequenceSummary::getRFPObservedR)
		.method("setRFPObserved", &SequenceSummary::setRFPObservedR)
		.method("getCodonPositions", &SequenceSummary::getCodonPositionsR)
		.method("clear", &SequenceSummary::clear, "removes all data from object")
		.method("processSequence", &SequenceSummary::processSequence, "generates codon and amino acid count for sequence")



		;

	//static class functions
	function("CodonToAA", &SequenceSummary::CodonToAA, List::create(_["codon"]), "returns an amino acid for a given codon");
	function("GetNumCodonsForAA", &SequenceSummary::GetNumCodonsForAA, List::create(_["aa"], _["forParamVector"] = false), "returns the number of codons for a given amino acid");
	function("AAToCodon", &SequenceSummary::AAToCodon, List::create(_["aa"], _["forParamVector"] = false), "returns a vector of codons for a given amino acid");
	function("CodonToIndex", &SequenceSummary::CodonToIndex, List::create(_["codon"], _["forParamVector"] = false));


	//R wrapper function

	function("AAToAAIndex", &SequenceSummary::AAToAAIndexR);

		function("aminoAcids", &SequenceSummary::aminoAcids, "returns all Amino Acids as one letter code");
	function("codons", &SequenceSummary::codons, "returns all codons or all reference codons");
}
#endif

