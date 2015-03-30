#include "../include/SequenceSummary.h"

#include <ctype.h>
#include <algorithm>


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
    std::copy(std::begin(other.ncodons), std::end(other.ncodons), std::begin(ncodons));
    std::copy(std::begin(other.naa), std::end(other.naa), std::begin(naa));
}

SequenceSummary& SequenceSummary::operator=(const SequenceSummary& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    std::copy(std::begin(rhs.ncodons), std::end(rhs.ncodons), std::begin(ncodons));
    std::copy(std::begin(rhs.naa), std::end(rhs.naa), std::begin(naa));
    //assignment operator
    return *this;
}

void SequenceSummary::clear()
{
    // Why is "ncodons[64] = {};" not working?
    for(int k = 0; k < 64; k++) { ncodons[k] = 0; }
    for(int k = 0; k < 22; k++) { naa[k] = 0; }
}

void SequenceSummary::processSequence(const std::string& sequence)
{
    for(int i = 0; i < sequence.length(); i+=3)
    {
        std::string codon = sequence.substr(i, 3);
        int codonID = SequenceSummary::CodonToIndex( codon );
        int aaID = SequenceSummary::CodonToAAIndex( codon );
        ncodons[codonID]++;
        naa[aaID]++;
    }
}

/*
* STATIC FUNCTIONS
*/
const char SequenceSummary::AminoAcidArray[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', SequenceSummary::Ser2, 'X'};
const std::string SequenceSummary::codonArray[] = {"GCA", "GCC", "GCG", "GCT", "TGC", "TGT", "GAC", "GAT", "GAA", "GAG",
                                                    "TTC", "TTT", "GGA", "GGC", "GGG", "GGT", "CAC", "CAT", "ATA", "ATC",
                                                    "ATT", "AAA", "AAG", "CTA", "CTC", "CTG", "CTT", "TTA", "TTG", "ATG",
                                                    "AAC", "AAT", "CCA", "CCC", "CCG", "CCT", "CAA", "CAG", "AGA", "AGG",
                                                    "CGA", "CGC", "CGG", "CGT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC",
                                                    "ACG", "ACT", "GTA", "GTC", "GTG", "GTT", "TGG", "TAC", "TAT", "AGT",
                                                    "AGC", "TAA", "TAG", "TGA"};

const std::map<char, int> SequenceSummary::aaToIndex = {{'A', 0}, {'C', 1}, {'D', 2}, {'E', 3}, {'F', 4}, {'G', 5}, {'H', 6}, {'I', 7},
    {'K', 8}, {'L', 9}, {'M', 10}, {'N', 11}, {'P', 12}, {'Q', 13}, {'R', 14},
    {'S', 15}, {'T', 16}, {'V', 17}, {'W', 18}, {'Y', 19}, {SequenceSummary::Ser2, 20}, {'X', 21}};

unsigned* SequenceSummary::AAindexToCodonRange(unsigned aaIndex, bool forParamVector)
{
    unsigned startAAIndex = 0;
    unsigned numCodons = 4;

    for(unsigned i = 0; i <= aaIndex; i++)
    {
        char aa = AminoAcidArray[i];
        numCodons = GetNumCodonsForAA(aa, forParamVector);
        if(i != aaIndex)
        {
            startAAIndex += numCodons;
        }
    }
    unsigned endAAIndex = numCodons + startAAIndex;

    unsigned ret[2];
    ret[0] = startAAIndex;
    ret[1] = endAAIndex;
    return ret;
}
unsigned* SequenceSummary::AAToCodonRange(char aa, bool forParamVector)
{
    unsigned aaIndex = aaToIndex.find(aa)->second;
    return AAindexToCodonRange(aaIndex, forParamVector);
}

char SequenceSummary::IndexToAA(int aa)
{
    return AminoAcidArray[aa];
}


unsigned SequenceSummary::CodonToAAIndex(std::string& codon)
{
    char aa = SequenceSummary::CodonToAA(codon);
    return SequenceSummary::aaToIndex.find(aa)->second;
}

unsigned SequenceSummary::CodonToIndex(std::string& codon)
{
    //std::transform(codon.begin(), codon.end(), codon.begin(), ::toupper);
    int i = 0;
    for(; i < 64; i++)
    {
        if(!codon.compare(SequenceSummary::codonArray[i])) break;
    }
    return i;
}
std::string SequenceSummary::IndexToCodon(unsigned i)
{
    return SequenceSummary::codonArray[i];
}


unsigned SequenceSummary::GetNumCodonsForAA(const char& aa, bool forParamVector)
{
    std::toupper(aa);
    unsigned ncodon = -1;
    if(aa == 'M' || aa == 'W') ncodon = 1;
    else if(aa == 'C' || aa == 'D' || aa == 'E' || aa == 'F' || aa == 'H' || aa == 'K' || aa == 'N' || aa == 'Q' || aa == Ser2 || aa == 'Y') ncodon = 2;
    else if(aa == 'I' || aa == 'X') ncodon = 3;
    else if(aa == 'A' || aa == 'G' || aa == 'P' || aa == 'S' || aa == 'T' || aa == 'V') ncodon = 4;
    else if(aa == 'L' || aa == 'R') ncodon = 6;

    return (forParamVector ? (ncodon - 1) : ncodon);
}


char SequenceSummary::CodonToAA(std::string& codon)
{
    //std::transform(codon.begin(), codon.end(), codon.begin(), ::toupper);
    char aa = '#';
    //Phenylalanine
    if(!codon.compare("TTT") || !codon.compare("UUU") || !codon.compare("TTC") || !codon.compare("UUC")) aa = 'F';
    //Leucine
    else if(!codon.compare("TTA") || !codon.compare("UUA") || !codon.compare("TTG") || !codon.compare("UUG") ||
    !codon.compare("CTT") || !codon.compare("CUU") || !codon.compare("CTC") || !codon.compare("CUC") ||
    !codon.compare("CTA") || !codon.compare("CUA") || !codon.compare("CTG") || !codon.compare("CUG")) aa = 'L';
    //Isoleucine
    else if(!codon.compare("ATT") || !codon.compare("AUU") || !codon.compare("ATC") || !codon.compare("AUC") ||
    !codon.compare("ATA") || !codon.compare("AUA")) aa = 'I';
    //Methionine
    else if(!codon.compare("ATG") || !codon.compare("AUG")) aa = 'M';
    //Valine
    else if(!codon.compare("GTT") || !codon.compare("GUU") || !codon.compare("GTC") || !codon.compare("GUC") ||
    !codon.compare("GTA") || !codon.compare("GUA") || !codon.compare("GTG") || !codon.compare("GUG")) aa = 'V';
    //Serine4
    else if(!codon.compare("TCT") || !codon.compare("UCU") || !codon.compare("TCC") || !codon.compare("UCC") ||
    !codon.compare("TCA") || !codon.compare("UCA") || !codon.compare("TCG") || !codon.compare("UCG")) aa = 'S';
    //Proline
    else if(!codon.compare("CCT") || !codon.compare("CCU") || !codon.compare("CCC") ||
    !codon.compare("CCA") || !codon.compare("CCG")) aa = 'P';
    //Threonine
    else if(!codon.compare("ACT") || !codon.compare("ACU") || !codon.compare("ACC") ||
    !codon.compare("ACA") || !codon.compare("ACG")) aa = 'T';
    //Alanine
    else if(!codon.compare("GCT") || !codon.compare("GCU") || !codon.compare("GCC") ||
    !codon.compare("GCA") || !codon.compare("GCG")) aa = 'A';
    //Tyrosine
    else if(!codon.compare("TAT") || !codon.compare("UAU") || !codon.compare("TAC") || !codon.compare("UAC")) aa = 'Y';
    //Histidine
    else if(!codon.compare("CAT") || !codon.compare("CAU") || !codon.compare("CAC")) aa = 'H';
    //Glutamine
    else if(!codon.compare("CAA") || !codon.compare("CAG")) aa = 'Q';
    //Asparagine
    else if(!codon.compare("AAT") || !codon.compare("AAU") || !codon.compare("AAC")) aa = 'N';
    //Lysine
    else if(!codon.compare("AAA") || !codon.compare("AAG")) aa = 'K';
    //Aspartic Acid
    else if(!codon.compare("GAT") || !codon.compare("GAU") || !codon.compare("GAC")) aa = 'D';
    //Glutamic Acid
    else if(!codon.compare("GAA") || !codon.compare("GAG")) aa = 'E';
    //Cysteine
    else if(!codon.compare("TGT") || !codon.compare("TAT") || !codon.compare("TGC") || !codon.compare("UGC")) aa = 'C';
    //Tryptophan
    else if(!codon.compare("TGG") || !codon.compare("UGG")) aa = 'W';
    //Arginine
    else if(!codon.compare("CGT") || !codon.compare("CGU") || !codon.compare("CGC") || !codon.compare("CGA") ||
    !codon.compare("CGG") || !codon.compare("AGA") || !codon.compare("AGG")) aa = 'R';
    //Serine2
    else if(!codon.compare("AGT") || !codon.compare("AGU") || !codon.compare("AGC")) aa = SequenceSummary::Ser2;
    //Glycine
    else if(!codon.compare("GGT") || !codon.compare("GGU") || !codon.compare("GGC")  || !codon.compare("GGC") ||
    !codon.compare("GGA")  || !codon.compare("GGG")) aa = 'G';
    //Stop
    else if (!codon.compare("TAA") || !codon.compare("UAA") || !codon.compare("TAG") || !codon.compare("UAG") ||
    !codon.compare("TGA") || !codon.compare("UGA")) aa = 'X';

    return aa;
}
