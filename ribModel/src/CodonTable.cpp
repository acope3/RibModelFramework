#include "include/CodonTable.h"

#include <iostream>


//----------------------------------------------//
//----------Constructors & Destructors----------//
//----------------------------------------------//
CodonTable::CodonTable()
{
    tableId = 0; //standard codon table by NCBI
    splitAA = true;
}


CodonTable::CodonTable(unsigned _tableId, bool _splitAA) : tableId(_tableId), splitAA(_splitAA)
{
	if(tableId == 7 || tableId == 8 || tableId == 15 || tableId == 17 || tableId == 18 || tableId == 19 || tableId == 20)
	{
        std::cerr << "Invalid codon table: " << tableId << " using default codon table (NCBI codon table 1)\n";
        tableId = 1; //standard codon table by NCBI
    }
    tableId--; //Make it 0 indexed
}


CodonTable::~CodonTable()
{
    //dtor
}


CodonTable::CodonTable(const CodonTable& other)
{
    tableId = other.tableId;
    splitAA = other.splitAA;
}


CodonTable& CodonTable::operator=(const CodonTable& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    tableId = rhs.tableId;
    splitAA = rhs.splitAA;
    return *this;
}



//------------------------------------//
//----------Getter Functions----------//
//------------------------------------//
unsigned  CodonTable::getTableId()
{
    return tableId;
}


bool CodonTable::getSplitAA()
{
    return splitAA;
}


std::vector<std::vector<unsigned>> CodonTable::getCodonIndexListing()
{
    return codonIndexListing;
}


std::vector<std::vector<unsigned>> CodonTable::getCodonIndexListingWithoutReference()
{
    return codonIndexListingWithoutReference;
}


std::vector <std::string> CodonTable::getAAListing()
{
    return AAListing;
}


std::vector <std::string> CodonTable::getForParamVectorListing()
{
    return forParamVectorListing;
}


std::map <std::string, std::string> CodonTable::getCodonToAAMap()
{
    return codonToAAMap;
}


std::map <std::string, unsigned> CodonTable::getAAMap()
{
    return AAMap;
}


std::map <std::string, unsigned> CodonTable::getAAToNumCodonsMap()
{
    return AAToNumCodonsMap;
}


std::map <std::string, unsigned> CodonTable::getForParamVectorMap()
{
    return forParamVectorMap;
}


unsigned CodonTable::getNumCodonsForAA(std::string aa, bool forParamVector)
{
    unsigned numCodons = AAToNumCodonsMap.find(aa) -> second;
    return forParamVector ? numCodons - 1 : numCodons;
}


unsigned CodonTable::getNumCodonsForAAIndex(unsigned aaIndex, bool forParamVector)
{
    std::string aa = indexToAA(aaIndex);
    return getNumCodonsForAA(aa, forParamVector);
}


std::string CodonTable::getForParamVectorCodon(unsigned codonIndex)
{
    return forParamVectorListing[codonIndex];
}



//--------------------------------------//
//----------Mapping Operations----------//
//--------------------------------------//
unsigned CodonTable::AAToAAIndex(std::string aa)
{
        return AAMap.find(aa) -> second;
}


std::vector <unsigned> CodonTable::AAIndexToCodonRange(unsigned aaIndex, bool forParamVector)
{
    return forParamVector ? codonIndexListingWithoutReference[aaIndex] : codonIndexListing[aaIndex];
}


std::string CodonTable::indexToCodon(unsigned index, bool forParamVector)
{
    return forParamVector ? forParamVectorListing[index] : CodonTable::codonArray[index];
}


std::vector <unsigned> CodonTable::AAToCodonRange(std::string aa, bool forParamVector)
{
    unsigned aaIndex = AAToAAIndex(aa);
    return AAIndexToCodonRange(aaIndex, forParamVector);
}


std::vector<std::string> CodonTable::AAToCodon(std::string aa, bool forParamVector)
{
    std::vector <std::string> codons;
    std::vector <unsigned> aaRange = AAToCodonRange(aa, forParamVector);
    if (forParamVector)
    {
        for (unsigned i = 0; i < aaRange.size(); i++)
        {
            codons.push_back(indexToCodon(aaRange[i]));
        }
    }
    else
    {
        for (unsigned i = 0; i < aaRange.size(); i++)
        {
            codons.push_back(indexToCodon(aaRange[i]));
        }
    }
    return codons;
}


std::string CodonTable::codonToAA(std::string& codon)
{
    std::map <std::string, std::string>::iterator mit;

    mit = codonToAAMap.find(codon);
    return mit -> second;
}


unsigned CodonTable::codonToIndex(std::string& codon, bool forParamVector)
{
    return forParamVector ? forParamVectorMap.find(codon) -> second : codonToIndexWithReference.find(codon) -> second;
}


unsigned CodonTable::codonToAAIndex(std::string& codon)
{
    std::string AA = codonToAA(codon);
    return AAMap.find(AA) -> second;
}


std::string CodonTable::indexToAA(unsigned aaIndex)
{
    return AAListing[aaIndex];
}



//-----------------------------------//
//----------Other Functions----------//
//-----------------------------------//
void CodonTable::setupCodonTable()
{
	unsigned numAA = 21;
	if(tableId >= 0 && tableId <= 5  && !splitAA) numAA = 20; //If not splitting AAs, tables 1-6 all have 20 AA, excluding stop codes
    else if (tableId >= 8 && tableId <= 13  && !splitAA) numAA = 20;
    else if (tableId == 15 && !splitAA) numAA = 20;
    else if (tableId >= 20 && tableId <= 24 && !splitAA) numAA = 20;


    numAA++; //Account for the stop codon
	codonIndexListing.resize(numAA);

    unsigned aaIndex = 0;
    unsigned filled = 0;
    while (filled != numAA)
    {
        unsigned numCodons = numCodonsPerAAForTable[tableId][aaIndex];
        if (numCodons != 0)
        {
            if (aaIndex == 0) //A
            {
                for (unsigned i = 0; i < 4; i++)
                {
                    codonIndexListing[filled].push_back(i);
                }
                AAToNumCodonsMap.insert(std::make_pair("A", 4));
            }
            else if (aaIndex == 1) //C
            {
                AAToNumCodonsMap.insert(std::make_pair("C", 2));
                for (unsigned i = 4; i < 6; i++)
                {
                    codonIndexListing[filled].push_back(i);
                }
                if (tableId == 10)
                {
                    codonIndexListing[filled].push_back(63);
                    AAToNumCodonsMap["C"]++;
                }

            }
            else if (aaIndex == 2) //D
            {
                AAToNumCodonsMap.insert(std::make_pair("D", 2));
                for (unsigned i = 6; i < 8; i++)
                {
                    codonIndexListing[filled].push_back(i);
                }
            }
            else if (aaIndex == 3) //E
            {
                AAToNumCodonsMap.insert(std::make_pair("E", 2));
                for (unsigned i = 8; i < 10; i++)
                {
                    codonIndexListing[filled].push_back(i);
                }
            }
            else if (aaIndex == 4) //F
            {
                AAToNumCodonsMap.insert(std::make_pair("F", 2));
                for (unsigned i = 10; i < 12; i++)
                {
                    codonIndexListing[filled].push_back(i);
                }
            }
            else if (aaIndex == 5) //G
            {
                AAToNumCodonsMap.insert(std::make_pair("G", 4));
                for (unsigned i = 12; i < 16; i++)
                {
                    codonIndexListing[filled].push_back(i);
                }
                if (tableId == 12)
                {
                    AAToNumCodonsMap["G"] += 2;
                    codonIndexListing[filled].push_back(38);
                    codonIndexListing[filled].push_back(39);
                }
                else if (tableId == 24)
                {
                    AAToNumCodonsMap["G"]++;
                    codonIndexListing[filled].push_back(63);
                }
            }
            else if (aaIndex == 6) //H
            {
                AAToNumCodonsMap.insert(std::make_pair("H", 2));
                for (unsigned i = 16; i < 18; i++)
                {
                    codonIndexListing[filled].push_back(i);
                }
            }
            else if (aaIndex == 7) //I
            {
                AAToNumCodonsMap.insert(std::make_pair("I", 2));
                if (tableId == 0 || tableId == 3 || tableId == 5 || tableId == 8 || tableId == 9 || tableId == 10
                        || tableId == 11 || tableId == 13 || tableId == 15 || tableId >= 21)
                {
                    AAToNumCodonsMap["I"]++;
                    codonIndexListing[filled].push_back(18);
                }
                for (unsigned i = 19; i < 21; i++)
                {
                    codonIndexListing[filled].push_back(i);
                }
            }
            else if (aaIndex == 8) //K
            {

                if ((tableId >= 0 && tableId <= 5) || (tableId >= 9 && tableId <= 12) || (tableId == 15) ||
                        (tableId >= 21))
                {
                    AAToNumCodonsMap.insert(std::make_pair("K", 2));
                    for (unsigned i = 21; i < 23; i++)
                    {
                        codonIndexListing[filled].push_back(i);
                    }
                }
                else if (tableId == 8 || tableId == 13 || tableId == 20)
                {
                    AAToNumCodonsMap.insert(std::make_pair("K", 1));
                    codonIndexListing[filled].push_back(22);
                }

                if (tableId == 23)
                {
                    AAToNumCodonsMap["K"]++;
                    codonIndexListing[filled].push_back(39);
                }
            }
            else if (aaIndex == 9) //L
            {
                if (tableId != 2 && tableId != 22)
                {
                    if (tableId == 11) AAToNumCodonsMap.insert(std::make_pair("L", 5));
                    else AAToNumCodonsMap.insert(std::make_pair("L", 6));
                    for (unsigned i = 23; i < 29; i++)
                    {
                        if (i == 25 && tableId != 11)
                        {
                            codonIndexListing[filled].push_back(i);
                        }
                        else
                        {
                            codonIndexListing[filled].push_back(i);
                        }
                    }
                    if ((tableId == 15 || tableId == 21) && !splitAA)
                    {
                        AAToNumCodonsMap["L"]++;
                        codonIndexListing[filled].push_back(62);
                    }
                }
                else if (tableId == 2)
                {
                    AAToNumCodonsMap.insert(std::make_pair("L", 2));
                    for (unsigned i = 27; i < 29; i++)
                    {
                        codonIndexListing[filled].push_back(i);
                    }
                }
                else if (tableId == 22)
                {
                    AAToNumCodonsMap.insert(std::make_pair("L", 5));
                    for (unsigned i = 23; i < 27; i++)
                    {
                        codonIndexListing[filled].push_back(i);
                    }
                    codonIndexListing[filled].push_back(28);
                }
            }
            else if (aaIndex == 10) //Leu1
            {
                if (splitAA)
                {
                    AAToNumCodonsMap.insert(std::make_pair(CodonTable::Leu1, 1));
                    codonIndexListing[filled].push_back(62);
                }
                else
                {
                    filled--;
                }
            }
            else if (aaIndex == 11) //M
            {
                AAToNumCodonsMap.insert(std::make_pair("M", 1));
                codonIndexListing[filled].push_back(29);
                if (tableId == 1 || tableId == 2 || tableId == 4 || tableId == 12 || tableId == 20)
                {
                    AAToNumCodonsMap["M"]++;
                    codonIndexListing[filled].push_back(18);
                }
            }
            else if (aaIndex == 12) //N
            {
                AAToNumCodonsMap["N"] += 2;
                for (unsigned i = 30; i < 32; i++)
                {
                    codonIndexListing[filled].push_back(i);
                }
                if (tableId == 8 || tableId == 13 || tableId == 20)
                {
                    AAToNumCodonsMap["N"]++;
                    codonIndexListing[filled].push_back(21);
                }
            }
            else if (aaIndex == 13) //P
            {
                AAToNumCodonsMap.insert(std::make_pair("P", 4));
                for (unsigned i = 32; i < 36; i++)
                {
                    codonIndexListing[filled].push_back(i);
                }
            }
            else if (aaIndex == 14) //Q
            {
                AAToNumCodonsMap.insert(std::make_pair("Q", 2));
                for (unsigned i = 36; i < 38; i++)
                {
                    codonIndexListing[filled].push_back(i);
                }
                if (tableId == 5)
                {
                    AAToNumCodonsMap["Q"] += 2;
                    codonIndexListing[filled].push_back(61);
                    codonIndexListing[filled].push_back(62);
                }
            }
            else if (aaIndex == 15) //R
            {
                AAToNumCodonsMap.insert(std::make_pair("R", 0));
                if (tableId == 0 || tableId == 2 || tableId == 3 || tableId == 5 || (tableId >= 9 && tableId <= 11) || tableId == 15 ||
                        tableId == 21 || tableId == 22 || tableId == 24)
                {
                    AAToNumCodonsMap["R"] += 2;
                    codonIndexListing[filled].push_back(38);
                    codonIndexListing[filled].push_back(39);
                }

                if (tableId == 2)
                {
                    AAToNumCodonsMap["R"] += 2;
                    codonIndexListing[filled].push_back(42);
                    codonIndexListing[filled].push_back(43);
                }
                else
                {
                    AAToNumCodonsMap["R"] += 4;
                    for (unsigned i = 40; i < 44; i++)
                    {
                        codonIndexListing[filled].push_back(i);
                    }
                }
            }
            else if (aaIndex == 16) //Ser1
            {
                if(splitAA)
                {
                    AAToNumCodonsMap.insert(std::make_pair(CodonTable::Ser1, 1));
                    codonIndexListing[filled].push_back(25); //table 12
                }
                else
                {
                    filled--;
                }
            }
            else if (aaIndex == 17) //Ser2
            {
                if (splitAA)
                {
                    AAToNumCodonsMap.insert(std::make_pair(CodonTable::Ser2, 2));
                    codonIndexListing[filled].push_back(59);
                    codonIndexListing[filled].push_back(60);

                    if (tableId == 4 || tableId == 8 || tableId == 13 || tableId == 20)
                    {
                        AAToNumCodonsMap[CodonTable::Ser2] += 2;
                        codonIndexListing[filled].push_back(38);
                        codonIndexListing[filled].push_back(39);
                    }
                    else if (tableId == 23)
                    {
                        AAToNumCodonsMap[CodonTable::Ser2]++;
                        codonIndexListing[filled].push_back(38);
                    }
                }
                else
                {
                    filled--;
                }
            }
            else if (aaIndex == 18) //Ser (Ser4)
            {
                AAToNumCodonsMap.insert(std::make_pair("S", 3));
                if (tableId != 21)
                {
                    AAToNumCodonsMap["S"]++;
                    for(unsigned i = 44; i < 48; i++)
                    {
                        codonIndexListing[filled].push_back(i);
                    }
                }
                else
                {
                    for(unsigned i = 45; i < 48; i++)
                    {
                        codonIndexListing[filled].push_back(i);
                    }
                }

                if (!splitAA)
                {
                    AAToNumCodonsMap["S"] += 2;
                    if (tableId == 11)
                    {
                        AAToNumCodonsMap["S"]++;
                        codonIndexListing[filled].push_back(25);
                    }

                    codonIndexListing[filled].push_back(59);
                    codonIndexListing[filled].push_back(60);

                    if (tableId == 4 || tableId == 8 || tableId == 13 || tableId == 20)
                    {
                        AAToNumCodonsMap["S"] += 2;
                        codonIndexListing[filled].push_back(38);
                        codonIndexListing[filled].push_back(39);
                    }
                    else if (tableId == 23)
                    {
                        AAToNumCodonsMap["S"]++;
                        codonIndexListing[filled].push_back(38);
                    }
                }
            }
            else if (aaIndex == 19) //T
            {
                AAToNumCodonsMap.insert(std::make_pair("T", 0));
                if (tableId != 2)
                {
                    AAToNumCodonsMap["T"] += 4;
                    for (unsigned i = 48; i < 52; i++)
                    {
                        codonIndexListing[filled].push_back(i);
                    }
                }
                else if (tableId == 2 && !splitAA)
                {
                    AAToNumCodonsMap["T"] += 8;
                    for (unsigned i = 48; i < 52; i++)
                    {
                        codonIndexListing[filled].push_back(i);
                    }

                    for (unsigned i = 23; i < 27; i++)
                    {
                        codonIndexListing[filled].push_back(i);
                    }
                }
            }
            else if (aaIndex == 20) //Thr1
            {
                if (splitAA)
                {
                    AAToNumCodonsMap.insert(std::make_pair(CodonTable::Thr4_1, 4));
                    for (unsigned i = 48; i < 52; i++)
                    {
                        codonIndexListing[filled].push_back(i);
                    }
                }
                else
                {
                    filled--;
                }
            }
            else if (aaIndex == 21) //Thr2
            {
                if (splitAA)
                {
                    AAToNumCodonsMap.insert(std::make_pair(CodonTable::Thr4_2, 4));
                    for (unsigned i = 23; i < 27; i++)
                    {
                        codonIndexListing[filled].push_back(i);
                    }
                }
                else
                {
                    filled--;
                }
            }
            else if (aaIndex == 22) //V
            {
                AAToNumCodonsMap.insert(std::make_pair("V", 4));
                for (unsigned i = 52; i < 56; i++)
                {
                    codonIndexListing[filled].push_back(i);
                }
            }
            else if (aaIndex == 23) //W
            {
                AAToNumCodonsMap.insert(std::make_pair("W", 1));
                codonIndexListing[filled].push_back(56);
                if ((tableId >= 1 && tableId <= 4) || (tableId == 8) || (tableId == 12) || (tableId == 13) || (tableId == 20))
                {
                    AAToNumCodonsMap["W"]++;
                    codonIndexListing[filled].push_back(63);
                }
            }
            else if (aaIndex == 24) //Y
            {
                AAToNumCodonsMap.insert(std::make_pair("Y", 2));
                for (unsigned i = 57; i < 59; i++)
                {
                    codonIndexListing[filled].push_back(i);
                }
                if (tableId == 13){
                    AAToNumCodonsMap["Y"]++;
                    codonIndexListing[filled].push_back(61);
                }
            }
            else if (aaIndex == 25) //X
            {
                AAToNumCodonsMap.insert(std::make_pair("X", 1));
                if (tableId == 5)
                {
                    codonIndexListing[filled].push_back(63);
                }
                else if (tableId == 13)
                {
                    codonIndexListing[filled].push_back(62);
                }
                else if (tableId == 15)
                {
                    codonIndexListing[filled].push_back(61);
                    codonIndexListing[filled].push_back(63);
                    AAToNumCodonsMap["X"]++;
                }
                else if ((tableId >= 2 && tableId <= 4) || tableId == 8 || tableId == 9 || tableId == 12 ||
                        tableId == 20 || tableId == 23 || tableId == 24)
                {
                    codonIndexListing[filled].push_back(61);
                    codonIndexListing[filled].push_back(62);
                    AAToNumCodonsMap["X"]++;
                }
                else if (tableId == 0 || tableId == 10 || tableId == 11)
                {
                    codonIndexListing[filled].push_back(61);
                    codonIndexListing[filled].push_back(63);
                    codonIndexListing[filled].push_back(62);
                    AAToNumCodonsMap["X"] += 2;
                }
                else if (tableId == 21)
                {
                    codonIndexListing[filled].push_back(61);
                    codonIndexListing[filled].push_back(63);
                    codonIndexListing[filled].push_back(44);
                    AAToNumCodonsMap["X"] += 2;
                }
                else if (tableId == 1)
                {
                    codonIndexListing[filled].push_back(61);
                    codonIndexListing[filled].push_back(62);
                    codonIndexListing[filled].push_back(38);
                    codonIndexListing[filled].push_back(39);
                    AAToNumCodonsMap["X"] += 3;
                }
                else if (tableId == 22)
                {
                    codonIndexListing[filled].push_back(61);
                    codonIndexListing[filled].push_back(63);
                    codonIndexListing[filled].push_back(62);
                    codonIndexListing[filled].push_back(27);
                    AAToNumCodonsMap["X"] += 3;
                }
            }
            filled++;
        }
        aaIndex++;
    }

    //fill without refernce mapping
    unsigned val = 0;
    codonIndexListingWithoutReference.resize(numAA);
    for (unsigned i = 0; i < numAA; i++)
    {
        for (unsigned j = 0; j < codonIndexListing[i].size() - 1; j++)
        {
            codonIndexListingWithoutReference[i].push_back(codonIndexListing[i][j]);
            unsigned codonIndex = codonIndexListing[i][j];
            std::string codon = codonArray[codonIndex];
            forParamVectorMap.insert(std::make_pair(codon, val));
            forParamVectorListing.push_back(codon);
            val++;
        }
    }

    //set up aamapping
    unsigned index = 0;
    if (splitAA)
    {
        for (unsigned i = 0; i < 26; i++)
        {
            if (numCodonsPerAAForTable[tableId][i] != 0)
            {
                AAListing.push_back(aminoAcidArray[i]);
                AAMap.insert(std::make_pair(aminoAcidArray[i], index));
                index++;
            }
        }
    }
    else
    {
        for (unsigned i = 0; i < 21; i++)
        {
            AAListing.push_back(aminoAcidArrayWithoutSplit[i]);
            AAMap.insert(std::make_pair(aminoAcidArrayWithoutSplit[i], index));
            index++;
        }
    }




    //Build a map for the Codon To AA
    for(unsigned AAIndex = 0; AAIndex < codonIndexListing.size(); AAIndex++)
    {
        std::string AA = AAListing[AAIndex];
        for (unsigned codonIndex = 0; codonIndex < codonIndexListing[AAIndex].size(); codonIndex++)
        {
            std::string codon = codonArray[codonIndexListing[AAIndex][codonIndex]];
            codonToAAMap.insert(std::make_pair(codon, AA));
        }
    }
}


bool CodonTable::checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound)
{
    bool check = false;
    if (lowerbound <= index && index <= upperbound)
    {
        check = true;
    }
    else
    {
        std::cerr <<"Error with the index\nGIVEN: " << index <<"\n";
        std::cerr <<"MUST BE BETWEEN:	" << lowerbound << " & " << upperbound <<"\n";
    }
    return check;
}



//----------------------------------------------//
//---------Static Variables & Functions---------//
//----------------------------------------------//
CodonTable* CodonTable::codonTable;
const std::string CodonTable::Ser2 = "Z";
const std::string CodonTable::Ser1 = "J";
const std::string CodonTable::Thr4_1 = "T";
const std::string CodonTable::Thr4_2 = "B";
const std::string CodonTable::Leu1 = "U";

const std::vector <std::string> CodonTable::aminoAcidArray = {
        "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", CodonTable::Leu1, "M", "N", "P", "Q", "R",
        CodonTable::Ser1, CodonTable::Ser2, "S", "T", CodonTable::Thr4_1, CodonTable::Thr4_2, "V", "W", "Y", "X"};
const std::vector <std::string> CodonTable::aminoAcidArrayWithoutSplit = {
        "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X"};


const std::map<std::string, unsigned> CodonTable::codonToIndexWithReference = {{"GCA", 0}, {"GCC", 1}, {"GCG", 2},
   {"GCT", 3}, {"TGC", 4}, {"TGT", 5}, {"GAC", 6}, {"GAT", 7}, {"GAA", 8}, {"GAG", 9}, {"TTC", 10}, {"TTT", 11},
   {"GGA", 12}, {"GGC", 13}, {"GGG", 14}, {"GGT", 15}, {"CAC", 16}, {"CAT", 17}, {"ATA", 18}, {"ATC", 19}, {"ATT", 20},
   {"AAA", 21}, {"AAG", 22}, {"CTA", 23}, {"CTC", 24}, {"CTG", 25}, {"CTT", 26}, {"TTA", 27}, {"TTG", 28}, {"ATG", 29},
   {"AAC", 30}, {"AAT", 31}, {"CCA", 32}, {"CCC", 33}, {"CCG", 34}, {"CCT", 35}, {"CAA", 36}, {"CAG", 37}, {"AGA", 38},
   {"AGG", 39}, {"CGA", 40}, {"CGC", 41}, {"CGG", 42}, {"CGT", 43}, {"TCA", 44}, {"TCC", 45}, {"TCG", 46}, {"TCT", 47},
   {"ACA", 48}, {"ACC", 49}, {"ACG", 50}, {"ACT", 51}, {"GTA", 52}, {"GTC", 53}, {"GTG", 54}, {"GTT", 55}, {"TGG", 56},
   {"TAC", 57}, {"TAT", 58}, {"AGC", 59}, {"AGT", 60}, {"TAA", 61}, {"TAG", 62}, {"TGA", 63}};


const std::string CodonTable::codonArray[] =
        {"GCA", "GCC", "GCG", "GCT", "TGC", "TGT", "GAC", "GAT", "GAA", "GAG",
         "TTC", "TTT", "GGA", "GGC", "GGG", "GGT", "CAC", "CAT", "ATA", "ATC",
         "ATT", "AAA", "AAG", "CTA", "CTC", "CTG", "CTT", "TTA", "TTG", "ATG",
         "AAC", "AAT", "CCA", "CCC", "CCG", "CCT", "CAA", "CAG", "AGA", "AGG",
         "CGA", "CGC", "CGG", "CGT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC",
         "ACG", "ACT", "GTA", "GTC", "GTG", "GTT", "TGG", "TAC", "TAT", "AGC",
         "AGT", "TAA", "TAG", "TGA"};


// http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
const std::vector <std::string> CodonTable::codonTableDefinition = {"1. The Standard Code",
   "2. The Vertebrate Mitochondrial Code", "3. The Yeast Mitochondrial Code",
   "4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code",
   "5. The Invertebrate Mitochondrial Code", "6. The Ciliate, Dasycladacean and Hexamita Nuclear Code",
   "7. Invalid Codon Table", "8. Invalid Codon Table", "9. The Echinoderm and Flatworm Mitochondrial Code",
   "10. The Euplotid Nuclear Code", "11. The Bacterial, Archaeal and Plant Plastid Code",
   "12. The Alternative Yeast Nuclear Code", "13. The Ascidian Mitochondrial Code",
   "14. The Alternative Flatworm Mitochondrial Code", "15. Invalid Codon Table", "16. Chlorophycean Mitochondrial Code",
   "17. Invalid Codon Table", "18. Invalid Codon Table", "19. Invalid Codon Table", "20. Invalid Codon Table",
   "21. Trematode Mitochondrial Code", "22. Scenedesmus obliquus Mitochondrial Code",
   "23. Thraustochytrium Mitochondrial Code", "24. Pterobranchia Mitochondrial Code",
   "25. Candidate Division SR1 and Gracilibacteria Code"};


const unsigned CodonTable::numCodonsPerAAForTable[25][26] = {
        {4,2,2,2,2,4,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,3}, // 1. The Standard Code
        {4,2,2,2,2,4,2,2,2,6,0,2,2,4,2,4,0,2,4,4,0,0,4,2,2,4}, // 2. The Vertebrate Mitochondrial Code
        {4,2,2,2,2,4,2,2,2,2,0,2,2,4,2,6,0,2,4,0,4,4,4,2,2,2}, // 3. The Yeast Mitochondrial Code
        {4,2,2,2,2,4,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,2,2,2}, // 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
        {4,2,2,2,2,4,2,2,2,6,0,2,2,4,2,4,0,4,4,4,0,0,4,2,2,2}, // 5. The Invertebrate Mitochondrial Code
        {4,2,2,2,2,4,2,3,2,6,0,1,2,4,4,6,0,2,4,4,0,0,4,1,2,1}, // 6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 7. Invalid Codon Table
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 8. Invalid Codon Table
        {4,2,2,2,2,4,2,3,1,6,0,1,3,4,2,4,0,4,4,4,0,0,4,2,2,2}, // 9. The Echinoderm and Flatworm Mitochondrial Code
        {4,3,2,2,2,4,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,2}, // 10. The Euplotid Nuclear Code
        {4,2,2,2,2,4,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,3}, // 11. The Bacterial, Archaeal and Plant Plastid Codee
        {4,2,2,2,2,4,2,3,2,5,0,1,2,4,2,6,1,2,4,4,0,0,4,1,2,3}, // 12. The Alternative Yeast Nuclear Code
        {4,2,2,2,2,6,2,2,2,6,0,2,2,4,2,4,0,2,4,4,0,0,4,2,2,2}, // 13. The Ascidian Mitochondrial Code
        {4,2,2,2,2,4,2,3,1,6,0,1,3,4,2,4,0,4,4,4,0,0,4,2,3,1}, // 14. The Alternative Flatworm Mitochondrial Code
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 15. Invalid Codon Table
        {4,2,2,2,2,4,2,3,2,6,1,1,2,4,2,6,0,2,4,4,0,0,4,1,2,2}, // 16. Chlorophycean Mitochondrial Code
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 17. Invalid Codon Table
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 18. Invalid Codon Table
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 19. Invalid Codon Table
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 20. Invalid Codon Table
        {4,2,2,2,2,4,2,2,1,6,0,2,3,4,2,4,0,4,4,4,0,0,4,2,2,2}, // 21. Trematode Mitochondrial Code
        {4,2,2,2,2,4,3,2,2,6,1,1,2,4,2,6,0,2,3,4,0,0,4,1,2,3}, // 22. Scenedesmus obliquus Mitochondrial Code
        {4,2,2,2,2,4,2,3,2,5,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,4}, // 23. Thraustochytrium Mitochondrial Code
        {4,2,2,2,2,4,2,3,3,6,0,1,2,4,2,4,0,3,4,4,0,0,4,2,2,2}, // 24. Pterobranchia Mitochondrial Code
        {4,2,2,2,2,5,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,2}, // 25. Candidate Division SR1 and Gracilibacteria Code
};


void CodonTable::createCodonTable(unsigned tableId, bool split)
{
    codonTable = new CodonTable(tableId, split);
    codonTable -> setupCodonTable();
}


CodonTable* CodonTable::getInstance()
{
    return codonTable;
}



//--------------------------------------------------//
//--------------------R Wrappers--------------------//
//--------------------------------------------------//



//------------------------------------//
//----------Getter Functions----------//
//------------------------------------//
unsigned CodonTable::getTableIdR()
{
    return getTableId() + 1;
}


std::vector<std::vector<unsigned>> CodonTable::getCodonIndexListingR()
{
    std::vector<std::vector<unsigned>> RV = getCodonIndexListing();

    //Increment every element in the vector since R is 1 indexed.
    for (unsigned i = 0; i < RV.size(); i++)
    {
        for(unsigned j = 0; j < RV[i].size(); j++)
        {
            RV[i][j]++;
        }
    }

    return RV;
}


std::vector<std::vector<unsigned>> CodonTable::getCodonIndexListingWithoutReferenceR()
{
    std::vector<std::vector<unsigned>> RV = getCodonIndexListingWithoutReference();

    //Increment every element in the vector since R is 1 indexed.
    for (unsigned i = 0; i < RV.size(); i++)
    {
        for(unsigned j = 0; j < RV[i].size(); j++)
        {
            RV[i][j]++;
        }
    }

    return RV;
}


std::map <std::string, unsigned> CodonTable::getAAMapR()
{
    //Return the map where all the indices are start from 1.
    std::map <std::string, unsigned>::iterator mit;
    std::map <std::string, unsigned> AAMapCopy = getAAMap();
    for(mit = AAMapCopy.begin(); mit != AAMapCopy.end(); mit++)
    {
        mit -> second++;
    }

    return AAMapCopy;
}


std::map <std::string, unsigned> CodonTable::getForParamVectorMapR()
{
    //Return the map where all the indices are start from 1.
    std::map <std::string, unsigned>::iterator mit;
    std::map <std::string, unsigned> forParamVectorMapCopy = getForParamVectorMap();
    for(mit = forParamVectorMapCopy.begin(); mit != forParamVectorMapCopy.end(); mit++)
    {
        mit -> second++;
    }

    return forParamVectorMapCopy;
}


unsigned CodonTable::getNumCodonsForAAIndexR(unsigned aaIndex, bool forParamVector)
{
    unsigned numCodons = 0;
    bool check = checkIndex(aaIndex, 1, (unsigned)AAListing.size());
    if (check)
    {
        aaIndex--;
        numCodons = getNumCodonsForAAIndex(aaIndex, forParamVector);
    }
    return numCodons;
}


std::string CodonTable::getForParamVectorCodonR(unsigned codonIndex)
{
    std::string codon = "";
    bool check = checkIndex(codonIndex, 1, (unsigned)forParamVectorListing.size());
    if (check)
    {
        codonIndex--;
        codon = getForParamVectorCodon(codonIndex);
    }
    return codon;
}



//--------------------------------------//
//----------Mapping Operations----------//
//--------------------------------------//
unsigned CodonTable::AAToAAIndexR(std::string aa)
{
    unsigned aaIndex = 0;
    aa[0] = (char)std::toupper(aa[0]);
    if (AAMap.find(aa) == AAMap.end())
    {
        std::cerr <<"AA, " << aa <<" is not a valid AA for table " << getTableIdR() <<".\n";
        std::cerr <<"Returning a value of 0.\n";
    }
    else
    {
        aaIndex = AAToAAIndex(aa) + 1; //The plus 1 is for R (1 indexed).
    }

    return aaIndex;
}


std::vector <unsigned> CodonTable::AAIndexToCodonRangeR(unsigned aaIndex, bool forParamVector)
{
    std::vector <unsigned> RV;
    bool check = checkIndex(aaIndex, 1, (unsigned)AAListing.size());
    if (check)
    {
        aaIndex--;
        RV = AAIndexToCodonRange(aaIndex, forParamVector);
        for (unsigned i = 0; i < RV.size(); i++)
        {
            RV[i]++; //R needs 1 indexed answers.
        }
    }

    return RV;
}


std::string CodonTable::indexToCodonR(unsigned index, bool forParamVector)
{
    std::string codon = "";
    bool check;
    check = forParamVector ? checkIndex(index, 1, (unsigned)forParamVectorListing.size()) : checkIndex(index, 1, 64);
    if (check)
    {
        index--;
        codon = indexToCodon(index, forParamVector);
    }
    else
    {
        std::cerr <<"Returning a blank string\n";
    }

    return codon;
}


std::vector <unsigned> CodonTable::AAToCodonRangeR(std::string aa, bool forParamVector)
{
    std::vector <unsigned> RV;
    aa[0] = (char)std::toupper(aa[0]);
    if (AAMap.find(aa) == AAMap.end())
    {
        std::cerr <<"AA, " << aa <<" is not a valid AA for table " << getTableIdR() <<".\n";
        std::cerr <<"Returning an empty vector.\n";
    }
    else
    {
        RV = AAToCodonRange(aa, forParamVector);
        for (unsigned i = 0; i < RV.size(); i++)
        {
            RV[i]++; //R needs 1 indexed answers.
        }
    }

    return RV; 
}


std::vector<std::string> CodonTable::AAToCodonR(std::string aa, bool forParamVector)
{
    std::vector <std::string> RV;
    aa[0] = (char)std::toupper(aa[0]);
    if (AAMap.find(aa) == AAMap.end())
    {
        std::cerr <<"AA, " << aa <<" is not a valid AA for table " << getTableIdR() <<".\n";
        std::cerr <<"Returning an empty vector.\n";
    }
    else
    {
        RV = AAToCodon(aa, forParamVector);
    }

    return RV;
}


std::string CodonTable::codonToAAR(std::string& codon)
{
    std::string AA = "";
    codon[0] = (char) std::toupper(codon[0]);
    codon[1] = (char) std::toupper(codon[1]);
    codon[2] = (char) std::toupper(codon[2]);
    if (codonToIndexWithReference.find(codon) == codonToIndexWithReference.end())
    {
        std::cerr <<"Bad codon, " << codon <<", given. Returning a blank AA\n";
    }
    else
    {
        AA = codonToAA(codon);
    }

    return AA;
}


unsigned CodonTable::codonToIndexR(std::string& codon, bool forParamVector)
{
    unsigned index = 0;
    codon[0] = (char) std::toupper(codon[0]);
    codon[1] = (char) std::toupper(codon[1]);
    codon[2] = (char) std::toupper(codon[2]);

    if ((forParamVector ? forParamVectorMap.find(codon) : codonToIndexWithReference.find(codon)) ==
            (forParamVector ? forParamVectorMap.end() : codonToIndexWithReference.end()))
    {
        std::cerr <<"Error with codon, " << codon <<", returning 0.\n";
        index = 0;
    }
    else
    {
        index = codonToIndex(codon, forParamVector);
        index++; //For R index purposes
    }

    return index;
}


unsigned CodonTable::codonToAAIndexR(std::string& codon)
{
    unsigned AAIndex = 0;
    codon[0] = (char) std::toupper(codon[0]);
    codon[1] = (char) std::toupper(codon[1]);
    codon[2] = (char) std::toupper(codon[2]);

    if (codonToIndexWithReference.find(codon) == codonToIndexWithReference.end())
    {
        std::cerr <<"Error with codon, " << codon <<". Not a valid codon, returning 0.\n";
    }
    else
    {
        AAIndex = codonToAAIndex(codon) + 1;
    }

    return AAIndex;
}

std::string CodonTable::indexToAAR(unsigned aaIndex)
{
    std::string AA = "";
    bool check = checkIndex(aaIndex, 1, (unsigned)AAListing.size());
    if (check)
    {
        aaIndex--;
        AA = indexToAA(aaIndex);
    }

    return AA;
}



//----------------------------------------------//
//---------Static Variables & Functions---------//
//----------------------------------------------//
std::string CodonTable::getSer2R()
{
    return Ser2;
}


std::string CodonTable::getSer1R()
{
    return Ser1;
}


std::string CodonTable::getThr4_1R()
{
    return Thr4_1;
}


std::string CodonTable::getThr4_2R()
{
    return Thr4_2;
}


std::string CodonTable::getLeu1R()
{
    return Leu1;
}


std::vector <std::string> CodonTable::getAminoAcidArrayR()
{
    return aminoAcidArray;
}


std::vector <std::string> CodonTable::getAminoAcidArrayWithoutSplitR()
{
    return aminoAcidArrayWithoutSplit;
}


std::vector <std::vector <unsigned>> CodonTable::getNumCodonsPerAAForTableR()
{
    std::vector <std::vector <unsigned>> RV(25);

    for (unsigned index = 0; index < 25; index++)
    {
        for (unsigned j = 0; j < 26; j++)
        {
            RV[index].push_back(numCodonsPerAAForTable[index][j]);
        }
    }

    return RV;
}
//Hard coded values (25,26) are ok here since this is a
//constant and there are no functions to get the size of an array.


std::vector <std::string> CodonTable::getCodonTableDefinitionR()
{
    return codonTableDefinition;
}


std::vector<std::string> CodonTable::getCodonArrayR()
{
    std::vector<std::string> RV;
    for (unsigned i = 0; i < 64; i++) RV.push_back(CodonTable::codonArray[i]);
    return RV;
}


// ---------------------------------------------------------------------------
// ----------------------------- RCPP MODULE ---------------------------------
// ---------------------------------------------------------------------------
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
RCPP_MODULE(CodonTable_mod)
{
	class_<CodonTable>( "CodonTable" )
		.constructor("Empty constructor - sets to codon table 1 and that amino acids should be split.")
		.constructor<unsigned, bool>("Indicate what codon table to use and if amino acids should be split.")

        //Getter functions:
		.method("getTableId", &CodonTable::getTableIdR)
		.method("getSplitAA", &CodonTable::getSplitAA)
        .method("getCodonIndexListing", &CodonTable::getCodonIndexListingR)
        .method("getCodonIndexListingWithoutReference", &CodonTable::getCodonIndexListingWithoutReferenceR)
        .method("getAAListing", &CodonTable::getAAListing)
        .method("getForParamVectorListing", &CodonTable::getForParamVectorListing)
        .method("getCodonToAAMap", &CodonTable::getCodonToAAMap)
        .method("getAAMap", &CodonTable::getAAMapR)
        .method("getAAToNumCodonsMap", &CodonTable::getAAToNumCodonsMap)
        .method("getForParamVectorMap", &CodonTable::getForParamVectorMapR)

        .method("getNumCodonsForAA", &CodonTable::getNumCodonsForAA)
        .method("getNumCodonsForAAIndex", &CodonTable::getNumCodonsForAAIndexR)
        .method("getForParamVectorCodon", &CodonTable::getForParamVectorCodonR)


        //Mapping operations:
        .method("setupCodonTable", &CodonTable::setupCodonTable)
        .method("AAToAAIndex", &CodonTable::AAToAAIndexR)
        .method("AAIndexToCodonRange", &CodonTable::AAIndexToCodonRangeR)
        .method("AAToCodonRange", &CodonTable::AAToCodonRangeR)
        .method("AAToCodon", &CodonTable::AAToCodonR)
        .method("indexToCodon", &CodonTable::indexToCodonR)
        .method("codonToAA", &CodonTable::codonToAAR)
        .method("codonToIndex", &CodonTable::codonToIndexR)
        .method("codonToAAIndex", &CodonTable::codonToAAIndexR)
        .method("indexToAA", &CodonTable::indexToAAR)
		;

		//Static functions:
		function("getAminoAcidArray", &CodonTable::getAminoAcidArrayR);
        function("getAminoAcidArrayWithoutSplit", &CodonTable::getAminoAcidArrayWithoutSplitR);
		function("getNumCodonsPerAAForTable", &CodonTable::getNumCodonsPerAAForTableR);
		function("getCodonTableDefinition", &CodonTable::getCodonTableDefinitionR);
	    function("getSer2", &CodonTable::getSer2R);
	    function("getSer1", &CodonTable::getSer1R);
	    function("getThr4_1", &CodonTable::getThr4_1R);
	    function("getThr4_2", &CodonTable::getThr4_2R);
	    function("getLeu1", &CodonTable::getLeu1R);
	    function("getCodonArray", &CodonTable::getCodonArrayR);
	    function("getInstance", &CodonTable::getInstance);
	    function("createCodonTable", &CodonTable::createCodonTable);

}
#endif