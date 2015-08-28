#include "include/CodonTable.h"

#include <iostream>

//---- OPERATOR AND CONSTRUCTOR ------
CodonTable::CodonTable()
{
    tableId = 1; //standard codon table by NCBI
    splitAA = true;
}


CodonTable::CodonTable(unsigned _tableId, bool _splitAA) : tableId(_tableId), splitAA(_splitAA)
{
	if(tableId == 7 || tableId == 8 || tableId == 15 || tableId == 17 || tableId == 18 || tableId == 19 || tableId == 20)
	{
		tableId = 1; //standard codon table by NCBI
        std::cerr << "Invalid codon table: " << tableId << " using default codon table (NCBI codon table 1)\n";
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

    // STATIC MEMBERS

    const std::string CodonTable::Ser2 = "Z";
    const std::string CodonTable::Ser1 = "J";
    const std::string CodonTable::Thr4_1 = "O";
    const std::string CodonTable::Thr4_2 = "B";
    const std::string CodonTable::Leu1 = "U";

    const std::string CodonTable::AminoAcidArray[] = {
        "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", CodonTable::Leu1, "M", "N", "P", "Q", "R",
        CodonTable::Ser1, CodonTable::Ser2, "S", "T", CodonTable::Thr4_1, CodonTable::Thr4_2, "V", "W", "Y", "X"};
    const std::string CodonTable::AminoAcidArrayWithoutSplit[] = {
        "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X"};


    const std::map<std::string, unsigned> CodonTable::aaToIndex = {{"A", 0}, {"C", 1}, {"D", 2}, {"E", 3}, {"F", 4},
     {"G", 5}, {"H", 6}, {"I", 7}, {"K", 8}, {"L",9}, {CodonTable::Leu1, 10}, {"M", 11}, {"N", 12}, {"P", 13}, {"Q", 14}, {"R", 15},
                                                                   {CodonTable::Ser1, 16}, {CodonTable::Ser2, 17}, {"S", 18},
     {"T", 19}, {CodonTable::Thr4_1, 20}, {CodonTable::Thr4_2, 21}, {"V", 22}, {"W", 23}, {"Y", 24}, {"X", 25}};

const std::map<std::string, unsigned> CodonTable::aaToIndexWithoutSplit = {{"A", 0}, {"C", 1}, {"D", 2}, {"E", 3}, {"F", 4},
  {"G", 5}, {"H", 6}, {"I", 7}, {"K", 8}, {"L",9}, {"M", 10}, {"N", 11}, {"P", 12}, {"Q", 13}, {"R", 14}, {"S", 15}, {"T", 16}, {"V", 17}, {"W", 18}, {"Y", 19}, {"X", 20}};

    const std::map<std::string, unsigned> CodonTable::codonToIndexWithReference = {{"GCA", 0}, {"GCC", 1}, {"GCG", 2},
        {"GCT", 3}, {"TGC", 4}, {"TGT", 5}, {"GAC", 6}, {"GAT", 7}, {"GAA", 8}, {"GAG", 9}, {"TTC", 10}, {"TTT", 11},
        {"GGA", 12}, {"GGC", 13}, {"GGG", 14}, {"GGT", 15}, {"CAC", 16}, {"CAT", 17}, {"ATA", 18}, {"ATC", 19}, {"ATT", 20},
        {"AAA", 21}, {"AAG", 22}, {"CTA", 23}, {"CTC", 24}, {"CTG", 25}, {"CTT", 26}, {"TTA", 27}, {"TTG", 28}, {"ATG", 29},
        {"AAC", 30}, {"AAT", 31}, {"CCA", 32}, {"CCC", 33}, {"CCG", 34}, {"CCT", 35}, {"CAA", 36}, {"CAG", 37}, {"AGA", 38},
        {"AGG", 39}, {"CGA", 40}, {"CGC", 41}, {"CGG", 42}, {"CGT", 43}, {"TCA", 44}, {"TCC", 45}, {"TCG", 46}, {"TCT", 47},
        {"ACA", 48}, {"ACC", 49}, {"ACG", 50}, {"ACT", 51}, {"GTA", 52}, {"GTC", 53}, {"GTG", 54}, {"GTT", 55}, {"TGG", 56},
        {"TAC", 57}, {"TAT", 58}, {"AGC", 59}, {"AGT", 60}, {"TAA", 61}, {"TAG", 62}, {"TGA", 63}};

    const std::map<std::string, unsigned> CodonTable::codonToIndexWithoutReference = {{"GCA", 0}, {"GCC", 1},
        {"GCG", 2}, {"TGC", 3}, {"GAC", 4}, {"GAA", 5}, {"TTC", 6}, {"GGA", 7}, {"GGC", 8}, {"GGG", 9}, {"CAC", 10},
        {"ATA", 11}, {"ATC", 12}, {"AAA", 13}, {"CTA", 14}, {"CTC", 15}, {"CTG", 16}, {"CTT", 17}, {"TTA", 18}, {"AAC", 19},
        {"CCA", 20}, {"CCC", 21}, {"CCG", 22}, {"CAA", 23}, {"AGA", 24}, {"AGG", 25}, {"CGA", 26}, {"CGC", 27}, {"CGG", 28},
        {"TCA", 29}, {"TCC", 30}, {"TCG", 31}, {"ACA", 32}, {"ACC", 33}, {"ACG", 34}, {"GTA", 35}, {"GTC", 36}, {"GTG", 37},
        {"TAC", 38}, {"AGC", 39}};

    const std::string CodonTable::codonArray[] =
        {"GCA", "GCC", "GCG", "GCT", "TGC", "TGT", "GAC", "GAT", "GAA", "GAG",
        "TTC", "TTT", "GGA", "GGC", "GGG", "GGT", "CAC", "CAT", "ATA", "ATC",
        "ATT", "AAA", "AAG", "CTA", "CTC", "CTG", "CTT", "TTA", "TTG", "ATG",
        "AAC", "AAT", "CCA", "CCC", "CCG", "CCT", "CAA", "CAG", "AGA", "AGG",
        "CGA", "CGC", "CGG", "CGT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC",
        "ACG", "ACT", "GTA", "GTC", "GTG", "GTT", "TGG", "TAC", "TAT", "AGC",
        "AGT", "TAA", "TAG", "TGA"};

    const std::string CodonTable::codonArrayParameter[] =
        {"GCA", "GCC", "GCG", "TGC", "GAC",
        "GAA", "TTC", "GGA", "GGC", "GGG",
        "CAC", "ATA", "ATC", "AAA", "CTA",
        "CTC", "CTG", "CTT", "TTA", "AAC",
        "CCA", "CCC", "CCG", "CAA", "AGA",
        "AGG", "CGA", "CGC", "CGG", "TCA",
        "TCC", "TCG", "ACA", "ACC", "ACG",
        "GTA", "GTC", "GTG", "TAC", "AGC"};

    // TODO NOTE: THERE IS NO CODON TABLE 7, 8, 15, 17, 18, 19, 20 ACCORDING TO NCBI !
    // http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    const std::string CodonTable::codonTableDefinition[] = {"1. The Standard Code", "2. The Vertebrate Mitochondrial Code",
        "3. The Yeast Mitochondrial Code", "4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code",
        "5. The Invertebrate Mitochondrial Code", "6. The Ciliate, Dasycladacean and Hexamita Nuclear Code", "7. Invalid Codon Table", "8. Invalid Codon Table",
        "9. The Echinoderm and Flatworm Mitochondrial Code", "10. The Euplotid Nuclear Code",
        "11. The Bacterial, Archaeal and Plant Plastid Code", "12. The Alternative Yeast Nuclear Code", "13. The Ascidian Mitochondrial Code",
        "14. The Alternative Flatworm Mitochondrial Code", "15. Invalid Codon Table", "16. Chlorophycean Mitochondrial Code",
        "17. Invalid Codon Table", "18. Invalid Codon Table", "19. Invalid Codon Table", "20. Invalid Codon Table",
        "21. Trematode Mitochondrial Code", "22. Scenedesmus obliquus Mitochondrial Code", "23. Thraustochytrium Mitochondrial Code",
        "24. Pterobranchia Mitochondrial Code",	"25. Candidate Division SR1 and Gracilibacteria Code"};


    // {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "Z", "X"};
    const unsigned CodonTable::numCodonsPerAAForTable[25][26] = {
            {4,2,2,2,2,4,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,3}, // 1. The Standard Code
            {4,2,2,2,2,4,2,2,2,6,0,2,2,4,2,4,0,2,4,4,0,0,4,2,2,4}, // 2. The Vertebrate Mitochondrial Code
            {4,2,2,2,2,4,2,2,2,2,0,2,2,4,2,4,0,2,4,0,4,4,4,2,2,2}, // 3. The Yeast Mitochondrial Code
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
            {4,2,2,2,2,4,2,3,2,6,1,1,2,4,2,6,0,2,3,4,0,0,4,1,2,3}, // 22. Scenedesmus obliquus Mitochondrial Code
            {4,2,2,2,2,4,2,3,2,5,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,4}, // 23. Thraustochytrium Mitochondrial Code
            {4,2,2,2,2,4,2,3,3,6,0,1,2,4,2,4,0,3,4,4,0,0,4,2,2,2}, // 24. Pterobranchia Mitochondrial Code
            {4,2,2,2,2,5,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,2}, // 25. Candidate Division SR1 and Gracilibacteria Code
    };

    // --- CODON TABLE SPECIFIC MAPPER FUNCTIONS ------

unsigned CodonTable::AAToAAIndex(std::string aa)
{
        return CodonTable::AAMap.find(aa) -> second;
}

std::array<unsigned, 8> CodonTable::AAIndexToCodonRange(unsigned aaIndex, bool forParamVector)
{
    std::array<unsigned, 8> aaRange;
    aaRange.fill(100);
    if (forParamVector)
    {
        unsigned numCodons = getNumCodons(aaIndex) - 1;
        for (unsigned i = 0; i < numCodons; i++)
        {
            //TODO: Assuming the index corrosponds to the current mapping
            aaRange[i] = codon_mapping[aaIndex][i];
        }
        //use without reference mapping
    }
    else
    {
        unsigned numCodons = getNumCodons(aaIndex);
        for (unsigned i = 0; i < numCodons; i++)
        {
            //TODO: Assuming the index corrosponds to the current mapping
            aaRange[i] = codon_mapping[aaIndex][i];
        }
    }

    return aaRange;
}

std::array<unsigned, 8> CodonTable::AAToCodonRange(std::string aa, bool forParamVector)
{
    unsigned index = 0;
    for (unsigned i = 0; i < aa_mapping.size(); i++)
    {
        if (aa_mapping[i] == aa)
        {
            index = i;
            break;
        }
    }
    //std::array<unsigned, 8> aaRange;
    return AAIndexToCodonRange(index, forParamVector);
}


std::vector<std::string> CodonTable::AAToCodon(std::string aa, bool forParamVector)
{
    std::vector <std::string> codons;
    std::array <unsigned , 8> aaRange = AAToCodonRange(aa, forParamVector);

        for (unsigned i = 0; i < 8; i++)
        {
            if (aaRange[i] != 100)
            {
                unsigned index = codon_mapping_without_reference[AAToAAIndex(aa)][aaRange[i]];
                codons.push_back(indexToCodon(index));
            }
        }

    return codons;
}


std::string CodonTable::indexToCodon(unsigned index)
{
    return codonArray[index];
}


std::string CodonTable::codonToAA(std::string& codon)
{
    std::map <std::string, std::string>::iterator mit;

    mit = codonToAAMap.find(codon);
    return mit -> second;
}


unsigned CodonTable::codonToIndex(std::string& codon)
{
    return codonToIndexWithReference.find(codon) -> second;
}

unsigned CodonTable::codonToAAIndex(std::string& codon)
{
    std::string AA = codonToAA(codon);
    return AAMap.find(AA) -> second;
}

std::string CodonTable::indexToAA(unsigned aaIndex)
{
    return aa_mapping[aaIndex];
}

unsigned CodonTable::getNumCodons(std::string aa, bool forParamVector)
{
    unsigned numCodons = AAToNumCodonsMap.find(aa);
    return forParamVector ? numCodons - 1 : numCodons;
}

unsigned CodonTable::getNumCodons(unsigned aaIndex, bool forParamVector)
{
    std::string aa = indexToAA(aaIndex);
    return getNumCodons(aa, forParamVector);
}

void CodonTable::setupCodonTable()
{
	unsigned numAA = 21;
	if(tableId >= 0 && tableId <= 5  && !splitAA) numAA = 20; //If not splitting AAs, tables 1-6 all have 20 AA, excluding stop codes
    else if (tableId >= 8 && tableId <= 13  && !splitAA) numAA = 20;
    else if (tableId == 15 && !splitAA) numAA = 20;
    else if (tableId >= 20 && tableId <= 24 && !splitAA) numAA = 20;


	codon_mapping.resize(numAA);

    unsigned aaIndex = 0;
    unsigned filled = 0;
    while (filled != numAA)
    {
        unsigned numCodons = getNumCodons(aaIndex);
        if (numCodons != 0)
        {
            if (aaIndex == 0) //A
            {
                for (unsigned i = 0; i < 4; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
                AAToNumCodonsMap.insert(std::make_pair("A", 4));
            }
            else if (aaIndex == 1) //C
            {
                AAToNumCodonsMap.insert(std::make_pair("C", 2));
                for (unsigned i = 4; i < 6; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
                if (tableId == 10)
                {
                    codon_mapping[filled].push_back(63);
                    AAToNumCodonsMap["C"]++;
                }

            }
            else if (aaIndex == 2) //D
            {
                AAToNumCodonsMap.insert(std::make_pair("D", 2));
                for (unsigned i = 6; i < 8; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
            }
            else if (aaIndex == 3) //E
            {
                AAToNumCodonsMap.insert(std::make_pair("E", 2));
                for (unsigned i = 8; i < 10; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
            }
            else if (aaIndex == 4) //F
            {
                AAToNumCodonsMap.insert(std::make_pair("F", 2));
                for (unsigned i = 10; i < 12; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
            }
            else if (aaIndex == 5) //G
            {
                AAToNumCodonsMap.insert(std::make_pair("G", 4));
                for (unsigned i = 12; i < 16; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
                if (tableId == 12)
                {
                    AAToNumCodonsMap["G"] += 2;
                    codon_mapping[filled].push_back(38);
                    codon_mapping[filled].push_back(39);
                }
                else if (tableId == 24)
                {
                    AAToNumCodonsMap["G"]++;
                    codon_mapping[filled].push_back(63);
                }
            }
            else if (aaIndex == 6) //H
            {
                AAToNumCodonsMap.insert(std::make_pair("H", 2));
                for (unsigned i = 16; i < 18; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
            }
            else if (aaIndex == 7) //I
            {
                AAToNumCodonsMap.insert(std::make_pair("I", 2));
                if (tableId == 0 || tableId == 3 || tableId == 5 || tableId == 8 || tableId == 9 || tableId == 10
                        || tableId == 11 || tableId == 13 || tableId == 15 || tableId >= 21)
                {
                    AAToNumCodonsMap["I"]++;
                    codon_mapping[filled].push_back(18);
                }
                for (unsigned i = 19; i < 21; i++)
                {
                    codon_mapping[filled].push_back(i);
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
                        codon_mapping[filled].push_back(i);
                    }
                }
                else if (tableId == 8 || tableId == 13 || tableId == 20)
                {
                    AAToNumCodonsMap.insert(std::make_pair("K", 1));
                    codon_mapping[filled].push_back(22);
                }

                if (tableId == 23)
                {
                    AAToNumCodonsMap["K"]++;
                    codon_mapping[filled].push_back(39);
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
                            codon_mapping[filled].push_back(i);
                        }
                        else
                        {
                            codon_mapping[filled].push_back(i);
                        }
                    }
                    if ((tableId == 15 || tableId == 21) && !splitAA)
                    {
                        AAToNumCodonsMap["L"]++;
                        codon_mapping[filled].push_back(62);
                    }
                }
                else if (tableId == 2)
                {
                    AAToNumCodonsMap.insert(std::make_pair("L", 2));
                    for (unsigned i = 27; i < 29; i++)
                    {
                        codon_mapping[filled].push_back(i);
                    }
                }
                else if (tableId == 22)
                {
                    AAToNumCodonsMap.insert(std::make_pair("L", 5));
                    for (unsigned i = 23; i < 27; i++)
                    {
                        codon_mapping[filled].push_back(i);
                    }
                    codon_mapping[filled].push_back(28);
                }
            }
            else if (aaIndex == 10) //Leu1
            {
                if (splitAA)
                {
                    AAToNumCodonsMap.insert(std::make_pair(CodonTable::Leu1, 1));
                    codon_mapping[filled].push_back(62);
                }
                else
                {
                    filled--;
                }
            }
            else if (aaIndex == 11) //M
            {
                AAToNumCodonsMap.insert(std::make_pair("M", 1));
                codon_mapping[filled].push_back(29);
                if (tableId == 1 || tableId == 2 || tableId == 4 || tableId == 12 || tableId == 20)
                {
                    AAToNumCodonsMap["M"]++;
                    codon_mapping[filled].push_back(18);
                }
            }
            else if (aaIndex == 12) //N
            {
                AAToNumCodonsMap["M"] += 2;
                for (unsigned i = 30; i < 32; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
                if (tableId == 8 || tableId == 13 || tableId == 20)
                {
                    AAToNumCodonsMap["M"]++;
                    codon_mapping[filled].push_back(21);
                }
            }
            else if (aaIndex == 13) //P
            {
                AAToNumCodonsMap.insert(std::make_pair("P", 4));
                for (unsigned i = 32; i < 36; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
            }
            else if (aaIndex == 14) //Q
            {
                AAToNumCodonsMap.insert(std::make_pair("Q", 2));
                for (unsigned i = 36; i < 38; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
                if (tableId == 5)
                {
                    AAToNumCodonsMap["Q"] += 2;
                    codon_mapping[filled].push_back(61);
                    codon_mapping[filled].push_back(62);
                }
            }
            else if (aaIndex == 15) //R
            {
                AAToNumCodonsMap.insert(std::make_pair("R", 0));
                if (tableId == 0 || tableId == 2 || tableId == 3 || tableId == 5 || (tableId >= 9 && tableId <= 11) || tableId == 15 ||
                        tableId == 21 || tableId == 22 || tableId == 24)
                {
                    AAToNumCodonsMap["R"] += 2;
                    codon_mapping[filled].push_back(38);
                    codon_mapping[filled].push_back(39);
                }

                if (tableId == 2)
                {
                    AAToNumCodonsMap["R"] += 2;
                    codon_mapping[filled].push_back(42);
                    codon_mapping[filled].push_back(43);
                }
                else
                {
                    AAToNumCodonsMap["R"] += 4;
                    for (unsigned i = 40; i < 44; i++)
                    {
                        codon_mapping[filled].push_back(i);
                    }
                }
            }
            else if (aaIndex == 16) //Ser1
            {
                if(splitAA)
                {
                    AAToNumCodonsMap.insert(std::make_pair(CodonTable::Ser1, 1));
                    codon_mapping[filled].push_back(25); //table 12
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
                    codon_mapping[filled].push_back(59);
                    codon_mapping[filled].push_back(60);

                    if (tableId == 4 || tableId == 8 || tableId == 13 || tableId == 20)
                    {
                        AAToNumCodonsMap[CodonTable::Ser2] += 2;
                        codon_mapping[filled].push_back(38);
                        codon_mapping[filled].push_back(39);
                    }
                    else if (tableId == 23)
                    {
                        AAToNumCodonsMap[CodonTable::Ser2]++;
                        codon_mapping[filled].push_back(38);
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
                        codon_mapping[filled].push_back(i);
                    }
                }
                else
                {
                    for(unsigned i = 45; i < 48; i++)
                    {
                        codon_mapping[filled].push_back(i);
                    }
                }

                if (!splitAA)
                {
                    AAToNumCodonsMap["S"] += 2;
                    if (tableId == 11)
                    {
                        AAToNumCodonsMap["S"]++;
                        codon_mapping[filled].push_back(25);
                    }

                    codon_mapping[filled].push_back(59);
                    codon_mapping[filled].push_back(60);

                    if (tableId == 4 || tableId == 8 || tableId == 13 || tableId == 20)
                    {
                        AAToNumCodonsMap["S"] += 2;
                        codon_mapping[filled].push_back(38);
                        codon_mapping[filled].push_back(39);
                    }
                    else if (tableId == 23)
                    {
                        AAToNumCodonsMap["S"]++;
                        codon_mapping[filled].push_back(38);
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
                        codon_mapping[filled].push_back(i);
                    }
                }
                else if (tableId == 2 && !splitAA)
                {
                    AAToNumCodonsMap["T"] += 8;
                    for (unsigned i = 48; i < 52; i++)
                    {
                        codon_mapping[filled].push_back(i);
                    }

                    for (unsigned i = 23; i < 27; i++)
                    {
                        codon_mapping[filled].push_back(i);
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
                        codon_mapping[filled].push_back(i);
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
                        codon_mapping[filled].push_back(i);
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
                    codon_mapping[filled].push_back(i);
                }
            }
            else if (aaIndex == 23) //W
            {
                AAToNumCodonsMap.insert(std::make_pair("W", 1));
                codon_mapping[filled].push_back(56);
                if ((tableId >= 1 && tableId <= 4) || (tableId == 8) || (tableId == 12) || (tableId == 13) || (tableId == 20))
                {
                    AAToNumCodonsMap["W"]++;
                    codon_mapping[filled].push_back(63);
                }
            }
            else if (aaIndex == 24) //Y
            {
                AAToNumCodonsMap.insert(std::make_pair("Y", 2));
                for (unsigned i = 56; i < 58; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
                if (tableId == 13){
                    AAToNumCodonsMap["Y"]++;
                    codon_mapping[filled].push_back(61);
                }
            }
            //ignoring stop amino acid for now
            filled++;
        }
        aaIndex++;
    }

    //fill without refernce mapping
    codon_mapping_without_reference.resize(numAA);
    for (unsigned i = 0; i < numAA; i++)
    {
        for (unsigned j = 0; j < codon_mapping[i].size() - 1; j++)
        {
            codon_mapping_without_reference[i].push_back(codon_mapping[i][j]);
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
                aa_mapping.push_back(AminoAcidArray[i]);
                AAMap.insert(std::make_pair(AminoAcidArray[i], index));
                index++;
            }
        }
    }
    else
    {
        for (unsigned i = 0; i < 26; i++)
        {
            aa_mapping.push_back(AminoAcidArrayWithoutSplit[i]);
            if (numCodonsPerAAForTable[tableId][i] != 0 || AminoAcidArray[i] != CodonTable::Leu1 || AminoAcidArray[i] != CodonTable::Ser1
                    || AminoAcidArray[i] != Ser2 || AminoAcidArray[i] != CodonTable::Thr4_2)
            {
                if (AminoAcidArray[i] == CodonTable::Thr4_1) //this accounts for the table where Thr4_1 and Thr4_2 and filled, but T is the one actually used when not split
                {
                    aa_mapping.push_back("T");
                    AAMap.insert(std::make_pair("T", index));
                }
                aa_mapping.push_back(AminoAcidArray[i]);
                AAMap.insert(std::make_pair(AminoAcidArray[i], index));
                index++;
            }
        }
    }




    //Build a map for the Codon To AA
    for(unsigned AAIndex = 0; AAIndex < codon_mapping.size(); AAIndex++)
    {
        std::string AA = aa_mapping[AAIndex];
        for (unsigned codonIndex = 0; codonIndex < codon_mapping[AAIndex].size(); codonIndex++)
        {
            std::string codon = codonArray[codon_mapping[AAIndex][codonIndex]];
            codonToAAMap.insert(std::make_pair(codon, AA));
        }
    }


}



