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

	// TODO CAREFULL! tableId is 1 indexed, but we need it 0 indexed! ALSO NOT ALL IDs ARE USED, SEE COMMENT FOR CodonTable::codonTableDefinition[]
    //ctor
	switch (tableId)
	{
	case 1:
		break;

	}
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

const std::string CodonTable::AminoAcidArray[25][23] = {
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", CodonTable::Ser2, "X", ""}, // 1. The Standard Code
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", CodonTable::Ser2, "X", ""}, // 2. The Vertebrate Mitochondrial Code
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", CodonTable::Thr4_1, CodonTable::Thr4_2, "V", "W", "Y", CodonTable::Ser2, "X"}, // 3. The Yeast Mitochondrial Code
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", CodonTable::Ser2, "X", ""}, // 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", CodonTable::Ser2, "X", ""}, // 5. The Invertebrate Mitochondrial Code
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", CodonTable::Ser2, "X", ""}, // 6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
	{"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""}, // 7. Invalid Codon Table
	{"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""}, // 8. Invalid Codon Table
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", CodonTable::Ser2, "X", ""}, // 9. The Echinoderm and Flatworm Mitochondrial Code
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", CodonTable::Ser2, "X", ""}, // 10. The Euplotid Nuclear Code
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", CodonTable::Ser2, "X", ""}, // 11. The Bacterial, Archaeal and Plant Plastid Code
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", CodonTable::Ser2, CodonTable::Ser1, "X"}, // 12. The Alternative Yeast Nuclear Code
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", CodonTable::Ser2, "X", ""}, // 13. The Ascidian Mitochondrial Code
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", CodonTable::Ser2, "X", ""}, // 14. The Alternative Flatworm Mitochondrial Code
	{"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""}, // 15. Invalid Codon Table
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", CodonTable::Leu1, "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", CodonTable::Ser2, "X"}, // 16. Chlorophycean Mitochondrial Code
	{"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""}, // 17. Invalid Codon Table
	{"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""}, // 18. Invalid Codon Table
	{"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""}, // 19. Invalid Codon Table
	{"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""}, // 20. Invalid Codon Table
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", CodonTable::Ser2, "X", ""}, // 21. Trematode Mitochondrial Code
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", CodonTable::Leu1, "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", CodonTable::Ser2, "X"}, // 22. Scenedesmus obliquus Mitochondrial Code
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", CodonTable::Ser2, "X", ""}, // 23. Thraustochytrium Mitochondrial Code
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", CodonTable::Ser2, "X", ""}, // 24. Pterobranchia Mitochondrial Code
	{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", CodonTable::Ser2, "X", ""}  // 25. Candidate Division SR1 and Gracilibacteria Code
};

const std::map<std::string, unsigned> CodonTable::aaToIndex = {{"A", 0}, {"C", 1}, {"D", 2}, {"E", 3}, {"F", 4},
	{"G", 5}, {"H", 6}, {"I", 7}, {"K", 8}, {"L", 9}, {"M", 10}, {"N", 11}, {"P", 12}, {"Q", 13}, {"R", 14}, {"S", 15},
	{"T", 16}, {"V", 17}, {"W", 18}, {"Y", 19}, {CodonTable::Ser2, 20}, {"X", 21}};

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
const unsigned CodonTable::numCodonsPerAAForTable[25][23] = {
		{4,2,2,2,2,4,2,3,2,6,1,2,4,2,6,4,4,4,1,2,2,3,0}, // 1. The Standard Code
		{4,2,2,2,2,4,2,2,2,6,2,2,4,2,4,4,4,4,2,2,2,4,0}, // 2. The Vertebrate Mitochondrial Code
		{4,2,2,2,2,4,2,2,2,2,2,2,4,2,4,4,4,4,4,2,2,2,2}, // 3. The Yeast Mitochondrial Code
		{4,2,2,2,2,4,2,3,2,6,1,2,4,2,6,4,4,4,2,2,2,2,0}, // 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
		{4,2,2,2,2,4,2,2,2,6,2,2,4,2,4,4,4,4,2,2,4,2,0}, // 5. The Invertebrate Mitochondrial Code
		{4,2,2,2,2,4,2,3,2,6,1,2,4,4,6,4,4,4,1,2,2,1,0}, // 6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
		{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 7. Invalid Codon Table
		{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 8. Invalid Codon Table
		{4,2,2,2,2,4,2,3,1,6,1,3,4,2,4,4,4,4,2,2,4,2,0}, // 9. The Echinoderm and Flatworm Mitochondrial Code
		{4,3,2,2,2,4,2,3,2,6,1,2,4,2,6,4,4,4,1,2,2,2,0}, // 10. The Euplotid Nuclear Code
		{4,2,2,2,2,4,2,3,2,6,1,2,4,2,6,4,4,4,1,2,2,3,0}, // 11. The Bacterial, Archaeal and Plant Plastid Codee
		{4,2,2,2,2,4,2,3,2,5,1,2,4,2,6,4,4,4,1,2,2,1,3}, // 12. The Alternative Yeast Nuclear Code
		{4,2,2,2,2,6,2,2,2,6,2,2,4,2,4,4,4,4,2,2,2,2,0}, // 13. The Ascidian Mitochondrial Code
		{4,2,2,2,2,4,2,3,1,6,1,3,4,2,4,4,4,4,2,3,4,1,0}, // 14. The Alternative Flatworm Mitochondrial Code
		{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 15. Invalid Codon Table
		{4,2,2,2,2,4,2,3,2,6,1,1,2,4,2,6,4,4,4,1,2,2,2}, // 16. Chlorophycean Mitochondrial Code
		{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 17. Invalid Codon Table
		{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 18. Invalid Codon Table
		{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 19. Invalid Codon Table
		{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 20. Invalid Codon Table
		{4,2,2,2,2,4,2,2,1,6,2,3,4,2,4,4,4,4,2,2,4,2,0}, // 21. Trematode Mitochondrial Code
		{4,2,2,2,2,4,2,3,2,6,1,1,2,4,2,6,3,4,4,1,2,2,3}, // 22. Scenedesmus obliquus Mitochondrial Code
		{4,2,2,2,2,4,2,3,2,5,1,2,4,2,6,4,4,4,1,2,2,4,0}, // 23. Thraustochytrium Mitochondrial Code
		{4,2,2,2,2,4,2,3,3,6,1,2,4,2,4,4,4,4,2,2,3,2,0}, // 24. Pterobranchia Mitochondrial Code
		{4,2,2,2,2,5,2,3,2,6,1,2,4,2,6,4,4,4,1,2,2,2,0}, // 25. Candidate Division SR1 and Gracilibacteria Code
};

// --- CODON TABLE SPECIFIC MAPPER FUNCTIONS ------

unsigned CodonTable::AAToAAIndex(std::string aa)
{
	return CodonTable::aaToIndex.find(aa) -> second;
}

unsigned CodonTable::getNumCodons(std::string aa)
{
	return getNumCodons(AAToAAIndex(aa));
}

unsigned CodonTable::getNumCodons(unsigned aa)
{
	return numCodonsPerAAForTable[tableId][aa];
}

void CodonTable::setupCodonTable()
{
	unsigned numAA = 21;
	if(tableId == 1 && splitAA == false) numAA = 20;
	// TODO add sizes for other codon tables

	codon_mapping.resize(numAA);

	for(unsigned aaIndex = 0; aaIndex < numAA; aaIndex++)
	{
		codon_mapping[aaIndex].resize(getNumCodons(aaIndex));
		// TODO fill with codon index
	}
}



