#ifndef CodonTable_H
#define CodonTable_H

#include <string>
#include <map>
#include <algorithm>
#include <cctype>
#include <vector>
#include <array>

class CodonTable
{
    private:

        static CodonTable *codonTable;
        unsigned tableId;
        bool splitAA;
        std::vector<std::vector<unsigned>> codon_mapping; // dim: AA, codon
        std::vector<std::vector<unsigned>> codon_mapping_without_reference; //AA, codon
        std::vector<std::string> aa_mapping;
        std::map <std::string, std::string> codonToAAMap;
        std::map <std::string, unsigned> AAMap;
        std::map <std::string, unsigned> AAToNumCodonsMap;
        std::map <std::string, unsigned> forParamVectorMap;



    public:
        static const std::string Ser2;
        static const std::string Ser1; // necessary for codon table 12
        static const std::string Thr4_1; // necessary for codon table 3
        static const std::string Thr4_2; // necessary for codon table 3
        static const std::string Leu1; // necessary for codon table 16, 22

		static const std::string AminoAcidArray[26]; // dim:AA
        static const std::string AminoAcidArrayWithoutSplit[21];
		static const unsigned numCodonsPerAAForTable[25][26]; // dim: table, AA
		static const std::string codonTableDefinition[25];
		static const std::string codonArray[];
		static const std::string codonArrayParameter[];
		static const std::map<std::string, unsigned> aaToIndex;
        static const std::map<std::string, unsigned> aaToIndexWithoutSplit;
		static const std::map<std::string, unsigned> codonToIndexWithReference;
		static const std::map<std::string, unsigned> codonToIndexWithoutReference;
        static void createCodonTable(unsigned tableId, bool split = true);
        static CodonTable* getInstance();


        //Constructors & destructors:
        explicit CodonTable();
        CodonTable(unsigned _tableId, bool _splitAA);
        virtual ~CodonTable();
        CodonTable(const CodonTable& other);
        CodonTable& operator=(const CodonTable& other);

        void setupCodonTable();
        unsigned AAToAAIndex(std::string aa);
        std::vector <unsigned> AAIndexToCodonRange(unsigned aaIndex, bool forParamVector = false);
        std::vector <unsigned> AAToCodonRange(std::string aa, bool forParamVector = false);
        std::vector<std::string> AAToCodon(std::string aa, bool forParamVector = false);
        std::string indexToCodon(unsigned index);
        std::string codonToAA(std::string& codon);
        unsigned codonToIndex(std::string& codon, bool forParamVector = false);
        unsigned codonToAAIndex(std::string& codon);
        std::string indexToAA(unsigned aaIndex);
        unsigned getNumCodons(std::string aa, bool forParamVector = false);
        unsigned getNumCodons(unsigned aaIndex, bool forParamVector = false);

};

#endif
