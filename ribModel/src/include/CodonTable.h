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
        std::vector<std::vector<unsigned>> codonIndexListing; //Stored by AA index then codon index.
        std::vector<std::vector<unsigned>> codonIndexListingWithoutReference;
        //Stored by AA index then codon index. Excludes the last codon index in every AA grouping.

        std::vector <std::string> AAListing; //List of all AAs for the current tableId and split condition.
        std::vector <std::string> forParamVectorListing; //List of all codons without the last codon in every AA group.
        std::map <std::string, std::string> codonToAAMap; //Maps ALL codons for current conditions to AAs.
        std::map <std::string, unsigned> AAMap; //Maps currently used AAs to indices.
        std::map <std::string, unsigned> AAToNumCodonsMap;
        //Maps currently used AAs to the number of codons that code for them.

        std::map <std::string, unsigned> forParamVectorMap;
        //Maps codons to indices (not including last in each AA grouping).


    public:

        //Constructors & destructors:
        explicit CodonTable(); //Defaults to table 1 and splitting AA
        CodonTable(unsigned _tableId, bool _splitAA);
        virtual ~CodonTable();
        CodonTable(const CodonTable& other); //Todo: Need? If so update the function.
        CodonTable& operator=(const CodonTable& other); //Todo: Need? if so update the function.



        //Getter functions:
        unsigned getTableId();
        bool getSplitAA();
        std::vector<std::vector<unsigned>> getCodonIndexListing();
        std::vector<std::vector<unsigned>> getCodonIndexListingWithoutReference();
        std::vector<std::string> getAAListing();
        std::vector <std::string> getForParamVectorListing(); //List of all codons without the last codon in every AA group.
        std::map <std::string, std::string> getCodonToAAMap(); //Maps ALL codons for current conditions to AAs.
        std::map <std::string, unsigned> getAAMap(); //Maps currently used AAs to indices.
        std::map <std::string, unsigned> getAAToNumCodonsMap();
        std::map <std::string, unsigned> getForParamVectorMap();

        unsigned getNumCodonsForAA(std::string aa, bool forParamVector = false);
        unsigned getNumCodonsForAAIndex(unsigned aaIndex, bool forParamVector = false);
        std::string getForParamVectorCodon(unsigned codonIndex);



        //Mapping operations:
        unsigned AAToAAIndex(std::string aa);
        std::vector <unsigned> AAIndexToCodonRange(unsigned aaIndex, bool forParamVector = false);
        std::string indexToCodon(unsigned index, bool forParamVector = false);
        std::vector <unsigned> AAToCodonRange(std::string aa, bool forParamVector = false);
        std::vector<std::string> AAToCodon(std::string aa, bool forParamVector = false);
        std::string codonToAA(std::string& codon);
        unsigned codonToIndex(std::string& codon, bool forParamVector = false);
        unsigned codonToAAIndex(std::string& codon);
        std::string indexToAA(unsigned aaIndex);



        //Other functions:
        void setupCodonTable(); //Sets up the private variables that do all the mappings.
        bool checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound);


        //Static variables & functions:
        static const std::string Ser2;
        static const std::string Ser1; //Necessary for codon table 12
        static const std::string Thr4_1; //Necessary for codon table 3
        static const std::string Thr4_2; //Necessary for codon table 3
        static const std::string Leu1; //Necessary for codon table 16, 22


        static const std::vector <std::string> aminoAcidArray; //Index = AA
        static const std::vector <std::string> aminoAcidArrayWithoutSplit; //Array containing all non-split AAs.
        static const std::map<std::string, unsigned> codonToIndexWithReference; //Map of indices to all codons.
        static const std::string codonArray[]; //List of codons.
        static const std::vector <std::string> codonTableDefinition;
        //Description title for each codon table according to NCBI.

        static const unsigned numCodonsPerAAForTable[25][26]; //Sized on tableId and AA.


        static void createCodonTable(unsigned tableId, bool split = true); //Used to create the singleton instance.
        static CodonTable* getInstance(); //Get the singleton instance.



        //--------------------R WRAPPERS--------------------//
        //Getter functions:
        unsigned getTableIdR();
        std::vector<std::vector<unsigned>> getCodonIndexListingR();
        std::vector<std::vector<unsigned>> getCodonIndexListingWithoutReferenceR();
        std::map <std::string, unsigned> getAAMapR();
        std::map <std::string, unsigned> getForParamVectorMapR();

        unsigned getNumCodonsForAAIndexR(unsigned aaIndex, bool forParamVector = false);
        std::string getForParamVectorCodonR(unsigned codonIndex);



        //Mapping operations:
        unsigned AAToAAIndexR(std::string aa);
        std::vector <unsigned> AAIndexToCodonRangeR(unsigned aaIndex, bool forParamVector = false);
        std::string indexToCodonR(unsigned index, bool forParamVector = false);
        std::vector <unsigned> AAToCodonRangeR(std::string aa, bool forParamVector = false);
        std::vector<std::string> AAToCodonR(std::string aa, bool forParamVector = false);
        std::string codonToAAR(std::string& codon);
        unsigned codonToIndexR(std::string& codon, bool forParamVector = false);
        unsigned codonToAAIndexR(std::string& codon);
        std::string indexToAAR(unsigned aaIndex);



        //Static getter functions:
        static std::string getSer2R();
        static std::string getSer1R();
        static std::string getThr4_1R();
        static std::string getThr4_2R();
        static std::string getLeu1R();

        static std::vector <std::string> getAminoAcidArrayR();
        static std::vector <std::string> getAminoAcidArrayWithoutSplitR();
        static std::vector <std::vector <unsigned>> getNumCodonsPerAAForTableR();
        static std::vector <std::string> getCodonTableDefinitionR();

        static std::vector<std::string> getCodonArrayR();
};

#endif
