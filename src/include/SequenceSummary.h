#ifndef SequenceSummary_H
#define SequenceSummary_H


#include "Utility.h"


#include <string>
#include <map>
#include <algorithm>
#include <cctype>
#include <vector>
#include <array>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

class SequenceSummary
{
	private:

		std::array<unsigned, 64> ncodons; //64 for the number of codons.
		std::array<unsigned, 64> RFPObserved; //64 for the number of codons.
		std::array<unsigned, 22> naa; //22 for the number of amino acids.
		std::vector <std::vector <unsigned>> codonPositions; // used in FONSEModel.
        // index is the codonID, size of 64 for number of codons
        // subindex is the position of each occurrence of the codonID specified.

        // TODO: Probably remove
        std::vector <std::vector <unsigned>> RFP_count;
		//index is the RFP_count for the category specified via index
		//subindex is number of position

	public:

		//Static Member Variables:
		static const std::string Ser2;
		static const std::vector<std::string> AminoAcidArray;
		static const std::string codonArray[];
		static const std::string codonArrayParameter[];
		static const std::map<std::string, unsigned> aaToIndex;
		static const std::map<std::string, unsigned> codonToIndexWithReference;
		static const std::map<std::string, unsigned> codonToIndexWithoutReference;


		//Constructors & Destructors:
		explicit SequenceSummary();
		SequenceSummary(const std::string& sequence);
		SequenceSummary(const SequenceSummary& other);
		SequenceSummary& operator=(const SequenceSummary& other);
		bool operator==(const SequenceSummary& other) const;
		virtual ~SequenceSummary(); //TODO:Why is this virtual????


		//Data Manipulation Functions (All tested):
		unsigned getAACountForAA(std::string aa);
		unsigned getAACountForAA(unsigned aaIndex);
		unsigned getCodonCountForCodon(std::string& codon);
		unsigned getCodonCountForCodon(unsigned codonIndex);
		unsigned getRFPObserved(std::string codon);
		unsigned getRFPObserved(unsigned codonIndex);
		void setRFPObserved(unsigned codonIndex, unsigned value);
		std::vector <unsigned> *getCodonPositions(std::string codon);
		std::vector <unsigned> *getCodonPositions(unsigned index);


		//RFP Functions (for PA and PANSE models) (All tested):
		std::vector <unsigned> getRFP_count(unsigned categoryIndex);
		void initRFP_count(unsigned numCategories);
		void setRFP_count(unsigned categoryIndex, std::vector <unsigned> arg);


		//Other Functions (All tested):
		void clear();
		bool processSequence(const std::string& sequence);
        bool processPA(std::vector <std::vector <unsigned>> table);


		//Static Functions:
		static unsigned AAToAAIndex(std::string aa); //Moving to CT
		static void AAIndexToCodonRange(unsigned aaIndex, unsigned& start, unsigned& end, bool forParamVector = false); //Moving to CT
		static void AAToCodonRange(std::string aa, unsigned& start, unsigned& end, bool forParamVector = false); //Moving to CT
		static std::vector<std::string> AAToCodon(std::string aa, bool forParamVector = false); //Moving to CT, but used in R currently
		static std::string codonToAA(std::string& codon); //Moving to CT
		static unsigned codonToIndex(std::string& codon, bool forParamVector = false); //Moving to CT
		static unsigned codonToAAIndex(std::string& codon); //Moving to CT
		static std::string indexToAA(unsigned aaIndex); //Moving to CT
		static std::string indexToCodon(unsigned index, bool forParamVector = false); //Moving to CT
		static unsigned GetNumCodonsForAA(std::string& aa, bool forParamVector = false); //Moving to CT
		static char complimentNucleotide(char ch); //Tested
		static std::vector<std::string> aminoAcids(); //Moving to CT, but used in R currently
		static std::vector<std::string> codons(); //Moving to CT, but used in R currently


	protected:
};

#endif // SequenceSummary_H


