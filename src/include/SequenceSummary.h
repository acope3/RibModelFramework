#ifndef SequenceSummary_H
#define SequenceSummary_H


#include <string>
#include <map>
#include <algorithm>
#include <cctype>
#include <vector>
#include <array>
#include <iostream>
#include "Utility.h"

#ifndef STANDALONE
#include <Rcpp.h>
#endif
#include <array>

class SequenceSummary
{
	private:

		std::array<unsigned, 64> ncodons; //64 for the number of codons.
		std::array<unsigned, 64> RFPObserved; //64 for the number of codons.
		std::array<unsigned, 22> naa; //22 for the number of amino acids.
		std::vector <std::vector <unsigned>> codonPositions;
		std::vector <unsigned> RFP_count; //index is number of position

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


		//Data Manipulation Functions:
		unsigned getAACountForAA(std::string aa);
		unsigned getAACountForAA(unsigned aaIndex);
		unsigned getCodonCountForCodon(std::string& codon);
		unsigned getCodonCountForCodon(unsigned codonIndex);
		unsigned getRFPObserved(std::string codon);
		unsigned getRFPObserved(unsigned codonIndex);
		void setRFPObserved(unsigned codonIndex, unsigned value);
		std::vector <unsigned> *getCodonPositions(std::string codon);
		std::vector <unsigned> *getCodonPositions(unsigned index);
		std::vector <unsigned> getRFP_count();
		void setRFP_count(std::vector <unsigned> arg);


		//Other Functions:
		void clear(); //Tested
		bool processSequence(const std::string& sequence);  //Tested TODO: WHY return a bool


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
		static char complimentNucleotide(char ch); //TODO: Testing (c++)
		static std::vector<std::string> aminoAcids(); //Moving to CT, but used in R currently
		static std::vector<std::string> codons(); //Moving to CT, but used in R currently


	protected:
};

#endif // SequenceSummary_H


