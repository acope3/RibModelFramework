#ifndef SequenceSummary_H
#define SequenceSummary_H


#include <string>
#include <map>
#include <algorithm>
#include <cctype>
#include <vector>
#include <array>


class SequenceSummary
{
	private:

		unsigned ncodons[64];
		unsigned RFPObserved[64];
		unsigned naa[22];
		std::vector <std::vector <unsigned>> codonPositions;


	public:

		//Static member variables:
		static const std::string Ser2;
		static const std::vector<std::string> AminoAcidArray;
		static const std::string codonArray[];
		static const std::string codonArrayParameter[];
		static const std::map<std::string, unsigned> aaToIndex;
		static const std::map<std::string, unsigned> codonToIndexWithReference;
		static const std::map<std::string, unsigned> codonToIndexWithoutReference;


		//Constructors & destructors:
		explicit SequenceSummary();
		SequenceSummary(const std::string& sequence);
		virtual ~SequenceSummary();
		SequenceSummary(const SequenceSummary& other);
		SequenceSummary& operator=(const SequenceSummary& other);


		//AA, codon, RFP, & position functions:
		unsigned getAACountForAA(std::string aa);
		unsigned getAACountForAA(unsigned aaIndex);
		unsigned getCodonCountForCodon(std::string& codon);
		unsigned getCodonCountForCodon(unsigned codonIndex);
		unsigned getRFPObserved(std::string codon);
		unsigned getRFPObserved(unsigned codonIndex);
		void setRFPObserved(unsigned codonIndex, unsigned value);
		std::vector <unsigned> getCodonPositions(std::string codon);
		std::vector <unsigned> getCodonPositions(unsigned index);



		//Functions to manage stored data:
		void clear();
		bool processSequence(const std::string& sequence);  //TODO: WHY return a bool


		//Static functions:
		static unsigned AAToAAIndex(std::string aa);
		static std::array<unsigned, 2> AAIndexToCodonRange(unsigned aaIndex, bool forParamVector = false);
		static std::array<unsigned, 2> AAToCodonRange(std::string aa, bool forParamVector = false);
		static std::vector<std::string> AAToCodon(std::string aa, bool forParamVector = false);
		static std::string codonToAA(std::string& codon);
		static unsigned codonToIndex(std::string& codon, bool forParamVector = false);
		static unsigned codonToAAIndex(std::string& codon);
		static std::string indexToAA(unsigned aaIndex);
		static std::string indexToCodon(unsigned index, bool forParamVector = false);
		static unsigned GetNumCodonsForAA(std::string& aa, bool forParamVector = false);
		static char complimentNucleotide(char ch); //TODO: impliment for test that!
		static std::vector<std::string> aminoAcids();
		static std::vector<std::string> codons();



		//R wrapper functions:
		unsigned getAACountForAAR(std::string aa);
		unsigned getAACountForAAIndexR(unsigned aaIndex); //TEST THAT ONLY!
		unsigned getCodonCountForCodonR(std::string& codon);
		unsigned getCodonCountForCodonIndexR(unsigned codonIndex); //TEST THAT ONLY!
		unsigned getRFPObservedForCodonR(std::string codon);
		unsigned getRFPObservedForCodonIndexR(unsigned codonIndex); //TEST THAT ONLY!
		std::vector <unsigned> getCodonPositionsForCodonR(std::string codon);
		std::vector <unsigned> getCodonPositionsForCodonIndexR(unsigned codonIndex); //TEST THAT ONLY!


	protected:
};

#endif // SequenceSummary_H

/*--------------------------------------------------------------------------------------------------
 *                                   !!!RCPP NOTE!!!
 * All functions are exposed to R. However, some functions are only exposed for the purpose of
 * unit testing (Test That). These functions can be accessed in R, but do not check indices or
 * potential problems with user input (such as codon strings being lower case). The two functions
 * that have R wrappers that only call the C++ function had to be implemented because RCPP cannot
 * deal with overloaded functions.
 -------------------------------------------------------------------------------------------------*/
