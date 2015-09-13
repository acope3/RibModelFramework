#ifndef SequenceSummary_H
#define SequenceSummary_H

#include "CodonTable.h"
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
		std::vector <unsigned> naa;
		std::vector <std::vector <unsigned>> codonPositions;


	public:

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
		static char complimentNucleotide(char ch);



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
