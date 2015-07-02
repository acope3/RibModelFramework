#ifndef SequenceSummary_H
#define SequenceSummary_H

#include <string>
#include <map>
#include <algorithm>
#include <cctype>
#include <vector>
class SequenceSummary
{
	private:

		int ncodons[64];
		int naa[22];
		unsigned RFPObserved[64];
		unsigned numCodonsInMRNA[64];
	public:

		//static member variables
		static const std::string Ser2;
		static const std::vector<std::string> AminoAcidArray;
		static const std::string codonArray[];
		static const std::string codonArrayParameter[];
		static const std::map<std::string, int> aaToIndex;

		//constructors and destructors
		explicit SequenceSummary();
		SequenceSummary(const std::string& sequence);
		virtual ~SequenceSummary();
		SequenceSummary(const SequenceSummary& other);
		SequenceSummary& operator=(const SequenceSummary& other);

		bool processSequence(const std::string& sequence);
		void clear();

		int getAAcountForAA(std::string aa) {return naa[aaToIndex.find(aa)->second];}
		int getAAcountForAA(int aa) {return naa[aa];}
		int getCodonCountForCodon(std::string& codon) {return ncodons[SequenceSummary::CodonToIndex(codon)];}
		int getCodonCountForCodon(int codon) {return ncodons[codon];}

		unsigned getRFPObserved(unsigned codonIndex) {return RFPObserved[codonIndex];}
		unsigned getNumCodonsInMRNA(unsigned codonIndex) {return numCodonsInMRNA[codonIndex];}
		void setRFPObserved(unsigned codonIndex, unsigned value) {RFPObserved[codonIndex] = value;}
		void setNumCodonsInMRNA(unsigned codonIndex, unsigned value) {numCodonsInMRNA[codonIndex] = value;}

		//R Wrapper Functions
		int getAAcount(std::string aa) {aa[0] = (char) std::toupper(aa[0]);	return getAAcountForAA(aa);}
		int getCodonCount(std::string& codon)
		{
			codon[0] = (char) std::toupper(codon[0]);
			codon[1] = (char) std::toupper(codon[1]);
			codon[2] = (char) std::toupper(codon[2]);

			return (codon.length() != 3) ? -1 : getCodonCountForCodon(codon);
		}



		//static functions
		static unsigned AAToAAIndex(std::string aa);
		static std::string CodonToAA(std::string& codon);
		static unsigned GetNumCodonsForAA(std::string& aa, bool forParamVector = false);
		static unsigned CodonToIndex(std::string& codon, bool forParamVector = false);
		static std::string IndexToCodon(unsigned i, bool forParamVector = false);
		static unsigned CodonToAAIndex(std::string& codon);
		static std::string IndexToAA(int aa);
		static void AAindexToCodonRange(unsigned aaIndex, bool forParamVector = false, unsigned aaRange[] = nullptr);
		static void AAToCodonRange(std::string aa, bool forParamVector = false, unsigned aaRange[] = nullptr);
		static std::vector<std::string> AAToCodon(std::string aa, bool forParamVector = false);

		//static R wrapper functions
		static std::vector<std::string> aminoAcids() {return AminoAcidArray; }

	protected:



};


#endif // SequenceSummary_H
