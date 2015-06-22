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

	public:

		//static member variables
		static const char Ser2;
		static const std::vector<char> AminoAcidArray;
		static const std::string codonArray[];
		static const std::string codonArrayParameter[];
		static const std::map<char, int> aaToIndex;

		//constructors and destructors
		explicit SequenceSummary();
		SequenceSummary(const std::string& sequence);
		virtual ~SequenceSummary();
		SequenceSummary(const SequenceSummary& other);
		SequenceSummary& operator=(const SequenceSummary& other);

		bool processSequence(const std::string& sequence);
		void clear();

		int getAAcountForAA(char aa) {return naa[aaToIndex.find(aa)->second];}
		int getAAcountForAA(int aa) 
		{
			return naa[aa];
		}
		int getCodonCountForCodon(std::string& codon) {return ncodons[SequenceSummary::CodonToIndex(codon)];}
		int getCodonCountForCodon(int codon) 
		{
			return ncodons[codon];	
		}

		//R Wrapper Functions
		int getAAcount(char aa) {aa = std::toupper(aa);	return getAAcountForAA(aa);}
		int getCodonCount(std::string& codon)
		{
			codon[0] = std::toupper(codon[0]);
			codon[1] = std::toupper(codon[1]);
			codon[2] = std::toupper(codon[2]);

			return (codon.length() != 3) ? -1 : getCodonCountForCodon(codon);
		}



		//static functions
		static unsigned AAToAAIndex(char aa);
		static char CodonToAA(std::string& codon);
		static unsigned GetNumCodonsForAA(char& aa, bool forParamVector = false);
		static unsigned CodonToIndex(std::string& codon, bool forParamVector = false);
		static std::string IndexToCodon(unsigned i, bool forParamVector = false);
		static unsigned CodonToAAIndex(std::string& codon);
		static char IndexToAA(int aa);
		static void AAindexToCodonRange(unsigned aaIndex, bool forParamVector = false, unsigned aaRange[] = nullptr);
		static void AAToCodonRange(char aa, bool forParamVector = false, unsigned aaRange[] = nullptr);
		static std::vector<std::string> AAToCodon(char aa, bool forParamVector = false);

		//static R wrapper functions
		static std::vector<char> aminoAcids() {return AminoAcidArray; }

	protected:



};


#endif // SequenceSummary_H
