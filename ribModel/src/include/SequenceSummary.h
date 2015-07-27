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

		unsigned ncodons[64];
		unsigned RFPObserved[64];
		unsigned naa[22];
		std::vector <std::vector <unsigned> > codonPositions;


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


		//Amino acid functions:
		unsigned getAACountForAA(std::string aa);  //DONE!
		unsigned getAACountForAA(unsigned aaIndex); //DONE!

		//Codon functions:
		unsigned getCodonCountForCodon(std::string& codon); //DONE!
		unsigned getCodonCountForCodon(unsigned codonIndex); //DONE!

		//Ribosome footprint functions:
		unsigned getRFPObserved(unsigned codonIndex); //DONE!
		void setRFPObserved(unsigned codonIndex, unsigned value); //DONE!

		//TODO: create RFP functions for a given codon


		//Posistion functions:
		std::vector <unsigned> getCodonPositions(unsigned index); //DONE!


		//Functions to manage stored data:
		void clear(); //DONE!
		bool processSequence(const std::string& sequence); //DONE! //TODO: WHY return a bool



		//Static functions:
		static unsigned AAToAAIndex(std::string aa); //DONE!
		static void AAIndexToCodonRange(unsigned aaIndex, bool forParamVector = false, unsigned aaRange[] = nullptr);
		static void AAToCodonRange(std::string aa, bool forParamVector = false, unsigned aaRange[] = nullptr);
		static std::vector<std::string> AAToCodon(std::string aa, bool forParamVector = false); //TODO: Not used????!


		static std::string CodonToAA(std::string& codon);
		static unsigned CodonToIndex(std::string& codon, bool forParamVector = false);
		static unsigned CodonToAAIndex(std::string& codon);


		static unsigned GetNumCodonsForAA(std::string& aa, bool forParamVector = false);
		static std::string IndexToCodon(unsigned i, bool forParamVector = false);
		static std::string IndexToAA(int aa);







		//R wrapper functions:
		unsigned getAACountForAAR(std::string aa) {aa[0] = (char) std::toupper(aa[0]);	return getAACountForAA(aa);}
		unsigned getAACountForAAIndexR(unsigned aaIndex) {return getAACountForAA(aaIndex);} //TEST THAT ONLY!

		unsigned getCodonCountForCodonR(std::string& codon)
		{
			codon[0] = (char) std::toupper(codon[0]);
			codon[1] = (char) std::toupper(codon[1]);
			codon[2] = (char) std::toupper(codon[2]);

			return (codon.length() != 3) ? -1 : getCodonCountForCodon(codon);
		}

		unsigned getCodonCountForCodonIndexR(unsigned codonIndex) //TEST THAT ONLY!
		{
			return getCodonCountForCodon(codonIndex);
		}


		unsigned getRFPObservedR(unsigned codonIndex) { return getRFPObserved(codonIndex);} //TEST THAT ONLY
		void setRFPObservedR(unsigned codonIndex, unsigned value) { setRFPObserved(codonIndex, value);} //TEST THAT ONLY

		std::vector <unsigned> getCodonPositionsR(unsigned index) { return getCodonPositions(index); } //TEST THAT ONLY




		//static R wrapper functions

		static unsigned AAToAAIndexR(std::string aa) { return AAToAAIndex(aa); } //TEST THAT ONLY


		static std::vector<std::string> aminoAcids() {return AminoAcidArray; }
		static std::vector<std::string> codons()
		{
			std::vector<std::string> RV;
			for (unsigned i = 0; i < 64; i++) RV.push_back(codonArray[i]);
			return RV;
		}
	protected:
};

//TODO: talk about why some functiosn have R wrappers while others don't

#endif // SequenceSummary_H
