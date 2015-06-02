#ifndef SequenceSummary_H
#define SequenceSummary_H

#include <string>
#include <map>
#include <Rcpp.h>
class SequenceSummary
{
    private:

        int ncodons[64];
        int naa[22];

    public:

        //static member variables
        static const char Ser2;
        static const char AminoAcidArray[];
        static const std::string codonArray[];
        static const std::map<char, int> aaToIndex;

        explicit SequenceSummary();
        SequenceSummary(const std::string& sequence);
        virtual ~SequenceSummary();
        SequenceSummary(const SequenceSummary& other);
        SequenceSummary& operator=(const SequenceSummary& other);


        void processSequence(const std::string& sequence);
        void clear();

        int getAAcountForAA(char aa) {return naa[aaToIndex.find(aa)->second];}
        int getAAcountForAA(int aa) {return naa[aa];}

        int getCodonCountForCodon(std::string& codon) {return ncodons[SequenceSummary::CodonToIndex(codon)];}
        int getCodonCountForCodon(int codon) {return ncodons[codon];}


				//R Wrapper Functions
				int getAAcount(char aa) {return getAAcountForAA(aa);}
				int getCodonCount(std::string& codon) {return getCodonCountForCodon(codon);}



        //statics
        static char CodonToAA(std::string& codon);
        static unsigned GetNumCodonsForAA(const char& aa, bool forParamVector = false);
        static unsigned CodonToIndex(std::string& codon);
        static std::string IndexToCodon(unsigned i);
        static unsigned CodonToAAIndex(std::string& codon);
        static char IndexToAA(int aa);
        static void AAindexToCodonRange(unsigned aaIndex, bool forParamVector = false, unsigned aaRange[] = nullptr);
        static void AAToCodonRange(char aa, bool forParamVector = false, unsigned aaRange[] = nullptr);


    protected:



};


RCPP_MODULE(SequenceSummary) 
{
		using namespace Rcpp;
    class_<SequenceSummary>( "SequenceSummary" )
    .constructor()
		.constructor<std::string>()
//		.constructor<SequenceSummary>() //custom object...How?
		//operator overloading????
				
    .method("getCodonCount", &SequenceSummary::getCodonCount)
		.method("getAAcount", &SequenceSummary::getAAcount)

		.method("processSequence", &SequenceSummary::processSequence)
		.method("clear", &SequenceSummary::clear)
		; 
}

#endif // SequenceSummary_H
