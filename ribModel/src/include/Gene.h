#ifndef GENE_H
#define GENE_H


#include "SequenceSummary.h"


#include <string>
#include <vector>


class Gene
{

	private:

		std::string seq;
		std::string id;
		std::string description;


		void cleanSeq(); // clean the sequence, remove non "AGCT" characters

	public:

		SequenceSummary geneData;  //TODO: might make private
		std::vector<double> observedPhiValues; //TODO: make private


		//Constructors & destructors:
		Gene();
		Gene(std::string _id, std::string _desc, std::string _seq);
		virtual ~Gene();
		Gene(const Gene& other);
		Gene& operator=(const Gene& rhs);


		//Stored data functions:
		std::string getId();
		void setId(std::string _id);
		std::string getDescription();
		void setDescription(std::string _desc);
		std::string getSequence();
		void setSequence(std::string _seq);
		SequenceSummary& getSequenceSummary();
		std::vector<double> getObservedPhiValues();
		void setObservedPhiValues(std::vector <double> values); //FOR TEST THAT ONLY!
		char getNucleotideAt(unsigned i);


		//Other functions:
		void clear(); // clear the content of object
		unsigned length();
		Gene reverseComplement(); // return the reverse compliment
		std::string toAASequence();



		//R wrapper functions:
		void cleanSeqR();
		unsigned getAACount(std::string aa);
		unsigned getCodonCount(std::string& codon);
		void setRFPObserved(unsigned index, unsigned value); // testthat only
		unsigned getRFPObserved(std::string codon);
		std::vector <unsigned> *getCodonPositions(std::string codon);

	protected:
};

#endif // GENE_H

/*--------------------------------------------------------------------------------------------------
 *                                   !!!RCPP NOTE!!!
 * The two R wrapper functions exist so the user does not have to know about the sequence summary
 * object. Ultimately, SequenceSummary does not need to be exposed - if it is however, these
 * functions could be removed.
 -------------------------------------------------------------------------------------------------*/