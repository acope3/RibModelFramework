#ifndef GENE_H
#define GENE_H


#include "SequenceSummary.h"


#include <string>
#include <vector>
#include <map>

class Gene
{

	private:

		std::string seq; //Gene sequence. Ex: "AATTCAGCT..."
		std::string id; //Gene id: Ex: "YALOC001"
		std::string description; //Additional information about the gene.

	public:


		SequenceSummary geneData;  //TODO: might make private
		std::vector<double> observedSynthesisRateValues; //TODO: make private


		///Constructors & Destructors:
		Gene();
		Gene(std::string _seq, std::string _id, std::string _desc);
		Gene(const Gene& other);
		Gene& operator=(const Gene& rhs);
		bool operator==(const Gene& other) const;
		virtual ~Gene();


		//Data Manipulation Functions:
		std::string getId();
		void setId(std::string _id);
		std::string getDescription();
		void setDescription(std::string _desc);
		std::string getSequence();
		void setSequence(std::string _seq);
		std::vector <unsigned> getRFP_count(); //Only for unit testing.
		void addRFP_count(std::vector <unsigned> RFP_counts);
		SequenceSummary *getSequenceSummary();
		std::vector<double> getObservedSynthesisRateValues(); //exposed to RCPP, tested in C++
		void setObservedSynthesisRateValues(std::vector <double> values); //Only for unit testing.
		double getObservedSynthesisRate(unsigned index);
		unsigned getNumObservedSynthesisSets();
		char getNucleotideAt(unsigned i);


		//Other functions:
		void clear(); // clear the content of object
		unsigned length(); //exposed to RCPP, tested in C++
		Gene reverseComplement(); // return the reverse compliment
		std::string toAASequence();


		//R Section:

#ifndef STANDALONE

		unsigned getAACount(std::string aa);
		unsigned getCodonCount(std::string& codon);
		unsigned getRFPObserved(std::string codon);
		std::vector <unsigned> getCodonPositions(std::string codon);
#endif

	protected:
};

#endif // GENE_H

/*--------------------------------------------------------------------------------------------------
 *                                   !!!RCPP NOTE!!!
 * The two R wrapper functions exist so the user does not have to know about the sequence summary
 * object. Ultimately, SequenceSummary does not need to be exposed - if it is however, these
 * functions could be removed.
 -------------------------------------------------------------------------------------------------*/
