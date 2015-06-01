#ifndef GENE_H
#define GENE_H

#include "../include/SequenceSummary.h"

#include <string>
#include <vector>


class Gene
{

    private:
        //member variables
        std::string seq;
        std::string id;
        std::string description;

        //double currentPhiValue;
        //double currentLikelihood;

        void cleanSeq(); // clean the sequence, remove non "AGCT" charactors

    public:
        SequenceSummary geneData;

        //constructor/destructors
        Gene();
        Gene(std::string _id, std::string _desc, std::string _seq);
        virtual ~Gene();
        Gene(const Gene& other);
        Gene& operator=(const Gene& rhs);

        void clear(); // clear the content of object
        int length() {return seq.size();}

        Gene reverseCompliment(); // return the reverse compliment
        std::string toAAsequence();

        std::string toString() {return ">" + id + "\n" + seq + "\n";}

        //getter/setter
        std::string getId() {return id;}
        std::string getDescription() {return description;}
        std::string getSequence() {return seq;}
        char getNucleotideAt(int i) {return seq[i];}

        void setId(std::string _id) { id = _id;}
        void setDescription(std::string _desc) {description = _desc;}
        void setSequence(std::string _seq);
        SequenceSummary& getSequenceSummary() {return geneData;}


    protected:

        //static member variables




};



#endif // GENE_H
