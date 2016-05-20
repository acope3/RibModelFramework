#include "include/Testing.h"

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

int testSequenceSummary()
{
    SequenceSummary SS("ATGCTCATTCTCACTGCTGCCTCGTAG");
    int error = 0;
    int globalError = 0;

    //----------------------------//
    //------ Clear Function ------//
    //----------------------------//
    SS.clear();
    for (unsigned i = 0; i < 64; i++)
    {
        if (0 != SS.getCodonCountForCodon(i))
        {
            std::cerr << "Problem with Sequence Summary \"clear\" function.\n";
            std::cerr << "Problem at codon index" << i << "\n";
            error = 1;
            globalError = 1;
        }
    }
    for (unsigned i = 0; i < 22; i++)
    {
        if (0 != SS.getAACountForAA(i))
        {
            std::cerr << "Problem with Sequence Summary \"clear\" function.\n";
            std::cerr << "Problem at amino acid index" << i << "\n";
            error = 1;
            globalError = 1;
        }
    }

    if (!error)
        std::cout << "Sequence Summary clear --- Pass\n";
    else
        error = 0; //Reset for next function.

    //--------------------------------------//
    //------ ProcessSequence Function ------//
    //--------------------------------------//
    SS.processSequence("ATGCTCATTCTCACTGCTGCCTCGTAG");

    if (1 != SS.getAACountForAA("I"))
    {
        std::cerr << "Problem with Sequence Summary \"processSequence\" function.\n";
        std::cerr << "Problem with amino acid \"I\".";
        std::cerr << "I is in the sequence once, but is returning " << SS.getAACountForAA("I") << "\n";
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA("T"))
    {
        std::cerr << "Problem with Sequence Summary \"processSequence\" function.\n";
        std::cerr << "Problem with amino acid \"T\".";
        std::cerr << "T is in the sequence once, but is returning " << SS.getAACountForAA("T") << "\n";
        error = 1;
        globalError = 1;
    }

    std::string codon = "ATT";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        std::cerr << "Problem with Sequence Summary \"processSequence\" function.\n";
        std::cerr << "Problem with codon \"ATT\".";
        std::cerr << "ATT is in the sequence once, but is returning " << SS.getCodonCountForCodon(codon) << "\n";
        error = 1;
        globalError = 1;
    }

    codon = "ACT";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        std::cerr << "Problem with Sequence Summary \"processSequence\" function.\n";
        std::cerr << "Problem with codon \"ACT\".";
        std::cerr << "ACT is in the sequence once, but is returning " << SS.getCodonCountForCodon(codon) << "\n";
        error = 1;
        globalError = 1;
    }

    codon = "GCT";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        std::cerr << "Problem with Sequence Summary \"processSequence\" function.\n";
        std::cerr << "Problem with codon \"GCT\".";
        std::cerr << "GCT is in the sequence once, but is returning " << SS.getCodonCountForCodon(codon) << "\n";
        error = 1;
        globalError = 1;
    }

    codon = "GCC";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        std::cerr << "Problem with Sequence Summary \"processSequence\" function.\n";
        std::cerr << "Problem with codon \"GCC\".";
        std::cerr << "GCC is in the sequence once, but is returning " << SS.getCodonCountForCodon(codon) << "\n";
        error = 1;
        globalError = 1;
    }

    std::vector <unsigned> *tmp;
    tmp = SS.getCodonPositions("CTC");
    if ((1 != tmp -> at(0)) && (3 != tmp -> at(1)))
    {
        std::cerr << "Codon CTC should be found at position 1 and 3(zero indexed), but is";
        std::cerr << "found at these locations:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonPositions("ATT");
    if (2 != tmp -> at(0))
    {
        std::cerr << "Codon ATT should be found at position 2(zero indexed), but is";
        std::cerr << "found at these locations:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Sequence Summary processSequence --- Pass\n";
    else
        error = 0; //Reset for next function.

    //------------------------------------------//
    //------ complimentNucleotide Function------//
    //------------------------------------------//

    if ('T' != SequenceSummary::complimentNucleotide('A'))
    {
        std::cerr << "The compliment of A should be T\n";
        error = 1;
        globalError = 1;
    }

    if ('A' != SequenceSummary::complimentNucleotide('T'))
    {
        std::cerr << "The compliment of T should be A\n";
        error = 1;
        globalError = 1;
    }

    if ('G' != SequenceSummary::complimentNucleotide('C'))
    {
        std::cerr << "The compliment of C should be G\n";
        error = 1;
        globalError = 1;
    }

    if ('C' != SequenceSummary::complimentNucleotide('G'))
    {
        std::cerr << "The compliment of G should be C\n";
        error = 1;
        globalError = 1;
    }

    if ('C' != SequenceSummary::complimentNucleotide('Q'))
    {
        std::cerr << "The compliment of Q should be C\n";
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Sequence Summary complimentNucleotide --- Pass\n";
    else
        error = 0; //Reset for next function.

    //----------------------------------------------//
    //------ getAACountForAA(string) Function ------//
    //----------------------------------------------//

    if (1 != SS.getAACountForAA("M"))
    {
        std::cerr << "Error with getAACountForAA(string) for amino acid M.\n";
        std::cerr << "Should return 1, returns " << SS.getAACountForAA("M") << "\n";
        error = 1;
        globalError = 1;
    }

    if (2 != SS.getAACountForAA("L"))
    {
        std::cerr << "Error with getAACountForAA(string) for amino acid L.\n";
        std::cerr << "Should return 2, returns " << SS.getAACountForAA("L") << "\n";
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA("I"))
    {
        std::cerr << "Error with getAACountForAA(string) for amino acid I.\n";
        std::cerr << "Should return 1, returns " << SS.getAACountForAA("I") << "\n";
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA("T"))
    {
        std::cerr << "Error with getAACountForAA(string) for amino acid T.\n";
        std::cerr << "Should return 1, returns " << SS.getAACountForAA("T") << "\n";
        error = 1;
        globalError = 1;
    }

    if (2 != SS.getAACountForAA("A"))
    {
        std::cerr << "Error with getAACountForAA(string) for amino acid A.\n";
        std::cerr << "Should return 2, returns " << SS.getAACountForAA("A") << "\n";
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA("S"))
    {
        std::cerr << "Error with getAACountForAA(string) for amino acid S.\n";
        std::cerr << "Should return 1, returns " << SS.getAACountForAA("S") << "\n";
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA("X"))
    {
        std::cerr << "Error with getAACountForAA(string) for amino acid X.\n";
        std::cerr << "Should return 1, returns " << SS.getAACountForAA("X") << "\n";
        error = 1;
        globalError = 1;
    }

    if (0 != SS.getAACountForAA("G"))
    {
        std::cerr << "Error with getAACountForAA(string) for amino acid G.\n";
        std::cerr << "Should return 0, returns " << SS.getAACountForAA("G") << "\n";
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Sequence Summary getAACountForAA(string) --- Pass\n";
    else
        error = 0; //Reset for next function.

    //---------------------------------------------//
    //------ getAACountForAA(index) Function ------//
    //---------------------------------------------//
    if (1 != SS.getAACountForAA(10))
    {
        std::cerr << "Error with getAACountForAA(index) for amino acid M (index 10).\n";
        std::cerr << "Should return 1, returns " << SS.getAACountForAA(10) << "\n";
        error = 1;
        globalError = 1;
    }

    if (2 != SS.getAACountForAA(9))
    {
        std::cerr << "Error with getAACountForAA(index) for amino acid L (index 9).\n";
        std::cerr << "Should return 2, returns " << SS.getAACountForAA(9) << "\n";
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA(7))
    {
        std::cerr << "Error with getAACountForAA(index) for amino acid I (index 7).\n";
        std::cerr << "Should return 1, returns " << SS.getAACountForAA(7) << "\n";
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA(16))
    {
        std::cerr << "Error with getAACountForAA(index) for amino acid T (index 16).\n";
        std::cerr << "Should return 1, returns " << SS.getAACountForAA(16) << "\n";
        error = 1;
        globalError = 1;
    }

    if (2 != SS.getAACountForAA(0))
    {
        std::cerr << "Error with getAACountForAA(index) for amino acid A (index 0).\n";
        std::cerr << "Should return 2, returns " << SS.getAACountForAA(0) << "\n";
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA(15))
    {
        std::cerr << "Error with getAACountForAA(index) for amino acid S (index 15).\n";
        std::cerr << "Should return 1, returns " << SS.getAACountForAA(15) << "\n";
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getAACountForAA(21))
    {
        std::cerr << "Error with getAACountForAA(index) for amino acid X (index 21).\n";
        std::cerr << "Should return 1, returns " << SS.getAACountForAA(21) << "\n";
        error = 1;
        globalError = 1;
    }

    if (0 != SS.getAACountForAA(2))
    {
        std::cerr << "Error with getAACountForAA(index) for amino acid D (index 2).\n";
        std::cerr << "Should return 0, returns " << SS.getAACountForAA(2) << "\n";
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Sequence Summary getAACountForAA(index) --- Pass\n";
    else
        error = 0; //Reset for next function.

    //--------------------------------------------//
    //------ getCodonCountsForCodon(string) ------//
    //--------------------------------------------//

    codon = "ATG";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        std::cerr << "Error with getCodonCountForCodon(string) for " << codon << ".\n";
        std::cerr << "Should return 1, but returns " << SS.getCodonCountForCodon(codon) << "\n";
        error = 1;
        globalError = 1;
    }

    codon = "CTC";
    if (2 != SS.getCodonCountForCodon(codon))
    {
        std::cerr << "Error with getCodonCountForCodon(string) for " << codon << ".\n";
        std::cerr << "Should return 2, but returns " << SS.getCodonCountForCodon(codon) << "\n";
        error = 1;
        globalError = 1;
    }

    codon = "ATT";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        std::cerr << "Error with getCodonCountForCodon(string) for " << codon << ".\n";
        std::cerr << "Should return 1, but returns " << SS.getCodonCountForCodon(codon) << "\n";
        error = 1;
        globalError = 1;
    }

    codon = "ACT";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        std::cerr << "Error with getCodonCountForCodon(string) for " << codon << ".\n";
        std::cerr << "Should return 1, but returns " << SS.getCodonCountForCodon(codon) << "\n";
        error = 1;
        globalError = 1;
    }

    codon = "GCT";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        std::cerr << "Error with getCodonCountForCodon(string) for " << codon << ".\n";
        std::cerr << "Should return 1, but returns " << SS.getCodonCountForCodon(codon) << "\n";
        error = 1;
        globalError = 1;
    }

    codon = "GCC";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        std::cerr << "Error with getCodonCountForCodon(string) for " << codon << ".\n";
        std::cerr << "Should return 1, but returns " << SS.getCodonCountForCodon(codon) << "\n";
        error = 1;
        globalError = 1;
    }

    codon = "TCG";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        std::cerr << "Error with getCodonCountForCodon(string) for " << codon << ".\n";
        std::cerr << "Should return 1, but returns " << SS.getCodonCountForCodon(codon) << "\n";
        error = 1;
        globalError = 1;
    }

    codon = "TAG";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        std::cerr << "Error with getCodonCountForCodon(string) for " << codon << ".\n";
        std::cerr << "Should return 1, but returns " << SS.getCodonCountForCodon(codon) << "\n";
        error = 1;
        globalError = 1;
    }

    codon = "AAA";
    if (0 != SS.getCodonCountForCodon(codon))
    {
        std::cerr << "Error with getCodonCountForCodon(string) for " << codon << ".\n";
        std::cerr << "Should return 0, but returns " << SS.getCodonCountForCodon(codon) << "\n";
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Sequence Summary getCodonCountsForCodon(string) --- Pass\n";
    else
        error = 0; //Reset for next function.

    //----------------------------------------------------//
    //------ getCodonCountsForCodon(index) Function ------//
    //----------------------------------------------------//

    if (1 != SS.getCodonCountForCodon(29))
    {
        std::cerr << "Error with getCodonCountForCodon(index) for codon \"ATG\" (index 29).\n";
        std::cerr << "Should return 1, but returns " << SS.getCodonCountForCodon(29) << "\n";
        error = 1;
        globalError = 1;
    }

    if (2 != SS.getCodonCountForCodon(24))
    {
        std::cerr << "Error with getCodonCountForCodon(index) for codon \"CTC\" (index 24).\n";
        std::cerr << "Should return 2, but returns " << SS.getCodonCountForCodon(24) << "\n";
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getCodonCountForCodon(20))
    {
        std::cerr << "Error with getCodonCountForCodon(index) for codon \"ATT\" (index 20).\n";
        std::cerr << "Should return 1, but returns " << SS.getCodonCountForCodon(20) << "\n";
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getCodonCountForCodon(51))
    {
        std::cerr << "Error with getCodonCountForCodon(index) for codon \"ACT\" (index 51).\n";
        std::cerr << "Should return 1, but returns " << SS.getCodonCountForCodon(51) << "\n";
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getCodonCountForCodon(3))
    {
        std::cerr << "Error with getCodonCountForCodon(index) for codon \"GCT\" (index 3).\n";
        std::cerr << "Should return 1, but returns " << SS.getCodonCountForCodon(3) << "\n";
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getCodonCountForCodon(1))
    {
        std::cerr << "Error with getCodonCountForCodon(index) for codon \"GCC\" (index 1).\n";
        std::cerr << "Should return 1, but returns " << SS.getCodonCountForCodon(1) << "\n";
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getCodonCountForCodon(46))
    {
        std::cerr << "Error with getCodonCountForCodon(index) for codon \"TCG\" (index 46).\n";
        std::cerr << "Should return 1, but returns " << SS.getCodonCountForCodon(46) << "\n";
        error = 1;
        globalError = 1;
    }

    if (1 != SS.getCodonCountForCodon(62))
    {
        std::cerr << "Error with getCodonCountForCodon(index) for codon \"TAG\" (index 62).\n";
        std::cerr << "Should return 1, but returns " << SS.getCodonCountForCodon(62) << "\n";
        error = 1;
        globalError = 1;
    }

    if (0 != SS.getCodonCountForCodon(2))
    {
        std::cerr << "Error with getCodonCountForCodon(index) for codon \"AAA\" (index 2).\n";
        std::cerr << "Should return 1, but returns " << SS.getCodonCountForCodon(2) << "\n";
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Sequence Summary getCodonCountsForCodon(index) --- Pass\n";
    else
        error = 0; //Reset for next function.

    //--------------------------------------------------------------//
    //------ getRFPObserved(string) & setRFPObserved Function ------//
    //--------------------------------------------------------------//

    SS.setRFPObserved(4, 35);
    SS.setRFPObserved(16,45);
    SS.setRFPObserved(54,2);
    SS.setRFPObserved(45,0);

    if (35 != SS.getRFPObserved("TGC"))
    {
        std::cerr << "Error in getRFPObserved(string) or setRFPObserved for codon \"TGC\".\n";
        std::cerr << "should return 35, but returns " << SS.getRFPObserved("TGC") << "\n";
        error = 1;
        globalError = 1;
    }

    if (45 != SS.getRFPObserved("CAC"))
    {
        std::cerr << "Error in getRFPObserved(string) or setRFPObserved for codon \"CAC\".\n";
        std::cerr << "should return 45, but returns " << SS.getRFPObserved("CAC") << "\n";
        error = 1;
        globalError = 1;
    }

    if (2 != SS.getRFPObserved("GTG"))
    {
        std::cerr << "Error in getRFPObserved(string) or set RFPObserved for codon \"GTG\".\n";
        std::cerr << "should return 2, but returns " << SS.getRFPObserved("GTG") << "\n";
        error = 1;
        globalError = 1;
    }

    if (0 != SS.getRFPObserved("TCC"))
    {
        std::cerr << "Error in getRFPObserved(string) or setRFPObserved for codon \"TCC\".\n";
        std::cerr << "should return 0, but returns " << SS.getRFPObserved("TCC") << "\n";
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Sequence Summary getRFPObserved(string) & setRFPObserved --- Pass\n";
    else
        error = 0; //Reset for next function.

    //--------------------------------------------//
    //------ getRFPObserved(index) Function ------//
    //--------------------------------------------//

    SS.setRFPObserved(0,45);
    SS.setRFPObserved(1,52);
    SS.setRFPObserved(2,63);
    SS.setRFPObserved(60,23);

    if (45 != SS.getRFPObserved(0))
    {
        std::cerr << "Error with getRFPObserved(index) for codon index 0.\n";
        std::cerr << "should return 45, but returns " << SS.getRFPObserved(0) << "\n";
        error = 1;
        globalError = 1;
    }

    if (52 != SS.getRFPObserved(1))
    {
        std::cerr << "Error with getRFPObserved(index) for codon index 1.\n";
        std::cerr << "should return 52, but returns " << SS.getRFPObserved(1) << "\n";
        error = 1;
        globalError = 1;
    }

    if (63 != SS.getRFPObserved(2))
    {
        std::cerr << "Error with getRFPObserved(index) for codon index 2.\n";
        std::cerr << "should return 63, but returns " << SS.getRFPObserved(2) << "\n";
        error = 1;
        globalError = 1;
    }

    if (23 != SS.getRFPObserved(60))
    {
        std::cerr << "Error with getRFPObserved(index) for codon index 60.\n";
        std::cerr << "should return 23, but returns " << SS.getRFPObserved(60) << "\n";
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Sequence Summary getRFPObserved(index) --- Pass\n";
    else
        error = 0; //Reset for next function.

    //------------------------------------------------//
    //------ getCodonPositions(string) Function ------//
    //------------------------------------------------//

    tmp = SS.getCodonPositions("ATG");
    if (tmp -> at(0) != 0 || tmp -> size() != 1)
    {
        std::cerr << "Error with getCodonPositions(string) for codon \"ATG\".\n";
        std::cerr << "Should return 0, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonPositions("CTC");
    if (tmp -> at(0) != 1 || tmp -> at(1) != 3|| tmp -> size() != 2)
    {
        std::cerr << "Error with getCodonPositions(string) for codon \"CTC\".\n";
        std::cerr << "Should return 1 and 3, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonPositions("ATT");
    if (tmp -> at(0) != 2 || tmp -> size() != 1)
    {
        std::cerr << "Error with getCodonPositions(string) for codon \"ATT\".\n";
        std::cerr << "Should return 2, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonPositions("ACT");
    if (tmp -> at(0) != 4 || tmp -> size() != 1)
    {
        std::cerr << "Error with getCodonPositions(string) for codon \"ACT\".\n";
        std::cerr << "Should return 4, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }


    tmp = SS.getCodonPositions("GCT");
    if (tmp -> at(0) != 5 || tmp -> size() != 1)
    {
        std::cerr << "Error with getCodonPositions(string) for codon \"GCT\".\n";
        std::cerr << "Should return 5, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonPositions("GCC");
    if (tmp -> at(0) != 6 || tmp -> size() != 1)
    {
        std::cerr << "Error with getCodonPositions(string) for codon \"GCC\".\n";
        std::cerr << "Should return 6, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonPositions("TCG");
    if (tmp -> at(0) != 7 || tmp -> size() != 1)
    {
        std::cerr << "Error with getCodonPositions(string) for codon \"TCG\".\n";
        std::cerr << "Should return 7, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonPositions("TAG");
    if (tmp -> at(0) != 8 || tmp -> size() != 1)
    {
        std::cerr << "Error with getCodonPositions(string) for codon \"TAG\".\n";
        std::cerr << "Should return 8, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonPositions("GTG");
    if (tmp -> size() != 0)
    {
        std::cerr << "Error with getCodonPositions(string) for codon \"GTG\".\n";
        std::cerr << "Should return an empty vector, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Sequence Summary getCodonPositions(string) --- Pass\n";
    else
        error = 0; //Reset for next function.

    //-----------------------------------------------//
    //------ getCodonPositions(index) Function ------//
    //-----------------------------------------------//

    tmp = SS.getCodonPositions(29);
    if (tmp -> at(0) != 0 || tmp -> size() != 1)
    {
        std::cerr << "Error with getCodonPositions(index) for codon index 29.\n";
        std::cerr << "Should return 0, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonPositions(24);
    if (tmp -> at(0) != 1 || tmp -> at(1) != 3 || tmp -> size() != 2)
    {
        std::cerr << "Error with getCodonPositions(index) for codon index 24.\n";
        std::cerr << "Should return 1 and 3, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonPositions(20);
    if (tmp -> at(0) != 2 || tmp -> size() != 1)
    {
        std::cerr << "Error with getCodonPositions(index) for codon index 20.\n";
        std::cerr << "Should return 2, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonPositions(51);
    if (tmp -> at(0) != 4 || tmp -> size() != 1)
    {
        std::cerr << "Error with getCodonPositions(index) for codon index 51.\n";
        std::cerr << "Should return 4, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonPositions(3);
    if (tmp -> at(0) != 5 || tmp -> size() != 1)
    {
        std::cerr << "Error with getCodonPositions(index) for codon index 3.\n";
        std::cerr << "Should return 4, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonPositions(1);
    if (tmp -> at(0) != 6 || tmp -> size() != 1)
    {
        std::cerr << "Error with getCodonPositions(index) for codon index 1.\n";
        std::cerr << "Should return 4, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonPositions(46);
    if (tmp -> at(0) != 7 || tmp -> size() != 1)
    {
        std::cerr << "Error with getCodonPositions(index) for codon index 46.\n";
        std::cerr << "Should return 7, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonPositions(62);
    if (tmp -> at(0) != 8 || tmp -> size() != 1)
    {
        std::cerr << "Error with getCodonPositions(index) for codon index 62.\n";
        std::cerr << "Should return 8, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp = SS.getCodonPositions(54);
    if (tmp -> size() != 0)
    {
        std::cerr << "Error with getCodonPositions(index) for codon index 54.\n";
        std::cerr << "Should return an empty vector, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) << "\n";
        }
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Sequence Summary getCodonPositions(index) --- Pass\n";
    else
        error = 0; //Reset for next function.

    //---------------------------------------//
    //------ get/setRFP_count Function ------//
    //---------------------------------------//
    std::vector <unsigned> tmp2 = SS.getRFP_count();

    if (0 != tmp2.size())
    {
        std::cerr << "Error with getRFP_count. Function should return an empty vector but returns:\n";
        for (unsigned i = 0; i < tmp2.size(); i++)
        {
            std::cerr << tmp2[i] << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp2 = {1, 2, 3, 4, 5};
    SS.setRFP_count(tmp2);

    if (SS.getRFP_count() != tmp2)
    {
        std::cerr << "Error in getRFP_count or setRFP_count. Function should return 1, 2, 3, 4, 5, but returns:\n";
        for (unsigned i = 0; i < tmp2.size(); i++)
        {
            std::cerr << tmp2[i] << "\n";
        }
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Sequence Summary get/setRFP_count --- Pass\n";

    // No need to reset error

    return globalError;
}


int testGene()
{
    Gene testGene;
    int error = 0;
    int globalError = 0;

    //--------------------------------//
    //------ get/setId Function ------//
    //--------------------------------//
    testGene.setId("testGene");

    if (testGene.getId() != "testGene") {
        std::cerr << "Error in setId or getId.\n";
        globalError = 1;
    }
    else
        std::cout << "Gene get/setId --- Pass\n";

    //-----------------------------------------//
    //------ get/setDescription Function ------//
    //-----------------------------------------//
    testGene.setDescription("Test Gene for Unit Testing");

    if (testGene.getDescription() != "Test Gene for Unit Testing")
    {
        std::cerr << "Error in setDescription or getDescription.\n";
        globalError = 1;
    }
    else
        std::cout << "Gene get/setDescription --- Pass\n";

    //--------------------------------------//
    //------ get/setSequence Function ------//
    //--------------------------------------//
    testGene.setSequence("ATGCTCATTCTCACTGCTGCCTCGTAG");

    if (testGene.getSequence() != "ATGCTCATTCTCACTGCTGCCTCGTAG")
    {
        std::cerr << "Error in setSequence or getSequence.\n";
        globalError = 1;
    }
    else
        std::cout << "Gene get/setSequence --- Pass\n";

    //---------------------------------------//
    //------ get/addRFP_count Function ------//
    //---------------------------------------//
    std::vector <unsigned> rfp_counts = {0, 1, 1};

    testGene.addRFP_count(rfp_counts);
    if (testGene.getRFP_count() != rfp_counts)
    {
        std::cerr << "Error in getRFP_count or addRFP_count.\n";
        globalError = 1;
    }
    else
        std::cout << "Gene get/addRFP_count --- Pass\n";

    //--------------------------------//
    //------ getSequenceSummary ------//
    //--------------------------------//
    SequenceSummary SS("ATGCTCATTCTCACTGCTGCCTCGTAG");
    SequenceSummary *GeneSS = testGene.getSequenceSummary();
    for (unsigned i = 0; i < 64; i++)
    {
        if (SS.getCodonCountForCodon(i) != GeneSS->getCodonCountForCodon(i))
        {
            std::cerr << "Error with getSequenceSummary. Codon counts are incorrect";
            std::cerr << " for codon " << i << ", " << SequenceSummary::codonArray[i] << ".\n";
            std::cerr << "Should return " << SS.getCodonCountForCodon(i) << ", but returns" << GeneSS->getCodonCountForCodon(i) << "\n";
            error = 1;
            globalError = 1;
        }
    }

    for (unsigned i = 0; i < 64; i++)
    {
        if (SS.getRFPObserved(i) != GeneSS->getRFPObserved(i))
        {
            std::cerr << "Error with getSequenceSummary. RFP observed is incorrect";
            std::cerr << " for codon " << i << ".\n";
            std::cerr << "Should return " << SS.getRFPObserved(i) << ", but returns" << GeneSS->getRFPObserved(i) << "\n";
            error = 1;
            globalError = 1;
        }
    }

    //This fails because this returns pointers to vectors and they need to be compared differently.
    std::vector <unsigned> *SSvec;
    std::vector <unsigned> *Gvec;
    for (unsigned i = 0; i < 64; i++)
    {
        SSvec = SS.getCodonPositions(i);
        Gvec = GeneSS->getCodonPositions(i);
        if (SSvec->size() != Gvec->size())
        {
            std::cerr << "Error with getSequenceSummary. Codon positions are incorrect.\n";
            std::cerr << "Information in compared vectors are not of equal size.\n";
            error = 1;
            globalError = 1;
        }
        else
        {
            for (unsigned j = 0; j < SSvec->size(); j++)
            {
                if (SSvec->at(j) != Gvec->at(j))
                {
                    std::cerr << "Error with getSequenceSummary. Codon positions are incorrect";
                    std::cerr << " for codon " << i << ".\n";
                    std::cerr << "Should return " << SSvec->at(j) << ", but returns" << Gvec->at(j) << "\n";
                    error = 1;
                    globalError = 1;
                }
            }
        }
    }

    unsigned AAListSize = (unsigned)SequenceSummary::aminoAcids().size();
    for (unsigned i = 0; i < AAListSize; i++)
    {
        if (SS.getAACountForAA(i) != GeneSS->getAACountForAA(i))
        {
            std::cerr << "Error with getSequenceSummary. AA counts are incorrect";
            std::cerr << " for amino acid " << i << ".\n";
            std::cerr << "Should return " << SS.getAACountForAA(i) << ", but returns" << GeneSS->getAACountForAA(i) << "\n";
            error = 1;
            globalError = 1;
        }
    }
    if (!error)
        std::cout << "Gene getSequenceSummary --- Pass\n";
    else
        error = 0; //Reset for next function.

    //---------------------------------------------------------//
    //------ get/setObservedSynthesisRateValues Function ------//
    //---------------------------------------------------------//
    std::vector <double> tmp;
    tmp = testGene.getObservedSynthesisRateValues();

    if (0 != tmp.size())
    {
        std::cerr << "Error with getObservedSynthesisRateValues. Function should return an empty vector but returns:\n";
        for (unsigned i = 0; i < tmp.size(); i++)
        {
            std::cerr << tmp[i] << "\n";
        }
        error = 1;
        globalError = 1;
    }

    tmp = {2.34, 3.234, 0.123};
    testGene.setObservedSynthesisRateValues(tmp);

    if (testGene.getObservedSynthesisRateValues() != tmp)
    {
        std::cerr << "Error in getObservedSynthesisRateValues or setObservedSynthesisRateValues. Function should return 2.34, 3.234, 0.123, but returns:\n";
        for (unsigned i = 0; i < tmp.size(); i++)
        {
            std::cerr << tmp[i] << "\n";
        }
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Gene get/setObservedSynthesisRateValues --- Pass\n";
    else
        error = 0; //Reset for next function.

    //--------------------------------------------------//
    //------ getNumObservedSynthesisSets Function ------//
    //--------------------------------------------------//
    if (3 != testGene.getNumObservedSynthesisSets())
    {
        std::cerr << "Error with getNumObservedSynthesisSets. Function should return 3, but returns ";
        std::cerr << testGene.getNumObservedSynthesisSets() << ".\n";
        globalError = 1;
    }
    else
        std::cout << "Gene getNumObservedSynthesisSets --- Pass\n";

    //-----------------------------------------------//
    //------ getObservedSynthesisRate Function ------//
    //-----------------------------------------------//

    // Declared above: tmp = {2.34, 3.234, 0.123}

    for (unsigned i = 0; i < 3; i++)
    {
        if (testGene.getObservedSynthesisRate(i) != tmp[i])
        {
            std::cerr << "Error with getObservedSynthesisRate. Function should return " << tmp[i] << "at index";
            std::cerr << i << ", but returns " << testGene.getObservedSynthesisRate(i) << ".\n";
            error = 1;
            globalError = 1;
        }
    }

    if (!error)
        std::cout << "Gene getObservedSynthesisRate --- Pass\n";
    else
        error = 0; //Reset for next function.

    //--------------------------------------//
    //------ getNucleotideAt Function ------//
    //--------------------------------------//
    if ('A' != testGene.getNucleotideAt(0))
    {
        std::cerr << "Error with getNucleotideAt. At index 0, the return value should be 'A', but is ";
        std::cerr << testGene.getNucleotideAt(0) << "\n";
        error = 1;
        globalError = 1;
    }
    if ('T' != testGene.getNucleotideAt(1))
    {
        std::cerr << "Error with getNucleotideAt. At index 1, the return value should be 'T', but is ";
        std::cerr << testGene.getNucleotideAt(1) << "\n";
        error = 1;
        globalError = 1;
    }
    if ('G' != testGene.getNucleotideAt(2))
    {
        std::cerr << "Error with getNucleotideAt. At index 2, the return value should be 'G', but is ";
        std::cerr << testGene.getNucleotideAt(2) << "\n";
        error = 1;
        globalError = 1;
    }
    if ('C' != testGene.getNucleotideAt(3))
    {
        std::cerr << "Error with getNucleotideAt. At index 3, the return value should be 'C', but is ";
        std::cerr << testGene.getNucleotideAt(3) << "\n";
        error = 1;
        globalError = 1;
    }
    if ('T' != testGene.getNucleotideAt(10))
    {
        std::cerr << "Error with getNucleotideAt. At index 10, the return value should be 'T', but is ";
        std::cerr << testGene.getNucleotideAt(10) << "\n";
        error = 1;
        globalError = 1;
    }
    if ('G' != testGene.getNucleotideAt(23))
    {
        std::cerr << "Error with getNucleotideAt. At index 23, the return value should be 'G', but is ";
        std::cerr << testGene.getNucleotideAt(23) << "\n";
        error = 1;
        globalError = 1;
    }
    if ('G' != testGene.getNucleotideAt(26))
    {
        std::cerr << "Error with getNucleotideAt. At index 26, the return value should be 'G', but is ";
        std::cerr << testGene.getNucleotideAt(26) << "\n";
        error = 1;
        globalError = 1;
    }
    if (!error)
        std::cout << "Gene getNucleotideAt --- Pass\n";
    else
        error = 0; //Reset for next function.

    //Todo: consider out of range test case?

    //-----------------------------//
    //------ length Function ------//
    //-----------------------------//
    if (testGene.length() == strlen("ATGCTCATTCTCACTGCTGCCTCGTAG"))
        std::cout << "Gene length --- Pass\n";
    else
    {
        std::cerr << "Error with length. Should return ";
        std::cerr << strlen("ATGCTCATTCTCACTGCTGCCTCGTAG") << " but returns: ";
        std::cerr << testGene.length() << "\n";
        globalError = 1;
    }

    //----------------------------------------//
    //------ reverseComplement Function ------//
    //----------------------------------------//
    Gene tmpGene;
    tmpGene = testGene.reverseComplement();
    if ("CTACGAGGCAGCAGTGAGAATGAGCAT" == tmpGene.getSequence())
        std::cout << "Gene reverseComplement --- Pass\n";
    else
    {
        std::cerr << "Error with reverseComplement. Should return \"CTACGAGGCAGCAGTGAGAATGAGCAT\" but returns: ";
        std::cerr << tmpGene.getSequence() << "\n";
        globalError = 1;
    }

    //-----------------------------------//
    //------ toAASequence Function ------//
    //-----------------------------------//
    if ("MLILTAASX" == testGene.toAASequence())
        std::cout << "Gene toAASequence --- Pass\n";
    else
    {
        std::cerr << "Error with toAASequence. Should return \"MLILTAASX\", but returns:" << testGene.toAASequence() <<
        "\n";
        globalError = 1;
    }

    //----------------------------//
    //------ clear Function ------//
    //----------------------------//
    testGene.clear();
    if ("" != testGene.getId())
    {
        std::cerr << "Error with clear. Gene Id should be blank, but is " << testGene.getId() << ".\n";
        error = 1;
        globalError = 1;
    }
    if ("" != testGene.getDescription())
    {
        std::cerr << "Error with clear. Gene description should be blank, but is " << testGene.getDescription() << ".\n";
        error = 1;
        globalError = 1;
    }
    if ("" != testGene.getSequence())
    {
        std::cerr << "Error with clear. Gene sequence should be blank, but is " << testGene.getSequence() << ".\n";
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Gene clear --- Pass\n";
    // No need to reset error

    return globalError;
}


int testGenome(std::string testFileDir)
{
    Genome genome;
    Genome testGenome;
    Gene g1("ATGGCCACTATTGGGTCTTAG", "TEST001", "TEST001 Test Gene");
    Gene g2("TGGGATTACCAA", "TEST002", "TEST002 Test Gene");
    Gene g3("TTGGAAACCACA", "TEST003", "TEST003 Test Gene");
    Gene g4("TGGGATTACCCC", "TEST004", "TEST004 Test Gene");
    Gene s1("TGGGATTACCAA", "TEST011", "TEST011 Test Gene"); //simulated gene
    int error = 0;
    int globalError = 0;

    /* Section 1:
     * Testing / Gene / Other Functions:
     * addGene, getGene, getGenes,
     * setNumGenesWithPhi, getNumGenesWithPhi, getNumGenesWithPhiForIndex,
     * getGenomeSize, getCodonCountsPerGene, clear
    */

    //TODO: should improper input be given (bad id/index)?

    //-----------------------------------//
    //------ get/addGene Functions ------//
    //-----------------------------------//
    genome.addGene(g1, false);
    genome.addGene(s1, true); //add the simulated gene s1

    Gene test = genome.getGene("TEST001", false);
    Gene test2 = genome.getGene(0, false);
    Gene test3 = genome.getGene("TEST011", true);
    Gene test4 = genome.getGene(0, true);

    if (!(test == g1 && test2 == g1)) //checking both by string and index
    {
        std::cerr << "Error in addGene or getGene with genes.\n";
        error = 1;
        globalError = 1;
    }

    if (!(test3 == s1 && test4 == s1)) //checking both by string and index
    {
        std::cerr << "Error in addGene or getGene with simulated genes.\n";
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Genome get/addGene --- Pass\n";
    else
        error = 0; //Reset for next function.

    //-------------------------------//
    //------ getGenes Function ------//
    //-------------------------------//
    std::vector<Gene> testVec;
    testVec.push_back(g1);

    if (!(testVec == genome.getGenes(false))) {
        std::cerr << "Error in getGenes(false).\n";
        error = 1;
        globalError = 1;
    }

    testVec.clear();
    testVec.push_back(s1);

    if (!(testVec == genome.getGenes(true))) {
        std::cerr << "Error in getGenes(true).\n";
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Genome getGenes --- Pass\n";
    else
        error = 0; //Reset for next function.

    //----------------------------------------------//
    //------ get/setNumGenesWithPhi Functions ------//
    //----------------------------------------------//
    genome.setNumGenesWithPhi({0, 1, 2, 3});

    std::vector<unsigned> uVector = {0, 1, 2, 3};

    if (genome.getNumGenesWithPhi() == uVector)
        std::cout << "Genome get/setNumGenesWithPhi --- Pass\n";
    else
    {
        std::cerr << "Error in setNumGenesWithPhi or getNumGenesWithPhi.\n";
        globalError = 1;
    }

    //-------------------------------------------------//
    //------ getNumGenesWithPhiForIndex Function ------//
    //-------------------------------------------------//
    for (unsigned i = 1; i < 4; i++) {
        if (genome.getNumGenesWithPhiForIndex(i) != i) {
            std::cerr << "Error in getNumGenesWithPhiForIndex with index ";
            std::cerr << i << ". Should return " << i << ", returns ";
            std::cerr << genome.getNumGenesWithPhiForIndex(i) << ".\n";
            error = 1;
            globalError = 1;
        }
    }

    if (!error)
        std::cout << "Genome getNumGenesWithPhiForIndex --- Pass\n";
    else
        error = 0; //Reset for next function.

    //------------------------------------//
    //------ getGenomeSize Function ------//
    //------------------------------------//
    if (1 != genome.getGenomeSize(false)) {
        std::cerr << "Error in getGenomesize(false). Should return 1, returns ";
        std::cerr << genome.getGenomeSize(false) << ".\n";
        error = 1;
        globalError = 1;
    }

    if (1 != genome.getGenomeSize(true)) {
        std::cerr << "Error in getGenomesize(true). Should return 1, returns ";
        std::cerr << genome.getGenomeSize(true) << ".\n";
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Genome getGenomeSize --- Pass\n";
    else
        error = 0; //Reset for next function.

    //--------------------------------------------//
    //------ getCodonCountsPerGene Function ------//
    //--------------------------------------------//

    //reuse generic vector of unsigned integers
    uVector = {1};
    if (uVector != genome.getCodonCountsPerGene("ATG")) {
        std::cerr << "Error in getCodonCountsPerGene with a single gene.\n";
        error = 1;
        globalError = 1;
    }

    genome.addGene(g2);
    genome.addGene(g4);

    uVector = {0, 1, 1};

    if (uVector != genome.getCodonCountsPerGene("GAT")) {
        std::cerr << "Error in getCodonCountsPerGene with three genes.\n";
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Genome getCodonCountsPerGene --- Pass\n";
    else
        error = 0; //Reset for next function.

    //----------------------------//
    //------ Clear Function ------//
    //----------------------------//

    // Empty Genome as a control variable
    Genome empty;

    // Test adding ObservedSynthesisRateValues
    Gene clear1("TTGATCGGGCAT", "TEST005", "TEST005 Test Gene");
    clear1.setObservedSynthesisRateValues({1, 2, 3, 4});
    genome.addGene(clear1, false);

    genome.clear();

    if (genome == empty)
        std::cout << "Genome clear --- Pass\n";
    else
    {
        std::cerr << "Error in clear. Genome is not empty.\n";
        globalError = 1;
    }

    /* Section 2:
     * Other and File I/O Functions:
     * getGenomeForGeneIndices
     * readFasta
     * readPANSEFile
     * readObservedPhiValues
    */

    //-----------------------------------------------//
    //------ getGenomeForGeneIndices Function ------//
    //-----------------------------------------------//

    // add more simulated and non-simulated genes
    genome.addGene(g1, false);
    genome.addGene(g2, false);
    genome.addGene(g3, false);
    genome.addGene(g4, false);

    //reuse generic vector of unsigned integers
    uVector = {0, 1, 2, 3};

    if (!(genome == genome.getGenomeForGeneIndices(uVector, false))) {
        std::cerr << "Error in getGenomeForGeneIndices with genes.\n";
        error = 1;
        globalError = 1;
    }

    genome.clear();

    Gene s2("TAGCATGATCCA", "TEST012", "TEST002 Test Gene"); //simulated gene
    Gene s3("TCATCAGGATTC", "TEST013", "TEST002 Test Gene"); //simulated gene
    Gene s4("AAACATGTCACG", "TEST014", "TEST002 Test Gene"); //simulated gene

    genome.addGene(s1, true);
    genome.addGene(s2, true);
    genome.addGene(s3, true);
    genome.addGene(s4, true);

    if (!(genome == genome.getGenomeForGeneIndices(uVector, true))) {
        std::cerr << "Error in getGenomeForGeneIndices with simulated genes.\n";
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Genome getGenomeForGeneIndices --- Pass\n";
    else
        error = 0; //Reset for next function.

    //--------------------------------//
    //------ readFasta Function ------//
    //--------------------------------//
    genome.clear();
    std::string file = testFileDir + "/" + "test.fasta";
    genome.readFasta(file, false);

    Gene fasta1("ATGACCGTAATTTTTTACTAG", "TEST002", "TEST002 Test Gene");
    Gene fasta2("ATGGTCTACTTTCTGACATAG", "TEST003", "TEST003 Test Gene");

    testGenome.addGene(g1, false);
    testGenome.addGene(fasta1, false);
    testGenome.addGene(fasta2, false);

    if (genome == testGenome)
        std::cout << "Genome readFasta --- Pass\n";
    else
    {
        std::cerr << "Error in readFasta. Genomes are not equivalent.\n";
        globalError = 1;
    }

    //---------------------------------//
    //------ writeFasta Function ------//
    //---------------------------------//

    // Now write a genome described above in readFasta to a file, read it in
    // again, and then compare its validity again.
    testGenome.clear();

    file = testFileDir + "/" + "testWrite.fasta";
    genome.writeFasta(file, false);
    testGenome.readFasta(file, false);

    if (!(genome == testGenome))
    {
        std::cerr << "Error in writeFasta with genes. Genomes are not equivalent.\n";
        error = 1;
        globalError = 1;
    }

    // Now, re-do writing check but with simulated genes.
    testGenome.clear();
    genome.clear();

    genome.addGene(g1, true);
    genome.addGene(fasta1, true);
    genome.addGene(fasta2, true);

    genome.writeFasta(file, true);

    // Note that while these genes were originally simulated, they are printed
    // as non-simulated genes.
    // It is up to the user to know that they were simulated, but they will
    // now be read in as non-simulated genes (and Unit Testing will compare their validity as such)

    testGenome.readFasta(file, true);

    genome.clear();
    genome.addGene(g1, false);
    genome.addGene(fasta1, false);
    genome.addGene(fasta2, false);

    if (!(genome == testGenome))
    {
        std::cerr << "Error in writeFasta with simulated genes. Genomes are not equivalent.\n";
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Genome writeFasta --- Pass\n";
    else
        error = 0; //Reset for next function. 

    //----------------------------------//
    //------ readRFPFile Function ------//
    //----------------------------------//
    genome.clear();
    testGenome.clear();

    file = testFileDir + "/" + "testReadRFP.csv";
    genome.readRFPFile(file);

    // These sequences are composed of those codons with RFPObserved values > 0
    // Each repetition of a codon denotes an incrementation in RFPObserved for that codon

    Gene rfp1("GCCGCCGCCGCCGCC", "TEST001", "No description for RFP Model");
    Gene rfp2("GCGGCGTTTTTTTTTTTT", "TEST002", "No description for RFP Model");
    Gene rfp3("ATGATGATGATGATGATGATGATGATGATGATGATGATG", "TEST003", "No description for RFP Model");

    testGenome.addGene(rfp1, false);
    testGenome.addGene(rfp2, false);
    testGenome.addGene(rfp3, false);

    if (genome == testGenome)
        std::cout << "Genome readRFPFile --- Pass\n";
    else
    {
        std::cerr << "Error in readRFPFile. Genomes are not equivalent.\n";
        globalError = 1;
    }

    //-----------------------------------//
    //------ writeRFPFile Function ------//
    //-----------------------------------//

    // Now write a genome described above in readRFPFile to a file, read it in
    // again, and then compare its validity again.
    testGenome.clear();

    file = testFileDir + "/" + "testWriteRFP.csv";
    genome.writeRFPFile(file, false);
    testGenome.readRFPFile(file);

    if (!(genome == testGenome))
    {
        std::cerr << "Error in writeRFPFile with genes. Genomes are not equivalent.\n";
        error = 1;
        globalError = 1;
    }

    // Now re-do writing check but with simulated genes.
    testGenome.clear();
    genome.clear();

    genome.addGene(rfp1, true);
    genome.addGene(rfp2, true);
    genome.addGene(rfp3, true);

    genome.writeRFPFile(file, true);

    // Note that while these genes were originally simulated, they are printed
    // as non-simulated genes.
    // It is up to the user to know that they were simulated, but they will
    // now be read in as non-simulated genes (and Unit Testing will compare their validity as such)

    genome.clear();
    genome.addGene(rfp1, false);
    genome.addGene(rfp2, false);
    genome.addGene(rfp3, false);

    testGenome.readRFPFile(file);

    if (!(genome == testGenome))
    {
        std::cerr << "Error in writeRFPFile with simulated genes. Genomes are not equivalent.\n";
        error = 1;
        globalError = 1;
    }

    if (!error)
        std::cout << "Genome writeRFPFile --- Pass\n";
    else
        error = 0; //Reset for next function.

    //------------------------------------//
    //------ readPANSEFile Function ------//
    //------------------------------------//
    genome.clear();
    testGenome.clear();

    file = testFileDir + "/" + "readPANSE.csv";
    genome.readPANSEFile(file);

    Gene panse1("CTTGCTATTTTT", "TEST001", "No description for PANSE Model");
    Gene panse2("CCTGTAATTTGG", "TEST002", "No description for PANSE Model");

    std::vector <unsigned> tmp1 = {0, 2, 0, 0};
    std::vector <unsigned> tmp2 = {0, 0, 1, 1};

    panse1.addRFP_count(tmp1);
    panse2.addRFP_count(tmp2);

    testGenome.addGene(panse1, false);
    testGenome.addGene(panse2, false);

    if (genome == testGenome)
        std::cout << "Genome readPANSE --- Pass\n";
    else
    {
        std::cerr << "Error in readPANSE. Genomes are not equivalent.\n";
        globalError = 1;
    }

    //-------------------------------------//
    //------ writePANSEFile Function ------//
    //-------------------------------------//

    /* readObservedPhiValues Testing Function
     *
     * Compares a genome with the readObservedPhiValues function's created genome.
     * Reads in "readObservedPhiValues.csv" and "readObservedPhiValuesError.csv" twice each,
     * once for byID and once for byIndex.
     *
     * Significant standard error output is produced by design: both files exhibit some errors.
    */
    //--------------------------------------------//
    //------ readObservedPhiValues Function ------//
    //--------------------------------------------//
    genome.clear();
    testGenome.clear();
    std::vector <double> emptyVec; // Empty vector used to clear ObservedSynthesisRateValues

    // Test 1: Test non-error file vs by ID readObservedPhiValues function
    genome.addGene(g1, false);
    g1.setObservedSynthesisRateValues({1, 2, 3, 4});
    testGenome.addGene(g1, false);

    genome.addGene(g2, false);
    g2.setObservedSynthesisRateValues({4, 3, 2, 1});
    testGenome.addGene(g2, false);

    genome.addGene(g3, false);
    g3.setObservedSynthesisRateValues({-1, -1, 4, 2});
    testGenome.addGene(g3, false);

    genome.addGene(g4, false);
    g4.setObservedSynthesisRateValues({2, 1, 4, -1});
    testGenome.addGene(g4, false);

    testGenome.setNumGenesWithPhi({3, 3, 4, 3});

    file = testFileDir + "/" + "readObservedPhiValues.csv";
    genome.readObservedPhiValues(file, true);

    if (!(genome == testGenome))
    {
        std::cerr << "Error comparing genomes: readObservedPhiValues.csv ";
        std::cerr << "by ID produces a different genome than expected.\n";
        error = 1;
        globalError = 1;
    }
    genome.clear();

    // Test 2: Test non-error file vs by index readObservedPhiValues function
    // Re-input genome as it was in the previous test, then run it by index instead
    g1.setObservedSynthesisRateValues(emptyVec);
    genome.addGene(g1, false);
    g2.setObservedSynthesisRateValues(emptyVec);
    genome.addGene(g2, false);
    g3.setObservedSynthesisRateValues(emptyVec);
    genome.addGene(g3, false);
    g4.setObservedSynthesisRateValues(emptyVec);
    genome.addGene(g4, false);

    genome.readObservedPhiValues(file, false);

    if (!(genome == testGenome))
    {
        std::cerr << "Error comparing genomes: readObservedPhiValues.csv ";
        std::cerr << "by index produces a different genome than expected.\n";
        error = 1;
        globalError = 1;
    }
    genome.clear();
    testGenome.clear();

    // Test 3: Test error file vs by ID readObservedPhiValues function
    // Since this file has an error in number of phi values, the ObservedSynthesisRateValues are cleared
    genome.addGene(g1, false);
    testGenome.addGene(g1, false);

    genome.addGene(g2, false);
    testGenome.addGene(g2, false);

    genome.addGene(g3, false);
    testGenome.addGene(g3, false);

    genome.addGene(g4, false);
    testGenome.addGene(g4, false);

    // As discussed in the documentation, however, NumGenesWithPhi is still initialized with 0's despite the error
    testGenome.setNumGenesWithPhi({0, 0, 0, 0});

    file = testFileDir + "/" + "readObservedPhiValuesError.csv";
    genome.readObservedPhiValues(file, true);

    if (!(genome == testGenome))
    {
        std::cerr << "Error comparing genomes: readObservedPhiValuesError.csv ";
        std::cerr << "by ID produces a different genome than expected.\n";
        error = 1;
        globalError = 1;
    }

    genome.clear();

    // Test 4: Test error file vs by index readObservedPhiValues function
    // Re-input genome as it was in the previous test, then run it by index instead
    genome.addGene(g1, false);
    genome.addGene(g2, false);
    genome.addGene(g3, false);
    genome.addGene(g4, false);

    genome.readObservedPhiValues(file, false);

    if (!(genome == testGenome))
    {
        std::cerr << "Error comparing genomes: readObservedPhiValuesError.csv ";
        std::cerr << "by index produces a different genome than expected.\n";
        error = 1;
        globalError = 1;
    }

    // If any errors are produced, reset variable for next function
    if (!error)
        std::cout << "Genome readObservedPhiValues --- Pass\n";
    // No need to reset error

    return globalError;
}


int testUtility()
{
    int error = 0;
    int globalError = 0;

    error = my_print("Product: %, Qty: %, Price is %", "Shampoo", 5, 100);

    if (error != 0)
    {
        std::cerr << "Error in my_print.\n";
        globalError = 1;
    }

    std::cout <<"\n";
    error = my_printError("Product: %, Qty: %, Price is %", "Shampoo", 5, 100);

    if (error != 0)
    {
        std::cerr << "Error in my_printError\n";
        globalError = 1;
    }

    return globalError;
}


int testCovarianceMatrix()
{
    CovarianceMatrix covM; //Default constructor sets numVariates to 2.
    int globalError = 0;

    //----------------------------------------------------------//
    //------ getCovMatrix & initCovarianceMatrix Function ------//
    //----------------------------------------------------------//

    // Currently, we can assume initCovarianceMatrix works since it is used in the default constructor.
    // TODO: Change this.

    std::vector <double> covM2 = {0.0025, 0, 0, 0, \
                                  0, 0.0025, 0, 0, \
                                  0, 0, 0.0025, 0, \
                                  0, 0, 0, 0.0025};
    std::vector <double> *covM2star = &covM2;

    covM.initCovarianceMatrix(4);

    // TODO: Should it really return a pointer?
    /*
    if (!(covM.getCovMatrix() == covM2star)) {
        std::cerr << "Error in getCovMatrix or initCovarianceMatrix.\n";
        globalError = 1;
    }
    else
        std::cout << "CovarianceMatrix getCovMatrix & initCovarianceMatrix --- Pass\n";
    */

    //------------------------------//
    //------ setDiag Function ------//
    //------------------------------//
    covM.setDiag(3.14);

    covM2 = {3.14, 0, 0, 0, \
             0, 3.14, 0, 0, \
             0, 0, 3.14, 0, \
             0, 0, 0, 3.14};

    /*
    if (!(covM.getCovMatrix() == covM2star)) {
        std::cerr << "Error in getCovMatrix or initCovarianceMatrix.\n";
        globalError = 1;
    }
    else
        std::cout << "CovarianceMatrix getCovMatrix & initCovarianceMatrix --- Pass\n";
    */

    // TODO: Write function to extract choleski matrix similar to covMatrix above
    //--------------------------------------------//
    //------ choleskiDecomposition Function ------//
    //--------------------------------------------//

    // TODO: Test print functions somehow
    //--------------------------------------------//
    //------ printCovarianceMatrix Function ------//
    //--------------------------------------------//
    //covM.printCovarianceMatrix();

    //------------------------------------------//
    //------ printCholeskiMatrix Function ------//
    //------------------------------------------//
    //covM.printCholeskiMatrix();

    //-------------------------------------//
    //------ getNumVariates Function ------//
    //-------------------------------------//
    covM.getNumVariates();

    if (covM.getNumVariates() != 4) {
        std::cerr << "Error in getNumVariates. Function should return 4, but returns " << covM.getNumVariates() << ".\n";
        globalError = 1;
    }
    else
        std::cout << "CovarianceMatrix getNumVariates --- Pass\n";

    //-------------------------------------------------------------//
    //------ transformIidNumersIntoCovaryingNumbers Function ------//
    //-------------------------------------------------------------//

    //------------------------------------------------//
    //------ calculateSampleCovariance Function ------//
    //------------------------------------------------//

    return globalError;
}


// -----------------------------------------------------------------------------------------------------//
// ---------------------------------------- R SECTION --------------------------------------------------//
// -----------------------------------------------------------------------------------------------------//

#ifndef STANDALONE
//---------------------------------//
//---------- RCPP Module ----------//
//---------------------------------//

RCPP_MODULE(Test_mod)
{
	function("testSequenceSummary", &testSequenceSummary);
	function("testGene", &testGene);
	function("testGenome", &testGenome);
	function("testUtility", &testUtility);
	function("testCovarianceMatrix", &testCovarianceMatrix);
}
#endif
