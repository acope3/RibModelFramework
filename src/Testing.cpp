#include "include/Testing.h"


void testSequenceSummary()
{
    SequenceSummary SS("ATGCTCATTCTCACTGCTGCCTCGTAG");
    int error = 0;
    //----------------------------//
    //------ Clear Function ------//
    //----------------------------//
    SS.clear();
    for (unsigned i = 0; i < 64; i++)
    {
        if (0 != SS.getCodonCountForCodon(i))
        {
            std::cerr <<"Problem with Sequence Summary \"clear\" function.\n";
            std::cerr <<"Problem at codon index" << i << "\n";
            error = 1;
        }
    }
    for (unsigned i = 0; i < 22; i++)
    {
        if (0 != SS.getAACountForAA(i))
        {
            std::cerr <<"Problem with Sequence Summary \"clear\" function.\n";
            std::cerr <<"Problem at amino acid index" << i << "\n";
            error = 1;
        }
    }

    if (!error)
    {
        std::cout <<"Sequence Summary clear --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }


    //--------------------------------------//
    //------ ProcessSequence Function ------//
    //--------------------------------------//
    SS.processSequence("ATGCTCATTCTCACTGCTGCCTCGTAG");

    if (1 != SS.getAACountForAA("I"))
    {
        std::cerr <<"Problem with Sequence Summary \"processSequence\" function.\n";
        std::cerr <<"Problem with amino acid \"I\".";
        std::cerr <<"I is in the sequence once, but is returning " << SS.getAACountForAA("I") << "\n";
        error = 1;
    }

    if (1 != SS.getAACountForAA("T"))
    {
        std::cerr <<"Problem with Sequence Summary \"processSequence\" function.\n";
        std::cerr <<"Problem with amino acid \"T\".";
        std::cerr <<"T is in the sequence once, but is returning " << SS.getAACountForAA("T") << "\n";
        error = 1;
    }

    std::string codon = "ATT";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        std::cerr <<"Problem with Sequence Summary \"processSequence\" function.\n";
        std::cerr <<"Problem with codon \"ATT\".";
        std::cerr <<"ATT is in the sequence once, but is returning " << SS.getCodonCountForCodon(codon) << "\n";
        error = 1;
    }

    codon = "ACT";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        std::cerr <<"Problem with Sequence Summary \"processSequence\" function.\n";
        std::cerr <<"Problem with codon \"ACT\".";
        std::cerr <<"ACT is in the sequence once, but is returning " << SS.getCodonCountForCodon(codon) << "\n";
        error = 1;
    }

    codon = "GCT";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        std::cerr <<"Problem with Sequence Summary \"processSequence\" function.\n";
        std::cerr <<"Problem with codon \"GCT\".";
        std::cerr <<"GCT is in the sequence once, but is returning " << SS.getCodonCountForCodon(codon) << "\n";
        error = 1;
    }

    codon = "GCC";
    if (1 != SS.getCodonCountForCodon(codon))
    {
        std::cerr <<"Problem with Sequence Summary \"processSequence\" function.\n";
        std::cerr <<"Problem with codon \"GCC\".";
        std::cerr <<"GCC is in the sequence once, but is returning " << SS.getCodonCountForCodon(codon) << "\n";
        error = 1;
    }

    std::vector <unsigned> *tmp;
    tmp = SS.getCodonPositions("CTC");
    if ((1 != tmp -> at(0)) && (3 != tmp -> at(1)))
    {
        std::cerr <<"Codon CTC should be found at position 1 and 3(zero indexed), but is";
        std::cerr <<"found at these locations:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    tmp = SS.getCodonPositions("ATT");
    if (2 != tmp -> at(0))
    {
        std::cerr <<"Codon ATT should be found at position 2(zero indexed), but is";
        std::cerr <<"found at these locations:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    if (!error)
    {
        std::cout <<"Sequence Summary processSequence --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }


    //------------------------------------------//
    //------ complimentNucleotide Function------//
    //------------------------------------------//

    if ('T' != SequenceSummary::complimentNucleotide('A'))
    {
        std::cerr <<"The compliment of A should be T\n";
        error = 1;
    }

    if ('A' != SequenceSummary::complimentNucleotide('T'))
    {
        std::cerr <<"The compliment of T should be A\n";
        error = 1;
    }

    if ('G' != SequenceSummary::complimentNucleotide('C'))
    {
        std::cerr <<"The compliment of C should be G\n";
        error = 1;
    }

    if ('C' != SequenceSummary::complimentNucleotide('G'))
    {
        std::cerr <<"The compliment of G should be C\n";
        error = 1;
    }

    if ('C' != SequenceSummary::complimentNucleotide('Q'))
    {
        std::cerr <<"The compliment of Q should be C\n";
        error = 1;
    }

    if (!error)
    {
        std::cout <<"Sequence Summary complimentNucleotide --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }


    //----------------------------------------------//
    //------ getAACountForAA(string) Function ------//
    //----------------------------------------------//
    SequenceSummary SS2("ATGCTCATTCTCACTGCTGCCTCGTAG");

    if (1 != SS2.getAACountForAA("M"))
    {
        std::cerr <<"Error with getAACountForAA(string) for amino acid M.\n";
        std::cerr <<"Should return 1, returns " << SS2.getAACountForAA("M") <<"\n";
        error = 1;
    }

    if (2 != SS2.getAACountForAA("L"))
    {
        std::cerr <<"Error with getAACountForAA(string) for amino acid L.\n";
        std::cerr <<"Should return 2, returns " << SS2.getAACountForAA("L") <<"\n";
        error = 1;
    }

    if (1 != SS2.getAACountForAA("I"))
    {
        std::cerr <<"Error with getAACountForAA(string) for amino acid I.\n";
        std::cerr <<"Should return 1, returns " << SS2.getAACountForAA("I") <<"\n";
        error = 1;
    }

    if (1 != SS2.getAACountForAA("T"))
    {
        std::cerr <<"Error with getAACountForAA(string) for amino acid T.\n";
        std::cerr <<"Should return 1, returns " << SS2.getAACountForAA("T") <<"\n";
        error = 1;
    }

    if (2 != SS2.getAACountForAA("A"))
    {
        std::cerr <<"Error with getAACountForAA(string) for amino acid A.\n";
        std::cerr <<"Should return 2, returns " << SS2.getAACountForAA("A") <<"\n";
        error = 1;
    }

    if (1 != SS2.getAACountForAA("S"))
    {
        std::cerr <<"Error with getAACountForAA(string) for amino acid S.\n";
        std::cerr <<"Should return 1, returns " << SS2.getAACountForAA("S") <<"\n";
        error = 1;
    }

    if (1 != SS2.getAACountForAA("X"))
    {
        std::cerr <<"Error with getAACountForAA(string) for amino acid X.\n";
        std::cerr <<"Should return 1, returns " << SS2.getAACountForAA("X") <<"\n";
        error = 1;
    }

    if (0 != SS2.getAACountForAA("G"))
    {
        std::cerr <<"Error with getAACountForAA(string) for amino acid G.\n";
        std::cerr <<"Should return 0, returns " << SS2.getAACountForAA("G") <<"\n";
        error = 1;
    }

    if (!error)
    {
        std::cout <<"Sequence Summary getAACountForAA(string) --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }

    //---------------------------------------------//
    //------ getAACountForAA(index) Function ------//
    //---------------------------------------------//
    if (1 != SS2.getAACountForAA(10))
    {
        std::cerr <<"Error with getAACountForAA(index) for amino acid M (index 10).\n";
        std::cerr <<"Should return 1, returns " << SS2.getAACountForAA(10) <<"\n";
        error = 1;
    }

    if (2 != SS2.getAACountForAA(9))
    {
        std::cerr <<"Error with getAACountForAA(index) for amino acid L (index 9).\n";
        std::cerr <<"Should return 2, returns " << SS2.getAACountForAA(9) <<"\n";
        error = 1;
    }

    if (1 != SS2.getAACountForAA(7))
    {
        std::cerr <<"Error with getAACountForAA(index) for amino acid I (index 7).\n";
        std::cerr <<"Should return 1, returns " << SS2.getAACountForAA(7) <<"\n";
        error = 1;
    }

    if (1 != SS2.getAACountForAA(16))
    {
        std::cerr <<"Error with getAACountForAA(index) for amino acid T (index 16).\n";
        std::cerr <<"Should return 1, returns " << SS2.getAACountForAA(16) <<"\n";
        error = 1;
    }

    if (2 != SS2.getAACountForAA(0))
    {
        std::cerr <<"Error with getAACountForAA(index) for amino acid A (index 0).\n";
        std::cerr <<"Should return 2, returns " << SS2.getAACountForAA(0) <<"\n";
        error = 1;
    }

    if (1 != SS2.getAACountForAA(15))
    {
        std::cerr <<"Error with getAACountForAA(index) for amino acid S (index 15).\n";
        std::cerr <<"Should return 1, returns " << SS2.getAACountForAA(15) <<"\n";
        error = 1;
    }

    if (1 != SS2.getAACountForAA(21))
    {
        std::cerr <<"Error with getAACountForAA(index) for amino acid X (index 21).\n";
        std::cerr <<"Should return 1, returns " << SS2.getAACountForAA(21) <<"\n";
        error = 1;
    }

    if (0 != SS2.getAACountForAA(2))
    {
        std::cerr <<"Error with getAACountForAA(index) for amino acid D (index 2).\n";
        std::cerr <<"Should return 0, returns " << SS2.getAACountForAA(2) <<"\n";
        error = 1;
    }

    if (!error)
    {
        std::cout <<"Sequence Summary getAACountForAA(index) --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }

    //--------------------------------------------//
    //------ getCodonCountsForCodon(string) ------//
    //--------------------------------------------//

    codon = "ATG";
    if (1 != SS2.getCodonCountForCodon(codon))
    {
        std::cerr <<"Error with getCodonCountForCodon(string) for " << codon <<".\n";
        std::cerr <<"Should return 1, but returns " << SS2.getCodonCountForCodon(codon) <<"\n";
        error = 1;
    }

    codon = "CTC";
    if (2 != SS2.getCodonCountForCodon(codon))
    {
        std::cerr <<"Error with getCodonCountForCodon(string) for " << codon <<".\n";
        std::cerr <<"Should return 2, but returns " << SS2.getCodonCountForCodon(codon) <<"\n";
        error = 1;
    }

    codon = "ATT";
    if (1 != SS2.getCodonCountForCodon(codon))
    {
        std::cerr <<"Error with getCodonCountForCodon(string) for " << codon <<".\n";
        std::cerr <<"Should return 1, but returns " << SS2.getCodonCountForCodon(codon) <<"\n";
        error = 1;
    }

    codon = "ACT";
    if (1 != SS2.getCodonCountForCodon(codon))
    {
        std::cerr <<"Error with getCodonCountForCodon(string) for " << codon <<".\n";
        std::cerr <<"Should return 1, but returns " << SS2.getCodonCountForCodon(codon) <<"\n";
        error = 1;
    }

    codon = "GCT";
    if (1 != SS2.getCodonCountForCodon(codon))
    {
        std::cerr <<"Error with getCodonCountForCodon(string) for " << codon <<".\n";
        std::cerr <<"Should return 1, but returns " << SS2.getCodonCountForCodon(codon) <<"\n";
        error = 1;
    }

    codon = "GCC";
    if (1 != SS2.getCodonCountForCodon(codon))
    {
        std::cerr <<"Error with getCodonCountForCodon(string) for " << codon <<".\n";
        std::cerr <<"Should return 1, but returns " << SS2.getCodonCountForCodon(codon) <<"\n";
        error = 1;
    }

    codon = "TCG";
    if (1 != SS2.getCodonCountForCodon(codon))
    {
        std::cerr <<"Error with getCodonCountForCodon(string) for " << codon <<".\n";
        std::cerr <<"Should return 1, but returns " << SS2.getCodonCountForCodon(codon) <<"\n";
        error = 1;
    }

    codon = "TAG";
    if (1 != SS2.getCodonCountForCodon(codon))
    {
        std::cerr <<"Error with getCodonCountForCodon(string) for " << codon <<".\n";
        std::cerr <<"Should return 1, but returns " << SS2.getCodonCountForCodon(codon) <<"\n";
        error = 1;
    }

    codon = "AAA";
    if (0 != SS2.getCodonCountForCodon(codon))
    {
        std::cerr <<"Error with getCodonCountForCodon(string) for " << codon <<".\n";
        std::cerr <<"Should return 0, but returns " << SS2.getCodonCountForCodon(codon) <<"\n";
        error = 1;
    }

    if (!error)
    {
        std::cout <<"Sequence Summary getCodonCountsForCodon(string) --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }

    //----------------------------------------------------//
    //------ getCodonCountsForCodon(index) Function ------//
    //----------------------------------------------------//

    if (1 != SS2.getCodonCountForCodon(29))
    {
        std::cerr <<"Error with getCodonCountForCodon(index) for codon \"ATG\" (index 29).\n";
        std::cerr <<"Should return 1, but returns " << SS2.getCodonCountForCodon(29) <<"\n";
        error = 1;
    }

    if (2 != SS2.getCodonCountForCodon(24))
    {
        std::cerr <<"Error with getCodonCountForCodon(index) for codon \"CTC\" (index 24).\n";
        std::cerr <<"Should return 2, but returns " << SS2.getCodonCountForCodon(24) <<"\n";
        error = 1;
    }

    if (1 != SS2.getCodonCountForCodon(20))
    {
        std::cerr <<"Error with getCodonCountForCodon(index) for codon \"ATT\" (index 20).\n";
        std::cerr <<"Should return 1, but returns " << SS2.getCodonCountForCodon(20) <<"\n";
        error = 1;
    }

    if (1 != SS2.getCodonCountForCodon(51))
    {
        std::cerr <<"Error with getCodonCountForCodon(index) for codon \"ACT\" (index 51).\n";
        std::cerr <<"Should return 1, but returns " << SS2.getCodonCountForCodon(51) <<"\n";
        error = 1;
    }

    if (1 != SS2.getCodonCountForCodon(3))
    {
        std::cerr <<"Error with getCodonCountForCodon(index) for codon \"GCT\" (index 3).\n";
        std::cerr <<"Should return 1, but returns " << SS2.getCodonCountForCodon(3) <<"\n";
        error = 1;
    }

    if (1 != SS2.getCodonCountForCodon(1))
    {
        std::cerr <<"Error with getCodonCountForCodon(index) for codon \"GCC\" (index 1).\n";
        std::cerr <<"Should return 1, but returns " << SS2.getCodonCountForCodon(1) <<"\n";
        error = 1;
    }

    if (1 != SS2.getCodonCountForCodon(46))
    {
        std::cerr <<"Error with getCodonCountForCodon(index) for codon \"TCG\" (index 46).\n";
        std::cerr <<"Should return 1, but returns " << SS2.getCodonCountForCodon(46) <<"\n";
        error = 1;
    }

    if (1 != SS2.getCodonCountForCodon(62))
    {
        std::cerr <<"Error with getCodonCountForCodon(index) for codon \"TAG\" (index 62).\n";
        std::cerr <<"Should return 1, but returns " << SS2.getCodonCountForCodon(62) <<"\n";
        error = 1;
    }

    if (0 != SS2.getCodonCountForCodon(2))
    {
        std::cerr <<"Error with getCodonCountForCodon(index) for codon \"AAA\" (index 2).\n";
        std::cerr <<"Should return 1, but returns " << SS2.getCodonCountForCodon(2) <<"\n";
        error = 1;
    }

    if (!error)
    {
        std::cout <<"Sequence Summary getCodonCountsForCodon(index) --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }

    //---------------------------------------------//
    //------ getRFPObserved(string) Function ------//
    //---------------------------------------------//

    SS2.setRFPObserved(4, 35);
    SS2.setRFPObserved(16,45);
    SS2.setRFPObserved(54,2);
    SS2.setRFPObserved(45,0);

    if (35 != SS2.getRFPObserved("TGC"))
    {
        std::cerr <<"Error with getRFPObserved(string) for codon \"TGC\".\n";
        std::cerr <<"should return 35, but returns " << SS2.getRFPObserved("TGC") <<"\n";
        error = 1;
    }

    if (45 != SS2.getRFPObserved("CAC"))
    {
        std::cerr <<"Error with getRFPObserved(string) for codon \"CAC\".\n";
        std::cerr <<"should return 45, but returns " << SS2.getRFPObserved("CAC") <<"\n";
        error = 1;
    }

    if (2 != SS2.getRFPObserved("GTG"))
    {
        std::cerr <<"Error with getRFPObserved(string) for codon \"GTG\".\n";
        std::cerr <<"should return 2, but returns " << SS2.getRFPObserved("GTG") <<"\n";
        error = 1;
    }

    if (0 != SS2.getRFPObserved("TCC"))
    {
        std::cerr <<"Error with getRFPObserved(string) for codon \"TCC\".\n";
        std::cerr <<"should return 0, but returns " << SS2.getRFPObserved("TCC") <<"\n";
        error = 1;
    }

    if (!error)
    {
        std::cout <<"Sequence Summary getRFPObserved(string) --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }

    //--------------------------------------------//
    //------ getRFPObserved(index) Function ------//
    //--------------------------------------------//

    SS2.setRFPObserved(0,45);
    SS2.setRFPObserved(1,52);
    SS2.setRFPObserved(2,63);
    SS2.setRFPObserved(60,23);

    if (45 != SS2.getRFPObserved(0))
    {
        std::cerr <<"Error with getRFPObserved(index) for codon index 0.\n";
        std::cerr <<"should return 45, but returns " << SS2.getRFPObserved(0) <<"\n";
        error = 1;
    }

    if (52 != SS2.getRFPObserved(1))
    {
        std::cerr <<"Error with getRFPObserved(index) for codon index 1.\n";
        std::cerr <<"should return 52, but returns " << SS2.getRFPObserved(1) <<"\n";
        error = 1;
    }

    if (63 != SS2.getRFPObserved(2))
    {
        std::cerr <<"Error with getRFPObserved(index) for codon index 2.\n";
        std::cerr <<"should return 63, but returns " << SS2.getRFPObserved(2) <<"\n";
        error = 1;
    }

    if (23 != SS2.getRFPObserved(60))
    {
        std::cerr <<"Error with getRFPObserved(index) for codon index 60.\n";
        std::cerr <<"should return 23, but returns " << SS2.getRFPObserved(60) <<"\n";
        error = 1;
    }

    if (!error)
    {
        std::cout <<"Sequence Summary getRFPObserved(index) --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }


    //------------------------------------------------//
    //------ getCodonPositions(string) Function ------//
    //------------------------------------------------//

    tmp = SS2.getCodonPositions("ATG");
    if (tmp -> at(0) != 0 || tmp -> size() != 1)
    {
        std::cerr <<"Error with getCodonPositions(string) for codon \"ATG\".\n";
        std::cerr <<"Should return 0, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    tmp = SS2.getCodonPositions("CTC");
    if (tmp -> at(0) != 1 || tmp -> at(1) != 3|| tmp -> size() != 2)
    {
        std::cerr <<"Error with getCodonPositions(string) for codon \"CTC\".\n";
        std::cerr <<"Should return 1 and 3, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    tmp = SS2.getCodonPositions("ATT");
    if (tmp -> at(0) != 2 || tmp -> size() != 1)
    {
        std::cerr <<"Error with getCodonPositions(string) for codon \"ATT\".\n";
        std::cerr <<"Should return 2, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    tmp = SS2.getCodonPositions("ACT");
    if (tmp -> at(0) != 4 || tmp -> size() != 1)
    {
        std::cerr <<"Error with getCodonPositions(string) for codon \"ACT\".\n";
        std::cerr <<"Should return 4, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }


    tmp = SS2.getCodonPositions("GCT");
    if (tmp -> at(0) != 5 || tmp -> size() != 1)
    {
        std::cerr <<"Error with getCodonPositions(string) for codon \"GCT\".\n";
        std::cerr <<"Should return 5, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    tmp = SS2.getCodonPositions("GCC");
    if (tmp -> at(0) != 6 || tmp -> size() != 1)
    {
        std::cerr <<"Error with getCodonPositions(string) for codon \"GCC\".\n";
        std::cerr <<"Should return 6, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    tmp = SS2.getCodonPositions("TCG");
    if (tmp -> at(0) != 7 || tmp -> size() != 1)
    {
        std::cerr <<"Error with getCodonPositions(string) for codon \"TCG\".\n";
        std::cerr <<"Should return 7, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }


    tmp = SS2.getCodonPositions("TAG");
    if (tmp -> at(0) != 8 || tmp -> size() != 1)
    {
        std::cerr <<"Error with getCodonPositions(string) for codon \"TAG\".\n";
        std::cerr <<"Should return 8, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    tmp = SS2.getCodonPositions("GTG");
    if (tmp -> size() != 0)
    {
        std::cerr <<"Error with getCodonPositions(string) for codon \"GTG\".\n";
        std::cerr <<"Should return an empty vector, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    if (!error)
    {
        std::cout <<"Sequence Summary getCodonPositions(string) --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }

    //-----------------------------------------------//
    //------ getCodonPositions(index) Function ------//
    //-----------------------------------------------//

    tmp = SS2.getCodonPositions(29);
    if (tmp -> at(0) != 0 || tmp -> size() != 1)
    {
        std::cerr <<"Error with getCodonPositions(index) for codon index 29.\n";
        std::cerr <<"Should return 0, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    tmp = SS2.getCodonPositions(24);
    if (tmp -> at(0) != 1 || tmp -> at(1) != 3 || tmp -> size() != 2)
    {
        std::cerr <<"Error with getCodonPositions(index) for codon index 24.\n";
        std::cerr <<"Should return 1 and 3, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    tmp = SS2.getCodonPositions(20);
    if (tmp -> at(0) != 2 || tmp -> size() != 1)
    {
        std::cerr <<"Error with getCodonPositions(index) for codon index 20.\n";
        std::cerr <<"Should return 2, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    tmp = SS2.getCodonPositions(51);
    if (tmp -> at(0) != 4 || tmp -> size() != 1)
    {
        std::cerr <<"Error with getCodonPositions(index) for codon index 51.\n";
        std::cerr <<"Should return 4, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    tmp = SS2.getCodonPositions(3);
    if (tmp -> at(0) != 5 || tmp -> size() != 1)
    {
        std::cerr <<"Error with getCodonPositions(index) for codon index 3.\n";
        std::cerr <<"Should return 4, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    tmp = SS2.getCodonPositions(1);
    if (tmp -> at(0) != 6 || tmp -> size() != 1)
    {
        std::cerr <<"Error with getCodonPositions(index) for codon index 1.\n";
        std::cerr <<"Should return 4, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    tmp = SS2.getCodonPositions(46);
    if (tmp -> at(0) != 7 || tmp -> size() != 1)
    {
        std::cerr <<"Error with getCodonPositions(index) for codon index 46.\n";
        std::cerr <<"Should return 7, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    tmp = SS2.getCodonPositions(62);
    if (tmp -> at(0) != 8 || tmp -> size() != 1)
    {
        std::cerr <<"Error with getCodonPositions(index) for codon index 62.\n";
        std::cerr <<"Should return 8, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    tmp = SS2.getCodonPositions(54);
    if (tmp -> size() != 0)
    {
        std::cerr <<"Error with getCodonPositions(index) for codon index 54.\n";
        std::cerr <<"Should return an empty vector, but returns:\n";
        for (unsigned i = 0; i < tmp -> size(); i++)
        {
            std::cerr << tmp -> at(i) <<"\n";
        }
        error = 1;
    }

    if (!error)
    {
        std::cout <<"Sequence Summary getCodonPositions(index) --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }
}




void testGene()
{
    int error = 0;

    //--------------------------------//
    //------ getSequenceSummary ------//
    //--------------------------------//

    Gene testGene("ATGCTCATTCTCACTGCTGCCTCGTAG", "1", "test Gene");
    SequenceSummary SS("ATGCTCATTCTCACTGCTGCCTCGTAG");
    SequenceSummary *GeneSS = testGene.getSequenceSummary();
    for (unsigned i = 0; i < 64; i++)
    {
        if (SS.getCodonCountForCodon(i) != GeneSS->getCodonCountForCodon(i))
        {
            std::cerr <<"Error with getSequenceSummary. Codon counts are incorrect";
            std::cerr <<" for codon " << i <<", " << SequenceSummary::codonArray[i] <<".\n";
            std::cerr <<"Should return " << SS.getCodonCountForCodon(i) <<", but returns" << GeneSS->getCodonCountForCodon(i) <<"\n";
            error = 1;
        }
    }


    for (unsigned i = 0; i < 64; i++)
    {
        if (SS.getRFPObserved(i) != GeneSS->getRFPObserved(i))
        {
            std::cerr <<"Error with getSequenceSummary. RFP observed is incorrect";
            std::cerr <<" for codon " << i <<".\n";
            std::cerr <<"Should return " << SS.getRFPObserved(i) <<", but returns" << GeneSS->getRFPObserved(i) <<"\n";
            error = 1;
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
            std::cerr <<"Error with getSequenceSummary. Codon positions are incorrect.\n";
            std::cerr <<"Information in compared vectors are not of equal size.\n";
            error = 1;
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
                }
            }
        }
    }


    unsigned AAListSize = (unsigned)SequenceSummary::aminoAcids().size();
    for (unsigned i = 0; i < AAListSize; i++)
    {
        if (SS.getAACountForAA(i) != GeneSS->getAACountForAA(i))
        {
            std::cerr <<"Error with getSequenceSummary. AA counts are incorrect";
            std::cerr <<" for amino acid " << i <<".\n";
            std::cerr <<"Should return " << SS.getAACountForAA(i) <<", but returns" << GeneSS->getAACountForAA(i) <<"\n";
            error = 1;
        }
    }
    if (!error)
    {
        std::cout <<"Gene getSequenceSummary --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }

    //-----------------------------------------------//
    //------ get/setObservedPhiValues Function ------//
    //-----------------------------------------------//
    std::vector <double> tmp;
    tmp = testGene.getObservedSynthesisRateValues();

    if (0 != tmp.size())
    {
        std::cerr <<"Error with getObservedPhiValues. Function should return an empty vector but returns:\n";
        for (unsigned i = 0; i < tmp.size(); i++)
        {
            std::cerr << tmp[i] <<"\n";
        }
        error = 1;
    }
    tmp.resize(3);
    tmp[0] = 2.34;
    tmp[1] = 3.234;
    tmp[2] = 0.123;
    testGene.setObservedSynthesisRateValues(tmp);
    tmp = testGene.getObservedSynthesisRateValues();
    if (3 != tmp.size() || 2.34 != tmp[0] || 3.234 != tmp[1] || 0.123 != tmp[2])
    {
        std::cerr <<"Error with getObservedPhiValues. Function should return 2.34, 3.234, 0.123, but returns:\n";
        for (unsigned i = 0; i < tmp.size(); i++)
        {
            std::cerr << tmp[i] <<"\n";
        }
        error = 1;
    }

    if (!error)
    {
        std::cout <<"Gene get/setObservedPhiValues --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }


    //--------------------------------------------------//
    //------ getNumObservedSynthesisSets Function ------//
    //--------------------------------------------------//

    if (3 != testGene.getNumObservedSynthesisSets())
    {
        std::cerr <<"Error with getNumObservedSynthesisSets. Function should return 3, but returns ";
        std::cerr << testGene.getNumObservedSynthesisSets() <<".\n";
    }
    else
    {
        std::cout <<"Gene getNumObservedSynthesisSets --- Pass\n";
    }




    //-----------------------------------------------//
    //------ getObservedSynthesisRate Function ------//
    //-----------------------------------------------//

    tmp = testGene.getObservedSynthesisRateValues();
    double trueValues[3] = {2.34, 3.234, 0.123};
    for (unsigned i = 0; i < 3; i++)
    {
        if (tmp[i] != trueValues[i])
        {
            std::cerr <<"Error with getObservedSynthesisRate. Function should return " << trueValues[i] <<"at index";
            std::cerr << i <<", but returns " << tmp[i] <<".\n";
            error = 1;
        }
    }

    if (!error)
    {
        std::cout <<"Gene getObservedSynthesisRate --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }




    //--------------------------------------//
    //------ getNucleotideAt Function ------//
    //--------------------------------------//


    if ('A' != testGene.getNucleotideAt(0))
    {
        std::cerr <<"Error with getNucleotideAt. At index 0, the return value should be 'A', but is ";
        std::cerr << testGene.getNucleotideAt(0) <<"\n";
        error = 1;
    }
    if ('T' != testGene.getNucleotideAt(1))
    {
        std::cerr <<"Error with getNucleotideAt. At index 1, the return value should be 'T', but is ";
        std::cerr << testGene.getNucleotideAt(1) <<"\n";
        error = 1;
    }
    if ('G' != testGene.getNucleotideAt(2))
    {
        std::cerr <<"Error with getNucleotideAt. At index 2, the return value should be 'G', but is ";
        std::cerr << testGene.getNucleotideAt(2) <<"\n";
        error = 1;
    }
    if ('C' != testGene.getNucleotideAt(3))
    {
        std::cerr <<"Error with getNucleotideAt. At index 3, the return value should be 'C', but is ";
        std::cerr << testGene.getNucleotideAt(3) <<"\n";
        error = 1;
    }
    if ('T' != testGene.getNucleotideAt(10))
    {
        std::cerr <<"Error with getNucleotideAt. At index 10, the return value should be 'T', but is ";
        std::cerr << testGene.getNucleotideAt(10) <<"\n";
        error = 1;
    }
    if ('G' != testGene.getNucleotideAt(23))
    {
        std::cerr <<"Error with getNucleotideAt. At index 23, the return value should be 'G', but is ";
        std::cerr << testGene.getNucleotideAt(23) <<"\n";
        error = 1;
    }
    if ('G' != testGene.getNucleotideAt(26))
    {
        std::cerr <<"Error with getNucleotideAt. At index 26, the return value should be 'G', but is ";
        std::cerr << testGene.getNucleotideAt(26) <<"\n";
        error = 1;
    }
    if (!error)
    {
        std::cout <<"Gene getNucleotideAt --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }

    //Todo: consider out of range test case?


    //-------------------------------//
    //------ reverseComplement ------//
    //-------------------------------//

    Gene tmpGene;
    tmpGene = testGene.reverseComplement();
    if ("CTACGAGGCAGCAGTGAGAATGAGCAT" != tmpGene.getSequence())
    {
        std::cerr <<"Error with reverseComplement. Should return \"CTACGAGGCAGCAGTGAGAATGAGCAT\" but returns: ";
        std::cerr << tmpGene.getSequence() <<"\n";
        error = 1;
    }
    if (!error)
    {
        std::cout <<"Gene reverseComplement --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }

    //-----------------------------------//
    //------ toAASequence Function ------//
    //-----------------------------------//

    if ("MLILTAASX" != testGene.toAASequence())
    {
        std::cerr <<"Error with toAASequence. Should return \"MLILTAASX\", but returns:" << testGene.toAASequence() <<"\n";
        error = 1;
    }

    if (!error)
    {
        std::cout <<"Gene toAASequence --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }



    //----------------------------//
    //------ clear Function ------//
    //----------------------------//

    testGene.clear();
    if ("" != testGene.getId())
    {
        std::cerr <<"Error with clear. Gene Id should be blank, but is " << testGene.getId() <<".\n";
        error = 1;
    }
    if ("" != testGene.getDescription())
    {
        std::cerr <<"Error with clear. Gene description should be blank, but is " << testGene.getDescription() <<".\n";
        error = 1;
    }
    if ("" != testGene.getSequence())
    {
        std::cerr <<"Error with clear. Gene sequence should be blank, but is " << testGene.getSequence() <<".\n";
        error = 1;
    }

    if (!error)
    {
        std::cout <<"Gene clear --- Pass\n";
    }
    else
    {
        error = 0; //Reset for next function.
    }



    //------ cleanSequence Function ------//
  //  testGene.setSequence("AAATTTNGCYRNKROTTN");
   // std::cout <<testGene.getSequence() <<"\n";
}


void testGenome(std::string testFileDir)
{

    int error = 0;
    //--------------------------------//
    //------ readFasta Function ------//
    //--------------------------------//


    Genome genome;
    std::string file = testFileDir + "/" + "test.fasta";
    genome.readFasta(file, false);

    Gene g1("TEST001", "TEST001 Test Gene", "ATGGCCACTATTGGGTCTTAG");
    Gene g2("TEST002", "TEST002 Test Gene", "ATGACCGTAATTTTTTACTAG");
    Gene g3("TEST003", "TEST003 Test Gene", "ATGGTCTACTTTCTGACATAG");

    Genome testGenome;
    testGenome.addGene(g1, false);
    testGenome.addGene(g2, false);
    testGenome.addGene(g3, false);

    if(genome == testGenome)
    {
        std::cout <<"Genome readFasta --- Pass\n";
    }
    else
    {
        std::cerr <<"Error in readFasta. Genomes are not equivelant.\n";
    }
}