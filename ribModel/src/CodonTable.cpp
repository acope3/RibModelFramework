#include "include/CodonTable.h"


CodonTable::CodonTable()
{
    tableId = 1; //standard codon table by NCBI
}


CodonTable::CodonTable(unsigned _tableId) : tableId(_tableId)
{
    //ctor
}


CodonTable::~CodonTable()
{
    //dtor
}


CodonTable::CodonTable(const CodonTable& other)
{
    tableId = other.tableId;
}


CodonTable& CodonTable::operator=(const CodonTable& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    tableId = rhs.tableId;
}