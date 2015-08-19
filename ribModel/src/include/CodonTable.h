#ifndef CodonTable_H
#define CodonTable_H


class CodonTable
{
    private:

        unsigned tableId;

    public:

        //Constructors & destructors:
        explicit CodonTable();
        CodonTable(unsigned _tableId);
        virtual ~CodonTable();
        CodonTable(const CodonTable& other);
        CodonTable& operator=(const CodonTable& other);
};

#endif
