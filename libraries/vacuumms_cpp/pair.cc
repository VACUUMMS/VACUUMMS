/* vacuumms/pair.cc */

#include <vacuumms/pair.hh>

#include <vacuumms/types.h>
#include <vacuumms/limits.h>

#include <vector>
#include <iostream>


IndexPair::IndexPair(int _A, int _B)
{
    A = _A;
    B = _B;
}

IndexPairList::IndexPairList(char *filename)
{
    records = std::vector<IndexPair>();
    FILE* instream=fopen(filename, "r");

    while (!feof(instream))
    {
        int A, B;
        fscanf(instream, "%d\t%d\n", &A, &B);
        records.push_back(IndexPair(A, B));
    }
}

IndexPairList::IndexPairList()
{
    records = std::vector<IndexPair>();
}

IndexPair IndexPairList::recordAt(int i)
{
    return records[i];
}

void IndexPairList::deleteRecordAt(int i)
{
    records.erase(records.begin() + i);
}

int IndexPairList::getSize()
{
    return records.size();
}

int IndexPairList::pushBack(IndexPair _index_pair)
{
    records.push_back(_index_pair);
    return records.size();
}

