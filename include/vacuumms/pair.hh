/* vacuumms/pair.hh */

#include <vacuumms/types.h>
#include <vacuumms/limits.h>

#include <vector>
#include <iostream>


class IndexPair
{
    public:

        int A;
        int B;

        IndexPair(int _A, int _B);

}; // end class IndexPair


class IndexPairList
{
    public:

        std::vector<IndexPair> records;

        IndexPairList();
        IndexPairList(char *filename);
        IndexPair recordAt(int i);
        void deleteRecordAt(int i);
        int getSize();
        int pushBack(IndexPair _index_pair);

}; // end class IndexPairList

