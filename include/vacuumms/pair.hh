/* vacuumms/pair.hh */

#pragma once

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


class PairCoefficient
{

    public:

        int index;
        vacuumms_float sigma;
        vacuumms_float epsilon;

        PairCoefficient();

        PairCoefficient(int _index, vacuumms_float _sigma, vacuumms_float _epsilon);

};


