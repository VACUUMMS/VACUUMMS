/* vacuumms/edge.hh */

#pragma once

#include <vector>
#include <iostream>

#include <vacuumms/types.h>
#include <vacuumms/limits.h>


class Edge
{
    public:

        int index;
        int A;
        int B;
        int foreign_key;

        Edge(int _index, int _A, int _B);

}; // end class Edge


class EdgeList
{
    public:

        std::vector<Edge> records;
        VertexList vertices;

        EdgeList();
        EdgeList(char *filename, VertexList vertices);
        Edge recordAt(int i);
        void deleteRecordAt(int i);
        int getSize();

}; // end class EdgeList

