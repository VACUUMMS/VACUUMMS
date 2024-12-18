/* vacuumms/vertex.hh */

#include <vacuumms/types.h>
#include <vacuumms/limits.h>

#include <vector>
#include <iostream>


class Vertex
{
    public:

        int index;
        vacuumms_float x;
        vacuumms_float y;
        vacuumms_float z;
        int foreign_key;

        Vertex(int _index, vacuumms_float _x, vacuumms_float _y, vacuumms_float _z);
        Vertex(int _index, vacuumms_float _x, vacuumms_float _y, vacuumms_float _z, int _foreign_key);

}; // end class Vertex


class VertexList
{
    public:

        std::vector<Vertex> records;

        vacuumms_float box_x;
        vacuumms_float box_y;
        vacuumms_float box_z;

        VertexList();
        VertexList(char *filename);
        VertexList(FILE *instream);
        void setBoxDimensions(vacuumms_float _box_x, vacuumms_float _box_y, vacuumms_float _box_z);
        Vertex recordAt(int i);
        void deleteRecordAt(int i);
        int getSize();
        int pushBack(Vertex _vertex);

}; // end class VertexList

