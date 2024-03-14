/* vacuumms/vertex.cc */

#include <vacuumms/vertex.hh>

#include <vacuumms/types.h>
#include <vacuumms/limits.h>

#include <vector>
#include <iostream>


Vertex::Vertex(int _index, vacuumms_float _x, vacuumms_float _y, vacuumms_float _z)
{
    index = _index;
    x = _x;
    y = _y;
    z = _z;
}

Vertex::Vertex(int _index, vacuumms_float _x, vacuumms_float _y, vacuumms_float _z, int _foreign_key)
{
    index = _index;
    x = _x;
    y = _y;
    z = _z;
    foreign_key = _foreign_key;
}

VertexList::VertexList()
{
    records = std::vector<Vertex>();
}

VertexList::VertexList(char *filename)
{
    FILE* instream=fopen(filename, "r");
    records = std::vector<Vertex>();
    vacuumms_float x, y, z; //, d, drift;
    int index=0;

    while (!feof(instream))
    {
        //fscanf(instream, "%d\t%f\t%f\t%f\t%f\t%f\n", &index, &x, &y, &z, &d, &drift);
        fscanf(instream, "%f\t%f\t%f\n", &x, &y, &z);
        records.push_back(Vertex(index++, x, y, z));
    }
}
    
VertexList::VertexList(FILE *instream)
{
    records = std::vector<Vertex>();
    vacuumms_float x, y, z;
    int index;

    while (!feof(instream))
    {
        fscanf(instream, "%f\t%f\t%f\n", &x, &y, &z);
        records.push_back(Vertex(index++, x, y, z));
    }
}
    
void VertexList::setBoxDimensions(vacuumms_float _box_x, vacuumms_float _box_y, vacuumms_float _box_z)
{
    box_x = _box_x;
    box_y = _box_y;
    box_z = _box_z;
}

Vertex VertexList::recordAt(int i)
{
    return records[i];
}

void VertexList::deleteRecordAt(int i)
{
    records.erase(records.begin() + i);
}

int VertexList::getSize()
{
    return records.size();
}

int VertexList::pushBack(Vertex _vertex)
{
    records.push_back(_vertex);
    return records.size();
}

