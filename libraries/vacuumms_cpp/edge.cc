/* vacuumms/vertex.cc */

#include <vacuumms/vertex.hh>
#include <vacuumms/edge.hh>

#include <vacuumms/types.h>
#include <vacuumms/limits.h>

#include <vector>
#include <iostream>


Edge::Edge(int _index, int _A, int _B)
{
    index = _index;
    A = _A;
    B = _B;
}

EdgeList::EdgeList()
{
    records = std::vector<Edge>();
}

EdgeList::EdgeList(char *filename, VertexList _vertices)
{
    vertices = _vertices;

    FILE* instream=fopen(filename, "r");
    records = std::vector<Edge>();
    int index = 0;
    vacuumms_float xA, yA, zA, xB, yB, zB;

    while (!feof(instream))
    {
        fscanf(instream, "%f\t%f\t%f\t%f\t%f\t%f\n", &xA, &yA, &zA, &xB, &yB, &zB); 
        //fscanf(instream, "%d\t%d\n", &A, &B);
        //records.push_back(Edge(index++, x, y, z));
// need to look up the record.. just a placeholder now
//printf("got %d\n", index);

        int A_found = 0;
        int B_found = 0;

        // Search the vertices to get A and B indices
        for (int i = 0; i < vertices.getSize(); i++)
        {
            if ( 
                  (xA == vertices.recordAt(i).x) 
               && (yA == vertices.recordAt(i).y) 
               && (zA == vertices.recordAt(i).z) 
               )
            {
                A_found = i + 1;
            }
            if ( 
                  (xB == vertices.recordAt(i).x) 
               && (yB == vertices.recordAt(i).y) 
               && (zB == vertices.recordAt(i).z) 
               )
            {
                B_found = i + 1;
            }

            if (A_found && B_found) 
            {
                records.push_back(Edge(index++, A_found, B_found));
                break;
            }
        }

        if (!A_found || !B_found) 
        {
            fprintf(stderr, "vacuumms/Edgelist()::EdgeList() record not found. Exiting. \n");
        }
        else
        {
            //printf("%d resolved.\n", index);
        }

    }
}
    
/*
EdgeList::EdgeList(FILE *instream)
{
    records = std::vector<Edge>();
    vacuumms_float x, y, z;
    int index;

    while (!feof(instream))
    {
        fscanf(instream, "%f\t%f\t%f\n", &x, &y, &z);
        records.push_back(Edge(index++, x, y, z));
    }
}
*/
    
Edge EdgeList::recordAt(int i)
{
    return records[i];
}

void EdgeList::deleteRecordAt(int i)
{
    records.erase(records.begin() + i);
}

int EdgeList::getSize()
{
    return records.size();
}

