/* vacuumms/edge.cc */

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
                A_found = i + 1; // Add one to index to return TRUE
            }
            if ( 
                  (xB == vertices.recordAt(i).x) 
               && (yB == vertices.recordAt(i).y) 
               && (zB == vertices.recordAt(i).z) 
               )
            {
                B_found = i + 1; // Add one to index to return TRUE
            }

            if (A_found && B_found) 
            {
                records.push_back(Edge(index++, A_found - 1, B_found - 1));
                break;
            }
        }

        if (!A_found || !B_found) 
        {
            fprintf(stderr, "vacuumms/Edgelist()::EdgeList() record not found. Exiting. \n");
            exit(1);
        }
        else
        {
            //printf("%d resolved.\n", index);
        }

    }
}
    
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

