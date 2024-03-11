/* edges2pairs.cc */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <vacuumms/param.hh>

#include <vacuumms/cavity.hh>
#include <vacuumms/vertex.hh>
#include <vacuumms/edge.hh>
#include <vacuumms/types.h>

double box_x=0.0, box_y=0.0, box_z=0.0;
char *verts_filename;
char *edges_filename;

int main(int argc, char *argv[])
{

  // Parse command line options

    setCommandLineParameters(argc, argv);
    if (getFlagParam((char*)"-usage"))
    {
        printf("usage:       edges2pairs");
        printf("             -box [ n.nn n.nn n.nn ] \n");
        printf("             -edges_file [ edges.edg ] \n");
        printf("             -verts_file [ vertices.vrt ] \n");
        printf("             < vcav (read formatted input from stdin) \n");
        printf("             [ where .vcav format is:\n");
        printf("               %%d\t%%f\t%%f\t%%f\t%%f\t%%f\n");
        printf("               i.e. index, x, y, z, d, drift ]\n\n"); 
        return 0;
    }

    getVectorParam((char*)"-box", &box_x, &box_y, &box_z);
    if (box_x * box_y * box_z <= 0.0) 
    {
        printf("vacuumms/edges2pairs: please specify box size.\n");
        exit(2);
    }

    getStringParam((char *)"-verts_file", &verts_filename);
    getStringParam((char *)"-edges_file", &edges_filename);

    if (verts_filename == NULL || edges_filename == NULL) 
    {
        printf("vacuumms/edges2pairs: please specify verts_file and edges_file.\n");
        exit(2);
    }

printf("got verts file %s, edges file %s\n", verts_filename, edges_filename);

    VertexList verts(verts_filename);

printf("verts size: %d\n", verts.getSize());

    EdgeList edges(edges_filename, verts);

printf("edges size: %d\n", edges.getSize());

    // we've got vertices and edges lists. Now we refine:

    // load the vcav output, so we can map verts to cavities. 

    CavityConfiguration cavs;

    int index;
    vacuumms_float x, y, z, d, drift;

    while (!feof(stdin))
    {
        fscanf(stdin, "%d%f\t%f\t%f\t%f%f\n", &index, &x, &y, &z, &d, &drift);
        int record_number = cavs.getSize();
        cavs.records.push_back(Cavity(x, y, z, d));
        cavs.recordAt(record_number).setForeignKey(index);
    }

printf("cavs size: %d\n", cavs.getSize());

// everything loaded in, so now start culling pairs for output

/*
// generate map of vertIDs to cavityIDs
//FTW BASED ON CAVITY START POINT!!!
    for (int i = 0; i < verts.getSize(); i++)
    {
        int matches = 0;
        for (int j = 0; j < cavs.getSize(); j++)
        {
            // does it match?
            if 
            (
                (verts.elementAt[i].x == cavs.elementAt[j].x) && 
                (verts.elementAt[i].y == cavs.elementAt[j].y) &&
                (verts.elementAt[i].z == cavs.elementAt[j].z)
            )
            {
                // matches, set the foreign key to cavity ID
                matches = 1;
                verts.elementAt[i].foreign_key = j;
                break; // process next vert
            }
        }
        if (!matches)
        {
            printf("no match found for verts[%d]\n", i);
            verts.elementAt[i].foreign_key = -1;
        } // next cav
    } // next vert

// loop over edges to generate list of vertID pairs

// remove duplicates from vertID pairs

// turn vertID pairs back into vert pairs

    for (int i = 0; i < pairs.getSize(); i++)
    {
        vacuumms_float x1 = cavs.elementAt[pairs[i].
        printf("%f\t%f\t%f\t%f\t%f\t%f\n", cavs.elementAt[
*/

    



  return 0;
}
