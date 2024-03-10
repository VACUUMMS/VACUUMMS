/* edges2pairs.cc */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Hacking to interface to C linkage;
extern "C"
{
#include <ftw_param.h>
}

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

    // we've got vertices and edge lists. Now we refine:

    // load the vcav output, so we can map verts to cavities. :x

//for (int i=0; i<edges.getSize(); i++) printf("(%d, %d)\n", edges.recordAt(i).A, edges.recordAt(i).B);
    



  return 0;
}
