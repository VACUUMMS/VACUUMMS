/* vuniq.cc */

extern "C"
{
#include <ftw_param.h>
}

#include <vacuumms/types.h>
#include <vacuumms/cavity.hh>
#include <vacuumms/vertex.hh>

#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <unistd.h>

double threshold = .01;
double min_diam = 0.0;

FILE *instream;

double box_x=0.0;
double box_y=0.0;
double box_z=0.0;

int main(int argc, char *argv[])
{
    double shift_x, shift_y, shift_z;

    instream=stdin;

    setCommandLineParameters(argc, argv);
    if (getFlagParam((char*)"-usage")) 
    {
      printf("\n");
      printf("usage:\t-box [ 10.0 10.0 10.0 ]\n");  
      printf("      \t-threshold [ .01 ]\n");  
      printf("      \t-min_diam [ 0.0 ]\n\n");
      exit(0);
    }

    getVectorParam((char*)"-box", &box_x, &box_y, &box_z);
    getDoubleParam((char*)"-threshold", &threshold);

    CavityConfiguration cavities;
    VertexList map;

    while (!feof(instream))
    {
        int vertex_index;
        vacuumms_float x, y, z, d, drift;
        fscanf(stdin, "%d\t%f\t%f\t%f\t%f\t%f\n", &vertex_index, &x, &y, &z, &d, &drift);

        int duplicate_found = 0;

        // seek all cavs up to the current record
        for (int i = 0; i < cavities.getSize(); i++)
        {

            for (shift_x=-box_x; shift_x<=box_x; shift_x += box_x)
            for (shift_y=-box_y; shift_y<=box_y; shift_y += box_y)
            for (shift_z=-box_z; shift_z<=box_z; shift_z += box_z)
            {
                vacuumms_float dsq = (shift_x + x - cavities.recordAt(i).x) * (shift_x + x - cavities.recordAt(i).x)
                                   + (shift_y + y - cavities.recordAt(i).y) * (shift_y + y - cavities.recordAt(i).y)
                                   + (shift_z + z - cavities.recordAt(i).z) * (shift_z + z - cavities.recordAt(i).z);
                if (dsq < threshold)
                {
                    duplicate_found = 1;
                    // we know that this vertex maps to the cavity at i so we add this to the map
//                    printf("###VERTEX2CAVITY\t%d\t%d\n", vertex_index, i);
                    map.pushBack(Vertex(vertex_index, x, y, z, i));
                    goto end_search;
                }
            }
        }

end_search:

        // no duplicate found, so add the new cavity and record the index
        if (!duplicate_found) 
        {        
            int cavity_index = cavities.pushBack(Cavity(x,y,z,d)) - 1;
            map.pushBack(Vertex(vertex_index, x, y, z, cavity_index));
//            printf("###VERTEX2CAVITY\t%d\t%d\n", vertex_index, cavity_index);
        }
    }

    fclose(instream);

    for (int i=0; i<cavities.getSize(); i++) printf("%d\t%lf\t%lf\t%lf\t%lf\n", i, cavities.recordAt(i).x, cavities.recordAt(i).y, cavities.recordAt(i).z, cavities.recordAt(i).d);

    for (int i=0; i<map.getSize(); i++) printf("Vertex-->Cavity map: %d --> %d\n", map.recordAt(i).index, map.recordAt(i).foreign_key);
}

