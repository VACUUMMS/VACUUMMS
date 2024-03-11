/* vuniq.cc */

#include <vacuumms/param.hh>

#include <vacuumms/types.h>
#include <vacuumms/cavity.hh>
#include <vacuumms/vertex.hh>

#include <stdio.h>

double threshold = .01;
double min_diam = 0.0;

FILE *instream;

double box_x=0.0;
double box_y=0.0;
double box_z=0.0;

int main(int argc, char *argv[])
{
    double shift_x, shift_y, shift_z;

    instream = stdin;

    setCommandLineParameters(argc, argv);
    if (getFlagParam((char*)"-usage")) 
    {
        printf("\n");
        printf("vuniq usage:    -box [ 10.0 10.0 10.0 ]\n");  
        printf("                -threshold [ 0.01 ]\n");  
        printf("                -min_diam [ 0.0 ]\n");
        printf("                -mapfile_name filename.map\n\n");
        printf("                IN:  read from stdin < .vcv as: %%d     %%f     %%f     %%f     %%f     %%f\n");
        printf("                                                idx     x      y      z      d      drift\n");
        printf("                OUT: write to stdout > .vnq as: %%d     %%f     %%f     %%f     %%f     %%f\n");
        printf("                                                idx     x      y      z      d      drift\n");

        exit(0);
    }

    getVectorParam((char*)"-box", &box_x, &box_y, &box_z);
    if (box_x * box_y * box_z <= 0.0)
    {
        printf("vacuumms/vuniq: box size specified is not appropriate: %f x %f x %f\n", box_x, box_y, box_z);
        exit(2);
    }

    getDoubleParam((char*)"-threshold", &threshold);

    char* mapfile_name = NULL;
    getStringParam((char*)"-mapfile_name", &mapfile_name);

    CavityConfiguration unique;
    VertexList map;

    // Read the input stream and generate unique list of cavs and map entries on the fly
    while (!feof(instream))
    {
        int vertex_index;
        vacuumms_float x, y, z, d, drift;
        fscanf(instream, "%d\t%f\t%f\t%f\t%f\t%f\n", &vertex_index, &x, &y, &z, &d, &drift);

        int duplicate_found = 0;

        // seek all cavs up to the current record
        for (int i = 0; i < unique.getSize(); i++)
        {

            for (shift_x=-box_x; shift_x<=box_x; shift_x += box_x)
            for (shift_y=-box_y; shift_y<=box_y; shift_y += box_y)
            for (shift_z=-box_z; shift_z<=box_z; shift_z += box_z)
            {
                vacuumms_float dsq = (shift_x + x - unique.recordAt(i).x) * (shift_x + x - unique.recordAt(i).x)
                                   + (shift_y + y - unique.recordAt(i).y) * (shift_y + y - unique.recordAt(i).y)
                                   + (shift_z + z - unique.recordAt(i).z) * (shift_z + z - unique.recordAt(i).z);
                if (dsq < threshold)
                {
                    duplicate_found = 1;
                    map.pushBack(Vertex(vertex_index, x, y, z, i));
                    goto end_search;
                }
            }
        }

end_search:

        // no duplicate found, so add the new cavity and record the index
        if (!duplicate_found) 
        {        
            unique.pushBack(Cavity(vertex_index,x,y,z,d, drift)) - 1;
            map.pushBack(Vertex(vertex_index, x, y, z, vertex_index));
        }
    }

    fclose(instream);

    for (int i=0; i<unique.getSize(); i++) 
        printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", unique.recordAt(i).index, 
                                                unique.recordAt(i).x, 
                                                unique.recordAt(i).y, 
                                                unique.recordAt(i).z, 
                                                unique.recordAt(i).d, 
                                                unique.recordAt(i).drift);

    if (mapfile_name != NULL)
    {
        FILE* mapfile = fopen(mapfile_name, "w");

        for (int i=0; i<map.getSize(); i++) 
            fprintf(mapfile, "%d\t%d\n", map.recordAt(i).index, map.recordAt(i).foreign_key);

        fflush(mapfile);
        fclose(mapfile);
    }
}


