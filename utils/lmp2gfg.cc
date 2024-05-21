// lmp2gfg.cc

#include <vacuumms/param.hh>
#include <vacuumms/types.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <map>
#include <vector>
#include <cstring>

class PairCoefficient
{
    
    public:

        int index;
        vacuumms_float sigma;
        vacuumms_float epsilon;

        PairCoefficient()
        {
        }

        PairCoefficient(int _index, vacuumms_float _sigma, vacuumms_float _epsilon)
        {
            index = _index;
            sigma = _sigma;
            epsilon = _epsilon;
        }

};


int main(int argc, char *argv[]) 
{
//    double box_x=10, box_y=10, box_z=10;
//    char* pair_filename = (char*)"pair.dat";
//    char* lmps_filename = (char*)"in.lmps";

    setCommandLineParameters(argc, argv);

    if (getFlagParam((char*)"-usage"))
    {
        printf("\n");
        printf("usage:     lmp2gfg        \n");
//        printf("                          +cram                 \n");
        printf("\n");
        exit(0);
    }

    int cram = getFlagParam((char*)"+cram");
//    getVectorParam((char*)"-box", &box_x, &box_y, &box_z);
//    getStringParam((char*)"-pair_file", &pair_filename);
//    getStringParam((char*)"-lmps_file", &lmps_filename);

//    FILE* pair_file = fopen(pair_filename, "r");
//    assert(pair_file);

//    FILE* lmps_file = fopen(lmps_filename, "r");
//    assert(lmps_file);

    FILE* lmps_file = stdin;

    char line[256];
    int index;
    vacuumms_float epsilon, sigma;

    std::map<int, PairCoefficient> pairs;
    std::vector<Atom> atoms;

    vacuumms_float xlo, xhi;
    vacuumms_float ylo, yhi;
    vacuumms_float zlo, zhi;

    // Parse the dang file, start to finish 
    while(fgets(line, 256, lmps_file))
    {
        char* retval;

        // Grab the box size params

        retval = strstr(line, "xlo");
        if (retval !=NULL)
        {
            sscanf(line, "%f %f xlo xhi", &xlo, &xhi);
            fprintf(stderr, "Got xlo, xhi = %f, %f\n", xlo, xhi);
        }

        retval = strstr(line, "ylo");
        if (retval !=NULL)
        {
            sscanf(line, "%f %f ylo yhi", &ylo, &yhi);
            fprintf(stderr, "Got ylo, yhi = %f, %f\n", ylo, yhi);
        }

        retval = strstr(line, "zlo");
        if (retval !=NULL)
        {
            sscanf(line, "%f %f zlo zhi", &zlo, &zhi);
            fprintf(stderr, "Got zlo, zhi = %f, %f\n", zlo, zhi);
        }

        // Grab the list of Atoms

        retval = strstr(line, "Atoms");
        if (retval != NULL) // found it, so do the work
        {
            // skip the blank line .. not needed?
            // fgets(line, 256, lmps_file);

            int dummy1, dummy2, dummy3, dummy4, dummy5, dummy6;
            int type;
            vacuumms_float f_dummy, x, y, z;

            // fscanf returns number of items read, so when it reaches 
            // the end of the list, it will return 0, dropping back to
            // the outer loop, with seek set to end of atom list.
            while(fscanf(lmps_file, "%d %d %d %f %f %f %f %d %d %d", 
                         &dummy1, &dummy2, &type, &f_dummy, 
                         &x, &y, &z, &dummy4, &dummy5, &dummy6) > 0)
            {
                Atom a;
                a.type = type;
                a.x = x;
                a.y = y;
                a.z = z;
                atoms.push_back(a);
            }
        }

        // now find the pairs

        retval = strstr(line, "Pair");
        if (retval != NULL) 
        {
            // found it, so do the work
            int atom_type;
            vacuumms_float sigma, epsilon;

            // Scanf will return 0 when nothing is matched, i.e. at the end of data
            while(fscanf(lmps_file, "%d %f %f", &atom_type, &epsilon, &sigma) > 0)
            {
                pairs[atom_type] = PairCoefficient(atom_type, sigma, epsilon);
            }
        }
    } // loop back to beginning of while

    vacuumms_float box_x = xhi - xlo;
    vacuumms_float box_y = yhi - ylo;
    vacuumms_float box_z = zhi - zlo;

    fprintf(stderr, "-box %f %f %f\n", box_x, box_y, box_z);

    // cram the box?
    if (cram)
    {
    }

    // use the pair data to finalize the atoms
    for (int i=0; i<atoms.size(); i++)
    {
        atoms[i].sigma = pairs[atoms[i].type].sigma;
        atoms[i].epsilon = pairs[atoms[i].type].epsilon;

        // write the box data 
        printf("%f\t%f\t%f\t%f\t%f\n", atoms[i].x, atoms[i].y, atoms[i].z, atoms[i].sigma, atoms[i].epsilon);
    }

    return 0;

}
