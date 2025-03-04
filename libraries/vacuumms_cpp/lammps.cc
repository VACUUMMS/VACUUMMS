/* libraries/vacuumms_cpp/lammps.cc */

#include <vacuumms/param.hh>
#include <vacuumms/types.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <map>
#include <vector>
#include <cstring>

#include <vacuumms/lammps.hh>
#include <vacuumms/pair.hh>


LAMMPSConfiguration::LAMMPSConfiguration(std::string filename)
{

    FILE* lmps_file = fopen(filename.c_str(), "r");

    char line[256];
    int index;
    vacuumms_float epsilon, sigma;


    // Parse the dang file, start to finish 
    while(fgets(line, 256, lmps_file))
    {
        char* retval;

        // Grab the box size params

        retval = strstr(line, "xlo");
        if (retval !=NULL)
        {
            sscanf(line, "%f %f xlo xhi", &xlo, &xhi);
        }

        retval = strstr(line, "ylo");
        if (retval !=NULL)
        {
            sscanf(line, "%f %f ylo yhi", &ylo, &yhi);
        }

        retval = strstr(line, "zlo");
        if (retval !=NULL)
        {
            sscanf(line, "%f %f zlo zhi", &zlo, &zhi);
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
                ConfigurationRecord record(x, y, z, 0, 0);
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

    // use the pair data to finalize sigma and epsilon values
    for (int i=0; i<records.size(); i++)
    {
        records[i].sigma = pairs[records[i].type].sigma;
        records[i].epsilon = pairs[records[i].type].epsilon;

        // write the box data 
//        printf("%f\t%f\t%f\t%f\t%f\n", atoms[i].x, atoms[i].y, atoms[i].z, atoms[i].sigma, atoms[i].epsilon);
    }
}

