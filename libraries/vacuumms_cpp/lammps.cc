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
                ConfigurationRecord record(x, y, z, pairs[type].sigma, pairs[type].epsilon);
                record.type = type;
                records.push_back(record);
            }
        }

        // now find the pairs

        retval = strstr(line, "Pair Coeffs");
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

    box_x = xhi - xlo;
    box_y = yhi - ylo;
    box_z = zhi - zlo;

/* This isn't necessary because Pair section is read first, and records 
 * are initialized with pair cofficients when Atom section is read.
 
    // use the pair data to finalize sigma and epsilon values
    for (int i=0; i<records.size(); i++)
    {
        records[i].sigma = pairs[records[i].type].sigma;
        records[i].epsilon = pairs[records[i].type].epsilon;
    }
*/

}


#ifdef BUILD_PYBIND_BINDINGS
pybind11::str LAMMPSConfiguration::__repr__()
{
// just dump superclass output
    pybind11::str retval = Configuration::__repr__();
    retval += pybind11::str("box_x: ");
    retval += pybind11::str(std::to_string(box_x));
    retval += pybind11::str("\n");
    retval += pybind11::str("box_y: ");
    retval += pybind11::str(std::to_string(box_y));
    retval += pybind11::str("\n");
    retval += pybind11::str("box_z: ");
    retval += pybind11::str(std::to_string(box_z));
    retval += pybind11::str("\n");
    return retval;
}
#endif

