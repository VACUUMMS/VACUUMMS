/* tcv.cc */

//#define MAX_CAVITIES 1310720
//#define N_SUCCESSES 10000

#include <vector>

#include <vacuumms/types.h>
#include <vacuumms/limits.h>
//#include <vacuumms/std.h>
//#include <vacuumms/rng2.h>

#include <vacuumms/cavity.hh>
#include <vacuumms/param.hh>

// Use functions from the vacuumms C library
extern "C"
{
#include <vacuumms/std.h>
#include <vacuumms/rng2.h>
}

// In:  1 or more records in .cav format, box size on command line
// Out: .dst (reports one single value for volume of all cavities) 


int main(int argc, char* argv[])
{
    // Should be updated to vacuumms_float after vacuumms/types.h is more widely implemented.
    double box_x, box_y, box_z;

    vacuumms_float volume;
    int successes=0, attempts=1000;

    setCommandLineParameters(argc, argv);

    getIntParam((char *)"-attempts", &attempts);
    if (attempts > VACUUMMS_INT_MAX) 
    {
        printf("%d attempts > VACUUMMS_INT_MAX = %d. Please reconsider. exiting.\n", attempts, VACUUMMS_INT_MAX);
        exit(1);
    }
    getVectorParam((char *)"-box", &box_x, &box_y, &box_z);

    if (getFlagParam((char *)"-usage"))
    {
        printf("usage:  tcv [ -attempts nnnnn] [ -box n.nnnnn n.nnnnn n.nnnnn ]\n\n");
        printf("        will return total volume of box which is enclosed by one or more cavities.\n");
        printf("\n");
        printf("// In:  one cluster; 1 or more records in .cav format \n");
        printf("// Out: .dst (reports one value of volume) \n");
        exit(1);
    }
   
    initializeRandomNumberGenerator2(0);

    CavityConfiguration configuration = CavityConfiguration(stdin);
    configuration.setBoxDimensions(box_x, box_y, box_z);

    successes=0;

    int attempts_remaining = attempts;
    while (attempts_remaining-- > 0)
    {
        /* Pick a point and check for inclusion:
         * Add offset of 0.5 to rnd so we sample from the box whose corners 
         * are defined by the centers of the eight mirror boxes. 
         */

        vacuumms_float tx = (0.5 + rnd2()) * box_x;
        vacuumms_float ty = (0.5 + rnd2()) * box_y;
        vacuumms_float tz = (0.5 + rnd2()) * box_z;

        successes += configuration.checkInclusion(tx, ty, tz);
    }

    volume = (box_x * box_y * box_z * successes) / attempts;

    printf("%f\n", volume);
    return 0;
}

