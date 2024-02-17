/* tcv.cc */

#define MAX_CAVITIES 1310720
#define N_SUCCESSES 10000

#include <vacuumms/types.h>
#include <vacuumms/limits.h>

#include <vector>

// Use functions from the vacuumms C library
extern "C"
{
#include <ftw_std.h>
#include <ftw_rng2.h>
#include <ftw_param.h>
}

// In:  1 or more records in .cav format, box size on command line
// Out: .dst (reports one single value for volume of all cavities) 


class Cavity
{
    public:

        vacuumms_float x;
        vacuumms_float y;
        vacuumms_float z;
        vacuumms_float d;

    Cavity(vacuumms_float _x, vacuumms_float _y, vacuumms_float _z, vacuumms_float _d)
    {
        x = _x;
        y = _y;
        z = _z;
        d = _d;
    }
}; // end class Cavity


class CavityConfiguration
{
    public:

        std::vector<Cavity> records;

        vacuumms_float box_x;
        vacuumms_float box_y;
        vacuumms_float box_z;

        int mirror_depth = 1;

        CavityConfiguration()
        {
            records = std::vector<Cavity>();
        }
    
        CavityConfiguration(FILE *instream)
        {
            records = std::vector<Cavity>();
            vacuumms_float x, y, z, d;

            while (!feof(instream))
            {
                fscanf(instream, "%f\t%f\t%f\t%f\n", &x, &y, &z, &d);
                records.push_back(Cavity(x, y, z, d));
            }
        }
    
        void setBoxDimensions(vacuumms_float _box_x, vacuumms_float _box_y, vacuumms_float _box_z)
        {
            box_x = _box_x;
            box_y = _box_y;
            box_z = _box_z;
        }

        void setMirrorDepth(int _mirror_depth)
        {
            mirror_depth = _mirror_depth;
        }

        Cavity recordAt(int i)
        {
            return records[i];
        }

        void deleteRecordAt(int i)
        {
            records.erase(records.begin() + i);
        }

        int getSize()
        {
            return records.size();
        }

        int checkInclusion(vacuumms_float tx, vacuumms_float ty, vacuumms_float tz)
        {
            int i;
            vacuumms_float dx, dy, dz, dd;

            for (i=0; i<getSize(); i++)
            {
                // Check inclusion in each of the eight mirror box images:

                // (0,0,0):
                dx = records[i].x - tx;
                dy = records[i].y - ty;
                dz = records[i].z - tz;
                dd = dx*dx + dy*dy + dz*dz;
                if (4*dd < (records[i].d * records[i].d)) return 1;

                // (0,0,1):
                dx = records[i].x - tx;
                dy = records[i].y - ty;
                dz = box_z + records[i].z - tz;
                dd = dx*dx + dy*dy + dz*dz;
                if (4*dd < (records[i].d * records[i].d)) return 1;

                // (0,1,0):
                dx = records[i].x - tx;
                dy = box_y + records[i].y - ty;
                dz = records[i].z - tz;
                dd = dx*dx + dy*dy + dz*dz;
                if (4*dd < (records[i].d * records[i].d)) return 1;

                // (0,1,1):
                dx = records[i].x - tx;
                dy = box_y + records[i].y - ty;
                dz = box_z + records[i].z - tz;
                dd = dx*dx + dy*dy + dz*dz;
                if (4*dd < (records[i].d * records[i].d)) return 1;

                // (1,0,0):
                dx = box_x + records[i].x - tx;
                dy = records[i].y - ty;
                dz = records[i].z - tz;
                dd = dx*dx + dy*dy + dz*dz;
                if (4*dd < (records[i].d * records[i].d)) return 1;

                // (1,0,1):
                dx = box_x + records[i].x - tx;
                dy = records[i].y - ty;
                dz = box_z + records[i].z - tz;
                dd = dx*dx + dy*dy + dz*dz;
                if (4*dd < (records[i].d * records[i].d)) return 1;

                // (1,1,0):
                dx = box_x + records[i].x - tx;
                dy = box_y + records[i].y - ty;
                dz = records[i].z - tz;
                dd = dx*dx + dy*dy + dz*dz;
                if (4*dd < (records[i].d * records[i].d)) return 1;

                // (1,1,1):
                dx = box_x + records[i].x - tx;
                dy = box_y + records[i].y - ty;
                dz = box_z + records[i].z - tz;
                dd = dx*dx + dy*dy + dz*dz;
                if (4*dd < (records[i].d * records[i].d)) return 1;
            }
 
            return 0;
        }
}; // end class CavityConfiguration


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

