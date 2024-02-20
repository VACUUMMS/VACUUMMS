/* vacuumms/cavity.hh */

//#define MAX_CAVITIES 1310720
//#define N_SUCCESSES 10000

#include <vacuumms/types.h>
#include <vacuumms/limits.h>

#include <vector>
#include <iostream>

/*
// Use functions from the vacuumms C library
extern "C"
{
#include <ftw_std.h>
#include <ftw_rng2.h>
#include <ftw_param.h>
}
*/


class Cavity
{
    public:

        vacuumms_float x;
        vacuumms_float y;
        vacuumms_float z;
        vacuumms_float d;

    Cavity(vacuumms_float _x, vacuumms_float _y, vacuumms_float _z, vacuumms_float _d);

}; // end class Cavity


class CavityConfiguration
{
    public:

        std::vector<Cavity> records;

        vacuumms_float box_x;
        vacuumms_float box_y;
        vacuumms_float box_z;

        int mirror_depth = 1;

        CavityConfiguration();
        CavityConfiguration(FILE *instream);
        void setBoxDimensions(vacuumms_float _box_x, vacuumms_float _box_y, vacuumms_float _box_z);
        void setMirrorDepth(int _mirror_depth);
        Cavity recordAt(int i);
        void deleteRecordAt(int i);
        int getSize();
        int checkInclusion(vacuumms_float tx, vacuumms_float ty, vacuumms_float tz);

}; // end class CavityConfiguration

