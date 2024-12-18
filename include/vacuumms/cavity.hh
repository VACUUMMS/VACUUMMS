/* vacuumms/cavity.hh */

#include <vacuumms/types.h>
#include <vacuumms/limits.h>

#include <vector>
#include <iostream>


class Cavity
{
    public:

        vacuumms_float x;
        vacuumms_float y;
        vacuumms_float z;
        vacuumms_float d;
        vacuumms_float drift;

        int index;
        int foreign_key;

        Cavity(vacuumms_float _x, vacuumms_float _y, vacuumms_float _z, vacuumms_float _d);
        Cavity(int _index, 
               vacuumms_float _x, 
               vacuumms_float _y, 
               vacuumms_float _z, 
               vacuumms_float _d, 
               vacuumms_float _drift);

        void setForeignKey(int _foreign_key);
        int getForeignKey();

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
        CavityConfiguration(char *filename);
        CavityConfiguration(FILE *instream);
        void setBoxDimensions(vacuumms_float _box_x, vacuumms_float _box_y, vacuumms_float _box_z);
        void setMirrorDepth(int _mirror_depth);
        Cavity recordAt(int i);
        void deleteRecordAt(int i);
        int getSize();
        int checkInclusion(vacuumms_float tx, vacuumms_float ty, vacuumms_float tz);
        int pushBack(Cavity _cavity);

}; // end class CavityConfiguration

