#include <vector>
#include <stdio.h>
#include <vacuumms/types.h>

class ConfigurationRecord
{
    public:

        vacuumms_float x;
        vacuumms_float y;
        vacuumms_float z;
        vacuumms_float sigma;
        vacuumms_float epsilon;
//        vacuumms_float sigma_6;
//        vacuumms_float sigma_12;

        ConfigurationRecord(vacuumms_float _x, vacuumms_float _y, vacuumms_float _z, vacuumms_float _sigma, vacuumms_float _epsilon);
};

class Configuration
{
        std::vector<ConfigurationRecord> records;

        vacuumms_float box_x;
        vacuumms_float box_y;
        vacuumms_float box_z;

        int mirror_depth = 1;
    
    public:

        Configuration(char *filename);
        Configuration(FILE *pipe); // allows stdin to be used to create pipeline
        Configuration();
        void dumpContents();
        vacuumms_float insertionEnergy(vacuumms_float x, vacuumms_float y, vacuumms_float z, vacuumms_float sigma, vacuumms_float epsilon);
        void setBoxDimensions(vacuumms_float _box_x, vacuumms_float _box_y, vacuumms_float _box_z);
        void setMirrorDepth(int _mirror_depth);

        ConfigurationRecord recordAt(int i);
        void deleteRecordAt(int i);
        int getSize();
};


