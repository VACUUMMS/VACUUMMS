// configuration.hh 

#pragma once

#include <vector>
#include <cstdio>

#include <vacuumms/types.h>

#include <vacuumms/exports.hh>


class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
ConfigurationRecord
{
    public:

        vacuumms_float x;
        vacuumms_float y;
        vacuumms_float z;
        vacuumms_float sigma;
        vacuumms_float epsilon;
        int type;
//        void set_x();
//        void set_y();
//        void set_z();
//        void set_sigma();
//        void set_epsilon();

        ConfigurationRecord();
        ConfigurationRecord(vacuumms_float _x, vacuumms_float _y, vacuumms_float _z, vacuumms_float _sigma, vacuumms_float _epsilon);

        std::vector<vacuumms_float> getXYZ();
};

class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
Configuration
{
    protected:

        int crammed = 0;
        int mirror_depth = 1;

    public:
   
        std::vector<vacuumms_float> box_dimensions = {0.0, 0.0, 0.0};
        std::vector<ConfigurationRecord> records;
        vacuumms_float temperature = 1.0f;

        Configuration(const char *filename);
        Configuration(const char *filename, std::vector<vacuumms_float>);
        Configuration(const char *filename, std::vector<vacuumms_float>, vacuumms_float);
        Configuration(FILE *pipe); // allows stdin to be used to create pipeline
        Configuration();
        void dumpContents();
        void writeToFile(const char*);

        vacuumms_float insertionEnergy(vacuumms_float x, vacuumms_float y, vacuumms_float z, vacuumms_float sigma, vacuumms_float epsilon);
        void setMirrorDepth(int _mirror_depth);
        void setTemperature(vacuumms_float);
        vacuumms_float getTemperature();
        void setBoxDimensions(std::vector<vacuumms_float> dims);
        std::vector<vacuumms_float> getBoxDimensions();
        void shift(std::vector<vacuumms_float>);
        ConfigurationRecord recordAt(int i);
        void deleteRecordAt(int i);
        void clear();
        int getSize();
        int pushBack(ConfigurationRecord);
        void cram();
        int isCrammed();
        // void replicate(int depth);
        void replicate(std::vector<int> depths);

#ifdef BUILD_PYBIND_BINDINGS
        pybind11::str __repr__();
#endif

};

