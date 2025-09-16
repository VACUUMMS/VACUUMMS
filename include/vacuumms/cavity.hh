// vacuumms/cavity.hh 

#pragma once

#include <vector>
#include <iostream>

#include <vacuumms/types.h>
#include <vacuumms/limits.h>

#include <vacuumms/parameters.hh>
#include <vacuumms/histogram.hh>
#include <vacuumms/exports.hh>


class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
Cavity
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


class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
CavityConfiguration
{
    public:

        CavityConfiguration();
        CavityConfiguration(const char *filename);
        CavityConfiguration(FILE *instream);
        void replicate(std::vector<int> depths);
        void setBoxDimensions(std::vector<vacuumms_float> dims);
        std::vector<vacuumms_float> getBoxDimensions();
        void setMirrorDepth(int _mirror_depth);
        void scrubDuplicates();
        std::vector<vacuumms_float> getDiameters();
        Cavity recordAt(int i);
        void deleteRecordAt(int i);
        int getSize();
        int checkInclusion(vacuumms_float tx, vacuumms_float ty, vacuumms_float tz);
        int pushBack(Cavity _cavity);
        void reset();

        std::vector<vacuumms_float> box_dimensions = {0.0, 0.0, 0.0};
        std::vector<Cavity> records;

#ifdef BUILD_PYBIND_BINDINGS
        pybind11::str __repr__();
#endif

    private:

        int mirror_depth = 1;
        vacuumms_float duplicate_threshold = 0.1f;

}; // end class CavityConfiguration


class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
CavitySizeDistribution : public Histogram
{
    public:

        CavitySizeDistribution(CavityConfiguration cc, Parameters p);
        CavitySizeDistribution(CavityConfiguration cc);

        std::vector<std::tuple<vacuumms_float, vacuumms_float>> getResult();
        
#ifdef BUILD_PYBIND_BINDINGS
        pybind11::str __repr__();
#endif

    private:

        Parameters p;
        CavityConfiguration cc;

}; // end class CavitySizeDistribution


