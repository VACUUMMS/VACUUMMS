/* vacuumms/types.hh */

/* Declaring some C types from types.h that have not already been implemented in C++ 
 * Implementation in types.cc
 */

#pragma once

#include <vector>

#include <vacuumms/limits.h>
#include <vacuumms/types.h>

#include <vacuumms/exports.hh>

class 
#ifdef PYBIND11_EXPORTS
PYBIND11_EXPORT
#endif
Histogram
{
    public: 
    
        Histogram();
        Histogram(int n_bins, vacuumms_float width);
        void bin(vacuumms_float value);
        int getMisses();
        void smooth(int);
        void normalize();
        void writeToFile(char* filename);
        void setWeightingExponent(vacuumms_float weight);
        void print();

#ifdef BUILD_PYBIND_BINDINGS 

        pybind11::str __repr__();

#endif
    
    protected:

        void setNumberOfBins(int n_bins);
        void setWidthOfBins(vacuumms_float width);

    private:

        int number_of_bins = 100;
        vacuumms_float width_of_bins = 1.0;
        std::vector<vacuumms_float> bins;
        int misses = 0;
        int scaler = 1;
        vacuumms_float weight = 1.0f;

};

