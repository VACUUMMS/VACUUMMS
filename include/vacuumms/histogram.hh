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
        void bin(vacuumms_float value); // add a new value
        void setNumberOfBins(int n_bins);
        void setBinWidth(vacuumms_float width);
        void setStartingValue(vacuumms_float);
        void setValueRange(vacuumms_float, vacuumms_float);

        void generate();

        int getMisses();
        void applyWeightExponent(int);
        void smooth(int);
        void normalize();
        void writeToFile(char* filename);
        std::vector<std::tuple<vacuumms_float, vacuumms_float>> getTuples();
        void print();

#ifdef BUILD_PYBIND_BINDINGS 

        pybind11::str __repr__();

#endif
    
    protected:

        std::vector<vacuumms_float> bins;

    private:

        int number_of_bins = 100;
        vacuumms_float width_of_bins = 1.0;
        vacuumms_float starting_value = 0.0;
        std::vector<vacuumms_float> values;
        int misses = 0;
        vacuumms_float scaler = 1.0f;

};

