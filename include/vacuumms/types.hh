/* vacuumms/types.hh */

/* Declaring some C types from types.h that have not already been implemented in C++ 
 * Implementation in types.cc
 */

#pragma once

#include <vacuumms/limits.h>
#include <vacuumms/types.h>

#include <vector>

class Histogram
{
    public: 
    
        Histogram();
        Histogram(int n_bins, vacuumms_float width);
        void bin(vacuumms_float value);
        int getMisses();
        void writeToFile(char* filename);
    
    protected:

        void setNumberOfBins(int n_bins);
        void setWidthOfBins(vacuumms_float width);

        int number_of_bins = 100;
        vacuumms_float width_of_bins = 1.0;
        std::vector<int> bins;
        int misses = 0;
};


/*
  
class EnergyArray
{
//  float energy[][][];
};


class FVI
{
//  float intensity[256][256][256];
};

*/


