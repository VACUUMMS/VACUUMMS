/* vacuumms/types.cc */

/* Implementation of types from types.hh, representing 
 * types that have not already been otherwise implemented 
 * in C++ 
 */

#include <vacuumms/limits.h>
#include <vacuumms/types.hh>

#include <cmath>
#include <cstdio>

Histogram::Histogram() : bins(100, 0)
{
}

Histogram::Histogram(int n_bins, vacuumms_float width) : bins(100, 0), number_of_bins(n_bins), width_of_bins(width)
{
//    number_of_bins = n_bins;
//    width_of_bins = width;
}

void Histogram::bin(vacuumms_float value)
{
    int bin = static_cast<int>(std::floor(value / width_of_bins));
    if (bin >= number_of_bins) misses++;
    else bins[bin]++;
}

int Histogram::getMisses()
{
    return misses;
}

        
void Histogram::writeToFile(char* filename)
{
    FILE *f = fopen("w", filename);
    for (int i = 0; i < number_of_bins; i++)
        fprintf(f, "%f\t%f\n", (vacuumms_float)(i * width_of_bins), bins[i]);
    fclose(f);
}

void Histogram::setNumberOfBins(int n_bins)
{
    number_of_bins = n_bins;
}

void Histogram::setWidthOfBins(vacuumms_float width)
{
    width_of_bins = width;
}

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

