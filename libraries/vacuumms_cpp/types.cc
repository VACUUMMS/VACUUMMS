/* vacuumms/types.cc */

/* Implementation of types from types.hh, representing 
 * types that have not already been otherwise implemented 
 * in C++ 
 */

#include <vacuumms/limits.h>
#include <vacuumms/types.h>

#include <vacuumms/types.hh>

#include <cmath>
#include <cstdio>


Histogram::Histogram() : bins(100, 0)
{
}

Histogram::Histogram(int n_bins, vacuumms_float width) : bins(100, 0), number_of_bins(n_bins), width_of_bins(width)
{
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

void Histogram::smooth(int iterations)
{
    std::vector<vacuumms_float> smoothed(bins);

    for (int iteration = 0; iteration < iterations; iteration++)
    {
        smoothed[0] = 0.75 * bins[0] + 0.25 * bins[1];
        for (int bin = 1; bin < number_of_bins - 1; bin++)
            smoothed[bin] = 0.25 * bins[bin - 1] + 0.5 * bins[bin] + 0.25 * bins[bin + 1];
        smoothed[number_of_bins - 1] = 0.75 * bins[number_of_bins - 1] + 0.25 * bins[number_of_bins - 2];
        std::vector<vacuumms_float> copy_of_smoothed(smoothed);
        bins = copy_of_smoothed;
    }
}

void Histogram::normalize()
{
    vacuumms_float total = 0.0f;

    for (int i = 0; i < number_of_bins; i++)
        total += pow(bins[i], weight);
    for (int i = 0; i < number_of_bins; i++)
        bins[i] /= total; 
    
    scaler = number_of_bins;
}

        
void Histogram::writeToFile(char* filename)
{
    FILE *f = fopen(filename, "w");
    for (int i = 0; i < number_of_bins; i++)
    {
        vacuumms_float value = pow(bins[i], weight);
        fprintf(f, "%f\t%f\n", (vacuumms_float)(i * width_of_bins), value);
    }
    fclose(f);
}

void Histogram::print()
{
    for (int i = 0; i < number_of_bins; i++)
    {
        vacuumms_float value = pow(bins[i], weight);
        printf("%f\t%f\n", (vacuumms_float)(i * width_of_bins), value);
    }
}

void Histogram::setNumberOfBins(int n_bins)
{
    number_of_bins = n_bins;
}

void Histogram::setWidthOfBins(vacuumms_float width)
{
    width_of_bins = width;
}


void Histogram::setWeightingExponent(vacuumms_float _weight)
{
    weight = _weight;
}

#ifdef BUILD_PYBIND_BINDINGS 

pybind11::str Histogram::__repr__()
{
    pybind11::str retval;

    for (int i = 0; i < number_of_bins; i++)
    {
        retval += pybind11::str(std::to_string(i * width_of_bins));
        retval += pybind11::str(":\t");
        for (int j = 0; j < scaler * pow(bins[i], weight); j++)
        {
            retval += pybind11::str("*");
        }
        retval += pybind11::str("\n");
    }
    return retval;
}

#endif

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

