/* libraries/libvacuumms_cpp/histogram.cc */

#include <vacuumms/limits.h>
#include <vacuumms/types.h>

#include <vacuumms/histogram.hh>

#include <cmath>
#include <cstdio>


Histogram::Histogram() : bins(100, 0.1)
{
}

Histogram::Histogram(int n_bins, vacuumms_float width) : bins(100, 0.1), number_of_bins(n_bins), width_of_bins(width)
{
}

void Histogram::bin(vacuumms_float value)
{
    values.push_back(value);
}


void Histogram::generate()
{
    // null out existing values
    for (int bin = 0; bin < number_of_bins; bin++) bins[bin] = 0;
    misses = 0;

    for (const auto& value : values) 
    {
        int bin = static_cast<int>(std::floor((value - starting_value) / width_of_bins));
        if ((bin >= number_of_bins) || (bin < 0)) misses++;
        else bins[bin]++;
    }
}


int Histogram::getMisses()
{
    return misses;
}


void Histogram::applyWeightExponent(int exponent)
{
    for (int bin = 0; bin < number_of_bins; bin++)
    {
        vacuumms_float weight = pow((starting_value + (bin * width_of_bins)), exponent);
        bins[bin] *= weight;
    }
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

    for (int bin = 0; bin < number_of_bins; bin++)
        total += bins[bin];
    for (int bin = 0; bin < number_of_bins; bin++)
        bins[bin] /= total; 
}

        
void Histogram::writeToFile(char* filename)
{
    FILE *f = fopen(filename, "w");
    for (int i = 0; i < number_of_bins; i++)
    {
        vacuumms_float value = bins[i];
        fprintf(f, "%f\t%f\n", (vacuumms_float)(i * width_of_bins), value);
    }
    fclose(f);
}


std::vector<std::tuple<vacuumms_float, vacuumms_float>> Histogram::getTuples()
{
    std::vector<std::tuple<vacuumms_float, vacuumms_float>> tuples;

    for (int i = 0; i < number_of_bins; i++)
    {
        vacuumms_float x = starting_value + i * width_of_bins;
        vacuumms_float y = bins[i];
        tuples.emplace_back(x, y);
    }
    return tuples;
}


void Histogram::print()
{
    for (int i = 0; i < number_of_bins; i++)
    {
        printf("%f\t%f\n", (vacuumms_float)(starting_value + (i * width_of_bins)), bins[i]);
    }
    fflush(stdout);
}


void Histogram::setNumberOfBins(int n_bins)
{
    number_of_bins = n_bins;
}


void Histogram::setBinWidth(vacuumms_float width)
{
    width_of_bins = width;
}


void Histogram::setStartingValue(vacuumms_float _starting_value)
{
    starting_value = _starting_value;
}


void Histogram::setValueRange(vacuumms_float _starting_value, vacuumms_float _end_value)
{
    starting_value = _starting_value;
    width_of_bins = (_end_value - _starting_value) / number_of_bins;
}


#ifdef BUILD_PYBIND_BINDINGS 

pybind11::str Histogram::__repr__()
{
    pybind11::str retval;

    // Scaler value used to make distribution output print nicely
    vacuumms_float max_bin_size = 0.0f;
    for (int i = 0; i < number_of_bins; i++) 
        if (bins[i] > max_bin_size) max_bin_size = bins[i];
    scaler = 100 / (max_bin_size);

    for (int i = 0; i < number_of_bins; i++)
    {
        retval += pybind11::str(std::to_string(starting_value + (i * width_of_bins)));
        retval += pybind11::str(":\t");

        for (int j = 1; j < (scaler * bins[i]); j++)
        {
            retval += pybind11::str("*");
        }
        retval += pybind11::str("\n");
    }
    return retval;
}

#endif

