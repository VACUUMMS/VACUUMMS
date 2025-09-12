#include <iostream>
#include <cmath>
#include <vacuumms/configuration.hh>
#include "vacuumms/types.h"


ConfigurationRecord::ConfigurationRecord(vacuumms_float _x, vacuumms_float _y, vacuumms_float _z, vacuumms_float _sigma, vacuumms_float _epsilon)
{
    x = _x;
    y = _y;
    z = _z;
    sigma = _sigma;
    epsilon = _epsilon;
}


Configuration::Configuration()
{
    // create the records object but don't populate
    records = std::vector<ConfigurationRecord>();    
}


Configuration::Configuration(const char *filename)
{
    FILE* infile = fopen(filename, "r");
    if (infile == NULL)
    {
        printf("Failed to open file: %s\n", strerror(errno)); // Print error message
        fflush(stdout);
    }
    else
    {
        vacuumms_float x, y, z, sigma, epsilon;
        records = std::vector<ConfigurationRecord>();    

        while (!feof(infile))
        {
            fscanf(infile, "%f\t%f\t%f\t%f\t%f\n", &x, &y, &z, &sigma, &epsilon);
            records.push_back(ConfigurationRecord(x, y, z, sigma, epsilon));
        }

        fclose(infile);
    }
}


Configuration::Configuration(const char *filename, std::vector<vacuumms_float> _box_dimensions) 
    : Configuration(filename)
{
    box_dimensions = _box_dimensions;
}    


Configuration::Configuration(const char *filename, std::vector<vacuumms_float> _box_dimensions, vacuumms_float _temperature) 
    : Configuration(filename, _box_dimensions)
{
    temperature = _temperature;
}


Configuration::Configuration(FILE *pipe)
{
    vacuumms_float x, y, z, sigma, epsilon;
    records = std::vector<ConfigurationRecord>();    

    while (!feof(pipe))
    {
        fscanf(pipe, "%f\t%f\t%f\t%f\t%f\n", &x, &y, &z, &sigma, &epsilon);
        records.push_back(ConfigurationRecord(x, y, z, sigma, epsilon));
    }
}


void Configuration::setTemperature(vacuumms_float _temperature)
{
    temperature = _temperature;
}


vacuumms_float Configuration::getTemperature()
{
    return temperature;
}


void Configuration::dumpContents()
{
    for (int i = 0; i < records.size(); i++)
        printf("%f\t%f\t%f\t%f\t%f\n", records[i].x, records[i].y, records[i].z, records[i].sigma, records[i].epsilon);
}


vacuumms_float Configuration::insertionEnergy(vacuumms_float x, vacuumms_float y, vacuumms_float z, vacuumms_float sigma, vacuumms_float epsilon)
{
    vacuumms_float total = 0.0;
    for (int box_i=-mirror_depth; box_i<=mirror_depth; box_i++)
        for (int box_j=-mirror_depth; box_j<=mirror_depth; box_j++)
            for (int box_k=-mirror_depth; box_k<=mirror_depth; box_k++)
                for (int i=0; i<records.size(); i++)
                {
                    vacuumms_float test_x = box_i * box_dimensions[0] + records[i].x;
                    vacuumms_float test_y = box_j * box_dimensions[1] + records[i].y;
                    vacuumms_float test_z = box_k * box_dimensions[2] + records[i].z;
                    vacuumms_float r_sq = (test_x - x) * (test_x - x)
                               + (test_y - y) * (test_y - y)
                               + (test_z - z) * (test_z - z);

                    // Lorentz-Berthelot combining rules for sigma and epsilon
                    vacuumms_float sigma_ij = 0.5 * (sigma + records[i].sigma);
                    vacuumms_float sigma_sq = sigma_ij * sigma_ij;
                    vacuumms_float epsilon_ij = sqrt(epsilon * records[i].epsilon);

                    vacuumms_float sigma_6 = sigma_sq * sigma_sq * sigma_sq;
                    vacuumms_float sigma_12 = sigma_6 * sigma_6;
                    vacuumms_float r_6 = r_sq * r_sq * r_sq;
                    vacuumms_float r_12 = r_6 * r_6;

                    total += 4 * epsilon_ij * (sigma_12/r_12 - sigma_6/r_6);
                }
    return total;
}


void Configuration::setBoxDimensions(std::vector<vacuumms_float> dims)
{
//    box_x = dims[0];
//    box_y = dims[1];
//    box_z = dims[2];

    box_dimensions = dims;
}


std::vector<vacuumms_float> Configuration::getBoxDimensions()
{
/*
    std::vector<float> arr(3);
    arr[0] = box_x;
    arr[1] = box_y;
    arr[2] = box_z;
    return arr;
*/
    return box_dimensions;
}

void Configuration::setMirrorDepth(int _mirror_depth)
{
    mirror_depth = _mirror_depth;
}


ConfigurationRecord Configuration::recordAt(int i)
{
    return records[i];
}


void Configuration::deleteRecordAt(int i)
{
    records.erase(records.begin() + i);
}


int Configuration::getSize()
{
    return records.size();
}

int Configuration::pushBack(ConfigurationRecord record)
{
    records.push_back(record);
    return records.size();
}

void Configuration::cram()
{
    for (int i=0; i<records.size(); i++)
    {
        // check for atoms above upper bound
        while (records[i].x > box_dimensions[0]) records[i].x -= box_dimensions[0];
        while (records[i].y > box_dimensions[1]) records[i].y -= box_dimensions[1];
        while (records[i].z > box_dimensions[2]) records[i].z -= box_dimensions[2];

        // check for atoms below lower bound
        while (records[i].x < 0.0f) records[i].x += box_dimensions[0];
        while (records[i].y < 0.0f) records[i].y += box_dimensions[1];
        while (records[i].z < 0.0f) records[i].z += box_dimensions[2];
    }
    crammed = 1;
}


int Configuration::isCrammed()
{
    return crammed;
}


/*
void Configuration::replicate(int depth)
{
    // Use size of original vector
    size_t size = records.size();

    for (int r = 0; r < size; r++)
    {
        for (int i=-depth; i<=depth; i++)
        for (int j=-depth; j<=depth; j++)
        for (int k=-depth; k<=depth; k++)
        {
            // skip the center box
            if (!((i == 0) && (j == 0) && (k == 0)))
                pushBack(ConfigurationRecord((box_dimensions[0] * i) + records[r].x, 
                                             (box_dimensions[1] * j) + records[r].y, 
                                             (box_dimensions[2] * k) + records[r].z, 
                                             records[r].sigma, 
                                             records[r].epsilon)); 
        }
    }
}
*/


void Configuration::replicate(std::vector<int> depths)
{
    // Use size of original vector
    size_t size = records.size();

    for (int r = 0; r < size; r++)
    {
        for (int i=0; i<=depths[0]; i++)
        for (int j=0; j<=depths[1]; j++)
        for (int k=0; k<=depths[2]; k++)
        {
            // skip the center box
            if (!((i == 0) && (j == 0) && (k == 0)))
                pushBack(ConfigurationRecord((box_dimensions[0] * i) + records[r].x, 
                                             (box_dimensions[1] * j) + records[r].y, 
                                             (box_dimensions[2] * k) + records[r].z, 
                                             records[r].sigma, 
                                             records[r].epsilon)); 
        }
    }
}


#ifdef BUILD_PYBIND_BINDINGS

pybind11::str Configuration::__repr__()
{
    pybind11::str retval("");

    for (int i=0; i<records.size(); i++)
        retval = retval + 
             pybind11::str(std::to_string(records[i].x)) +
             pybind11::str("\t") +
             pybind11::str(std::to_string(records[i].y)) +
             pybind11::str("\t") +
             pybind11::str(std::to_string(records[i].z)) +
             pybind11::str("\t") +
             pybind11::str(std::to_string(records[i].sigma)) +
             pybind11::str("\t") +
             pybind11::str(std::to_string(records[i].epsilon)) +
             pybind11::str("\n");

    retval = retval + pybind11::str("box dims: ")
             + pybind11::str(std::to_string(box_dimensions[0]))
             + pybind11::str("\n");
    retval = retval + pybind11::str("          ")
             + pybind11::str(std::to_string(box_dimensions[1]))
             + pybind11::str("\n");
    retval = retval + pybind11::str("          ")
             + pybind11::str(std::to_string(box_dimensions[2]))
             + pybind11::str("\n");

    retval = retval + pybind11::str("temperature: ")
             + pybind11::str(std::to_string(temperature))
             + pybind11::str("\n");

    return retval;
}

#endif

