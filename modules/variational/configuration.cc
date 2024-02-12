#include <iostream>
#include <cmath>
#include <vacuumms/variational/configuration.hh>
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

Configuration::Configuration(char *filename)
{
    FILE* infile = fopen(filename, "r");
    vacuumms_float x, y, z, sigma, epsilon;
    records = std::vector<ConfigurationRecord>();    

    while (!feof(infile))
    {
        fscanf(infile, "%f\t%f\t%f\t%f\t%f\n", &x, &y, &z, &sigma, &epsilon);
        records.push_back(ConfigurationRecord(x, y, z, sigma, epsilon));
    }
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

void Configuration::dumpContents()
{
    for (int i = 0; i < records.size(); i++)
        printf("%f\t%f\t%f\t%f\t%f\n", records[i].x, records[i].y, records[i].z, records[i].sigma, records[i].epsilon);
}

/*
float Configuration::insertionEnergy2D(float x, float y)
{
printf("Calculating insertion energy at %0.012f, %0.012f\n", x,y);
    float total = 0.0;
    for (int i=0; i<records.size(); i++)
    {
        float r_sq = (records[i].x - x) * (records[i].x - x) 
                   + (records[i].y - y) * (records[i].y - y);
        float r_6 = r_sq * r_sq * r_sq;
        float r_12 = r_6 * r_6;
        total += 4 * records[i].epsilon * (records[i].sigma_12/r_12 - records[i].sigma_6/r_6);
    }
    return total;
}
*/

//float Configuration::insertionEnergy3D(float x, float y, float z)
vacuumms_float Configuration::insertionEnergy(vacuumms_float x, vacuumms_float y, vacuumms_float z, vacuumms_float sigma, vacuumms_float epsilon)
{
//printf("IN: sigma = %f, epsilon = %f\n", sigma, epsilon);
    vacuumms_float total = 0.0;
    for (int box_i=-mirror_depth; box_i<=mirror_depth; box_i++)
        for (int box_j=-mirror_depth; box_j<=mirror_depth; box_j++)
            for (int box_k=-mirror_depth; box_k<=mirror_depth; box_k++)
                for (int i=0; i<records.size(); i++)
                {
                    vacuumms_float test_x = box_i * box_x + records[i].x;
                    vacuumms_float test_y = box_j * box_y + records[i].y;
                    vacuumms_float test_z = box_k * box_z + records[i].z;
                    vacuumms_float r_sq = (test_x - x) * (test_x - x)
                               + (test_y - y) * (test_y - y)
                               + (test_z - z) * (test_z - z);

                    // Lorentz-Berthelot combining rules for sigma and epsilon
                    vacuumms_float sigma_ij = 0.5 * (sigma + records[i].sigma);
                    vacuumms_float sigma_sq = sigma_ij * sigma_ij;
                    vacuumms_float epsilon_ij = sqrt(epsilon * records[i].epsilon);

//printf("sigma_ij = %f, epsilon_ij = %f\n", sigma_ij, epsilon_ij);

                    vacuumms_float sigma_6 = sigma_sq * sigma_sq * sigma_sq;
                    vacuumms_float sigma_12 = sigma_6 * sigma_6;
                    vacuumms_float r_6 = r_sq * r_sq * r_sq;
                    vacuumms_float r_12 = r_6 * r_6;
                    // total += (1.0/r_12 - 1.0/r_6);

                    total += 4 * epsilon_ij * (sigma_12/r_12 - sigma_6/r_6);
                }
    return total;
}

void Configuration::setBoxDimensions(vacuumms_float _box_x, vacuumms_float _box_y, vacuumms_float _box_z)
{
    box_x = _box_x;
    box_y = _box_y;
    box_z = _box_z;
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
