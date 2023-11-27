#include <iostream>
#include <cmath>
#include "configuration.hh"

ConfigurationRecord::ConfigurationRecord(float _x, float _y, float _z, float _sigma, float _epsilon)
{
    x = _x;
    y = _y;
    z = _z;
    sigma = _sigma;
    epsilon = _epsilon;
    sigma_6 = powf(sigma, (1.0/6.0));
    sigma_12 = sqrt(sigma_6);

    /* FTW: will add this in later, and include in energy calcs also
    int mirror_depth_x = 0;
    int mirror_depth_y = 0;
    int mirror_depth_z = 0;
    */
}

Configuration::Configuration()
{
    // create the records object but don't populate
    records = std::vector<ConfigurationRecord>();    
}

Configuration::Configuration(char *filename)
{
    FILE* infile = fopen(filename, "r");
    //std::cout << "constructing Configuration from " << filename << " as " << infile << std::endl;
    float x, y, z, sigma, epsilon;
    records = std::vector<ConfigurationRecord>();    
    //std::cout << "constructing Configuration.records as " << this << "." << &records << std::endl;

    while (!feof(infile))
    {
        fscanf(infile, "%f\t%f\t%f\t%f\t%f\n", &x, &y, &z, &sigma, &epsilon);
        records.push_back(ConfigurationRecord(x, y, z, sigma, epsilon));
    }
}

void Configuration::dumpContents()
{
    for (int i = 0; i < records.size(); i++)
        printf("%f\t%f\t%f\t%f\t%f\n", records[i].x, records[i].y, records[i].z, records[i].sigma, records[i].epsilon);
}

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

float Configuration::insertionEnergy3D(float x, float y, float z)
{
    float total = 0.0;
    for (int i=0; i<records.size(); i++)
    {
        float r_sq = (records[i].x - x) * (records[i].x - x) 
                   + (records[i].y - y) * (records[i].y - y)
                   + (records[i].z - z) * (records[i].z - z);
        float r_6 = r_sq * r_sq * r_sq;
        float r_12 = r_6 * r_6;
        total += (1.0/r_12 - 1.0/r_6);
    }
    return total;
}

