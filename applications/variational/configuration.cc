#include <iostream>
#include "configuration.hh"

ConfigurationRecord::ConfigurationRecord(float _x, float _y, float _z, float _sigma, float _epsilon)
{
    x = _x;
    y = _y;
    z = _z;
    sigma = _sigma;
    epsilon = _epsilon;
}

Configuration::Configuration(char *filename)
{
    FILE* infile = fopen(filename, "r");
    float x, y, z, sigma, epsilon;

    while (!feof(infile))
    {
        fscanf(infile, "%f\t%f\t%f\t%f\t%f\n", &x, &y, &z, &sigma, &epsilon);
        push_back(ConfigurationRecord(x, y, z, sigma, epsilon));
    }

}
