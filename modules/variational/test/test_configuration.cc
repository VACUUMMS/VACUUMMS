#include <vacuumms/variational/configuration.hh>
#include <vacuumms/types.h>
#include <iostream>

int main()
{
    char filename[] = "x.gfg";
    Configuration c = Configuration(filename);
    c.dumpContents();

    vacuumms_float x = 0.0;
    vacuumms_float y = 0.0;
    vacuumms_float z = 0.0;
    vacuumms_float sigma = 1.0;
    vacuumms_float epsilon = 1.0;

    std::cout << "insertion energy: " << c.insertionEnergy(x, y, z, sigma, epsilon) << std::endl;
}
