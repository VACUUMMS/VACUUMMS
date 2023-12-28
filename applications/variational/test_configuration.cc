#include "configuration.hh"
#include <iostream>

int main()
{
    char filename[] = "x.gfg";
    Configuration c = Configuration(filename);
    c.dumpContents();

    std::cout << "insertion energy: " << c.insertionEnergy(0.0, 0.0, 0.0) << std::endl;
}
