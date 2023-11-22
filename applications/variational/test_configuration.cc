#include "configuration.hh"
#include <iostream>

int main()
{
    char filename[] = "x.gfg";
    Configuration c = Configuration(filename);
//    std::cout << c << std::endl;
    c.dumpContents();

    std::cout << "insertion energy: " << c.insertionEnergy2D(0.0,0.0) << std::endl;
}
