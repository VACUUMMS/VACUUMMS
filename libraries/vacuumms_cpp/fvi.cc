/* libraries/vacuumms_cpp/fvi.cc */

#include <vacuumms/ddx.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>
#include <vacuumms/fvi.hh>

#include <vacuumms/limits.h>
#include <vacuumms/rng.h>

#include <math.h>


FVIX::FVIX(Configuration c, Parameters p) : 
    c{c}, p{p} 
{
    setParameters(p);
}

FVIX::FVIX()
{
}

FVIX::FVIX(Configuration _c)
{
    c = _c;
}


void FVIX::setParameters(Parameters _p)
{
    p = _p;

    //p.getIntParam((char*)"-resolution", &resolution);
    resolution = p.getIntParam("-resolution");

    p.getFloatParam((char*)"-attenuator", &attenuator);
    p.getFloatParam((char*)"-preexponential", &preexponential);
    p.getFloatParam((char*)"-sigma", &sigma);
    p.getFloatParam((char*)"-epsilon", &epsilon);
    p.getFloatParam((char*)"-temperature", &temperature);

    if (p.getFlagParam((char*)"-usage")) printUsage();

    p.getVectorParam((char*)"-box", &c.box_x, &c.box_y, &c.box_z);
}

Parameters FVIX::getParameters()
{
    return p;
}

void FVIX::setConfiguration(Configuration _c)
{
    c = _c;
}

Configuration FVIX::getConfiguration()
{
    return c;
}

void FVIX::setDimensions(std::vector<size_t> _dimensions)
{
    dimensions = _dimensions;
}

std::vector<size_t> FVIX::getDimensions()
{
    return dimensions;
}

#ifdef BUILD_PYBIND_BINDINGS
pybind11::array_t<vacuumms_float>FVIX::getRepulsion()
{
    return pybind11::array_t<vacuumms_float>(dimensions, repulsion.data());
}
#endif


void FVIX::execute()
{
/*
    p.getIntParam((char*)"-resolution", &resolution);
    p.getFloatParam((char*)"-attenuator", &attenuator);
    p.getFloatParam((char*)"-preexponential", &preexponential);
    p.getFloatParam((char*)"-sigma", &sigma);
    p.getFloatParam((char*)"-epsilon", &epsilon);
    p.getFloatParam((char*)"-temperature", &temperature);

    if (p.getFlagParam((char*)"-usage")) printUsage();

    p.getVectorParam((char*)"-box", &c.box_x, &c.box_y, &c.box_z);
*/
//    ea = calculateRepulsions(c);

    runKernel();

} // end FVI::execute()

void FVIX::printResult()
{
    for (int i=0; i<resolution; i++) 
    for (int j=0; j<resolution; j++) 
    for (int k=0; k<resolution; k++)
        printf("%f\t%f\t%f\t%f\n", 
               i*c.box_x / resolution, 
               j*c.box_y / resolution, 
               k*c.box_z / resolution, 
               preexponential * exp(ea->energy[i][j][k]/(-temperature * attenuator))); 
}

//CavityConfiguration DDX::getResult()
// Need to think about how to return the huge result...
// maybe getTIFF()?

void* FVIX::getResult()
{
    void* result;
    return result;   
}


void FVIX::printUsage()
{
    printf("usage:     gfg2fvi        -box [10.0 10.0 10.0]\n");
    printf("                          -potential [69] \n");
    printf("                          -resolution [256] \n");
    printf("                          -temperature [298.0] \n");
    printf("                          -attenuator [1.0] \n");
    printf("                          -preexponential [1.0] \n");
    printf("                          -sigma [0.0] \n");
    printf("                          -epsilon [1.0] \n");
    printf("                          -check_device \n");
}

#ifdef BUILD_PYBIND_BINDINGS
pybind11::str FVIX::__repr__()
{
    pybind11::str retval;
    retval += c.__repr__();
    retval += p.__repr__();
//    retval += result.__repr__();
    return retval;
}
#endif

FVIX::~FVIX()
{
    free(ea);
}
