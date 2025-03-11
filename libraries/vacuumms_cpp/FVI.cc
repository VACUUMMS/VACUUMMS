/* libraries/vacuumms_cpp/FVI.cc */

/*
#include <vacuumms/config_parser.h>
#include <vacuumms/param.h>
#include <vacuumms/types.h>
#include <vacuumms/gfg2fvi.h>

#include <stdio.h>
#include <math.h>
*/

//float attenuator = 1.0;

#include <vacuumms/ddx.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>

#include <vacuumms/limits.h>
#include <vacuumms/rng.h>

#include <math.h>


FVI::FVI(Configuration c, Parameters p) : 
    c{c}, p{p} 
{
}

FVI::FVI()
{
}


void FVI::setParameters(Parameters _p)
{
    p = _p;
}

void FVI::setConfiguration(Configuration _c)
{
    c = _c;
}


void FVI::execute()
{
/*
  int i,j,k;
  double box_x=10, box_y=10, box_z=10;
  int potential = 69;
  int resolution = 256;
  float temperature = 298.0;
*/
  float attenuator = 1.0;
  float preexponential = 1.0;
  float sigma=0.0;
  float epsilon=1.0;

  p.getFloatParam((char*)"-attenuator", &attenuator);
  p.getFloatParam((char*)"-preexponential", &preexponential);
  p.getFloatParam((char*)"-sigma", &sigma);
  p.getFloatParam((char*)"-epsilon", &epsilon);

  if (p.getFlagParam((char*)"-usage")) printUsage();


/* this is all implicit in Configuration object now 
  vacuumms_GFG65536 *gfg = readGFG65536(stdin);
  gfg->box_x = box_x;
  gfg->box_y = box_y;
  gfg->box_z = box_z;
*/

// we only care about this case:

    vacuumms_EnergyArray256 *ea = GFGToRepulsion256_612(gfg, sigma, epsilon);
    for (i=0; i<resolution; i++) 
    for (j=0; j<resolution; j++) 
    for (k=0; k<resolution; k++)
        printf("%f\t%f\t%f\t%f\n", 
                i*box_x / resolution, 
                j*box_y / resolution, 
                k*box_z / resolution, 
                preexponential * exp(ea->energy[i][j][k]/(-temperature * attenuator))); 

} // end FVI::execute()

//CavityConfiguration DDX::getResult()
// Need to think about how to return the huge result...
// maybe getTIFF()?

void* FVI::getResult()
{
    return result;   
}

void DDX::execute()
{

// was  loadConfiguration(); 
// now just copy over the records and run the old algorithm
  for (int i=0; i<c.getSize(); i++) 
  {
    ConfigurationRecord r = c.recordAt(i);
    x[i] = r.x;
    y[i] = r.y;
    z[i] = r.z;
    sigma[i] = r.sigma;
    epsilon[i] = r.epsilon;
  }
  
  number_of_molecules = c.getSize();

} // end DDX::execute()


//void DDX::findEnergyMinimum() { } // end DDX::findEnergyMinimum()
//double DDX::calculateRepulsion() { } // end DDX::calculateRepulsion()
//double DDX::calculateEnergy(double test_diameter) { } // end DDX::calculateEnergy()
//void DDX::expandTestParticle() { }

void FVI::printUsage()
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
pybind11::str DDX::__repr__()
{
    pybind11::str retval;
    retval += c.__repr__();
    retval += p.__repr__();
    retval += result.__repr__();
    return retval;
}
#endif

