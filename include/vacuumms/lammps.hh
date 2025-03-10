// lammps.hh 
#pragma once

#include <vector>
#include <map>

#include <stdio.h>
#include <vacuumms/types.h>

#include <vacuumms/exports.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/pair.hh>

#include <vacuumms/exports.hh>


class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
LAMMPSConfiguration : public Configuration
{
    public: 

        LAMMPSConfiguration(std::string filename);
      
#ifdef BUILD_PYBIND_BINDINGS
        pybind11::str __repr__();
#endif

    private:

        std::map<int, PairCoefficient> pairs;
// replace with std::vector<ConfigurationRecord>
//        std::vector<Atom> atoms;

        vacuumms_float xlo, xhi;
        vacuumms_float ylo, yhi;
        vacuumms_float zlo, zhi;

};
