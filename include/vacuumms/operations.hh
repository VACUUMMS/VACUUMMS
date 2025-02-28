// operations.hh
#pragma once

#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>
#include <vacuumms/parameters.hh>

#include <vacuumms/limits.h>
#define MAX_CLOSE (VACUUMMS_MAX_NUMBER_OF_MOLECULES)

#ifdef BUILD_PYBIND_BINDINGS 
    #include <pybind11/pybind11.h>
    #ifdef PYBIND11_EXPORTS 
        #define PYBIND11_EXPORT __attribute__((visibility("default")))
    #endif
#endif



class PYBIND11_EXPORT Operation
{
    public:

        virtual void execute() = 0;
//        virtual void printUsage() = 0;
        static void printUsage();
};


