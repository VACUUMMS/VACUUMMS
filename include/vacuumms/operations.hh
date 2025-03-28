// operations.hh
#pragma once

#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>
#include <vacuumms/parameters.hh>

#include <vacuumms/limits.h>
#define MAX_CLOSE (VACUUMMS_MAX_NUMBER_OF_MOLECULES)

#include <vacuumms/exports.hh>



class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
Operation
{
    public:

        virtual void execute() = 0;
        static void printUsage();

        virtual void setParameters(Parameters _p)
        {
            p = _p;
        }

        virtual Parameters getParameters() 
        { 
            return p;
        }

    private:

        Parameters p;
};


