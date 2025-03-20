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
//        virtual void printUsage() = 0;
        static void printUsage();

        void setParameters(Parameters _p)
//        virtual void setParameters(Parameters _p)
        {
  printf("setting params\n");
            p = _p;
        }

        virtual Parameters getParameters() 
        { 
  printf("getting params\n");
            return p;
        }

    private:

        Parameters p;
};


