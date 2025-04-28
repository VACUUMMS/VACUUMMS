// operations.hh

#pragma once

#include <vacuumms/limits.h>

#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>
#include <vacuumms/parameters.hh>

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


