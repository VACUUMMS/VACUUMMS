/* vacuumms/fvi.hh */

#pragma once

#include <vacuumms/operations.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/parameters.hh>
#include <vacuumms/types.hh>

#include <vacuumms/exports.hh>


class
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
FVIX : public Operation
{
    public:
        
        FVIX();
        FVIX(Configuration, Parameters);
        void printUsage();
        void setParameters(Parameters);
        Parameters getParameters();
        void setConfiguration(Configuration);
        Configuration getConfiguration();
        void execute();
        void* getResult();
        void printResult();

        template<size_t resolution> FVIArray<resolution>* calculateFVI(Configuration);

        // wrapper to CUDA kernel
        vacuumms_EnergyArray16* calculateRepulsions(Configuration gfg);

        ~FVIX();

#ifdef BUILD_PYBIND_BINDINGS
        pybind11::str __repr__();
#endif

    private:
        
        Parameters p;
        Configuration c;
        vacuumms_EnergyArray16* ea;

        int resolution = 16;
        float attenuator = 1.0;
        float preexponential = 1.0;
        float sigma=0.0;
        float epsilon=1.0;
        float temperature = 1.0;

};
 






