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
        FVIX(Configuration);
        FVIX(Configuration, Parameters);
        void setParameters(Parameters);
        Parameters getParameters();
        void setConfiguration(Configuration);
        Configuration getConfiguration();
        void setDimensions(std::vector<size_t>);
        std::vector<size_t> getDimensions();
        void execute();
        void* getResult();
        void printResult();

        void printUsage();

        template<size_t resolution> FVIArray<resolution>* calculateFVI(Configuration);

        // wrapper to CUDA kernel
        vacuumms_EnergyArray16* calculateRepulsions(Configuration gfg);
        void runKernel();

        ~FVIX();

#ifdef BUILD_PYBIND_BINDINGS
        pybind11::array_t<vacuumms_float>getRepulsion();
        pybind11::str __repr__();
#endif

    private:
        
        Parameters p;
        Configuration c;
        vacuumms_EnergyArray16* ea;
        
        std::vector<size_t> dimensions = {2,2,2};
        std::vector<vacuumms_float> attraction;
        std::vector<vacuumms_float> repulsion;
        std::vector<vacuumms_float> energy;

        int resolution = 16;
        vacuumms_float attenuator = 1.0;
        vacuumms_float preexponential = 1.0;
        vacuumms_float sigma=0.0;
        vacuumms_float epsilon=1.0;
        vacuumms_float temperature = 1.0;

};
 






