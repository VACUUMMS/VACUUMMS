/* vacuumms/fvi.hh */

#pragma once

#include <vacuumms/operations.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/parameters.hh>
#include <vacuumms/types.hh>

#include <vacuumms/exports.hh>

/* Mask bits for internal API */
#define FVIX_ATTRACTION 1
#define FVIX_REPULSION 2
#define FVIX_ENERGY 4
#define FVIX_FVI 8


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
        void* getResult();

        void execute();
        void executeMask(int mask);

        void printUsage();

        void calculateAll();
        std::vector<vacuumms_float> calculateAttraction();
        std::vector<vacuumms_float> calculateRepulsion();
        std::vector<vacuumms_float> calculateEnergy();
        std::vector<vacuumms_float> calculateFVI();
//        template<size_t resolution> FVIArray<resolution>* calculateFVI(Configuration);

        ~FVIX();

#ifdef BUILD_TIFF_UTILS
        void generateTIFF(char*);
#endif

#ifdef BUILD_PYBIND_BINDINGS
        pybind11::array_t<vacuumms_float>getAttraction();
        pybind11::array_t<vacuumms_float>getRepulsion();
        pybind11::array_t<vacuumms_float>getEnergy();
        pybind11::array_t<vacuumms_float>getFVI();
        pybind11::str __repr__();
#endif

    private:
        
        Parameters p;
        Configuration c;
        
        std::vector<size_t> dimensions = {2,2,2};
        std::vector<vacuumms_float> attraction;
        std::vector<vacuumms_float> repulsion;
        std::vector<vacuumms_float> energy;
        std::vector<vacuumms_float> FVI;

        int resolution = 16;
        vacuumms_float attenuator = 1.0;
        vacuumms_float preexponential = 1.0;
        vacuumms_float sigma=0.0;
        vacuumms_float epsilon=1.0;
        vacuumms_float temperature = 1.0;

};
 






