/* vacuumms/parameters.hh */

#pragma once

#include <string>
#include <vector>

#include <vacuumms/types.h>

#include <vacuumms/exports.hh>


#ifdef PYBIND11_EXPORTS 
class PYBIND11_EXPORT Parameters
#else
class Parameters
#endif
{
    private:

        int parameter_argc;
        std::vector<std::string> parameter_argv;

    public:

        Parameters();
        Parameters(int argc, char **argv);
        Parameters(std::vector<std::string>);
        Parameters(const char* filename);

#ifdef BUILD_PYBIND_BINDINGS 

        Parameters(const pybind11::list&);

        pybind11::str __str__();
        pybind11::str __repr__();

        pybind11::list getVectorStringParam(char* param_name);

#endif


        // pulling these out of pybind/boost world because it can live without
       
        int getIntParam(const char* param_name);
        int getIntParam(char* param_name);
        int getIntParam(std::string param_name);
        vacuumms_float getFloatParam(char* param_name);
        const char* getStringParam(char* param_name);
        std::vector<vacuumms_float> getVectorParam(std::string);
        std::vector<vacuumms_float> getVectorParam(const char*);
        std::vector<vacuumms_float> getVectorParam(char*);

        /* if a parameter is received, return a true value, otherwise return NULL) */
        int getIntParam(char *param_name, int *parameter);
        int getLongParam(char *param_name, long *parameter);
        int getFloatParam(char *param_name, vacuumms_float *parameter);
        int getDoubleParam(char *param_name, double *parameter);
        int getStringParam(char *param_name, const char **parameter);
        int getVectorParam(char *param_name, 
                                 double *parameter1,  
                                 double *parameter2, 
                                 double *parameter3);
        int getVectorParam(char *param_name, 
                                 vacuumms_float *parameter1,  
                                 vacuumms_float *parameter2, 
                                 vacuumms_float *parameter3);
        int getVectorStringParam(char *param_name, 
                                 const char **parameter1, 
                                 const char **parameter2, 
                                 const char **parameter3);

        int addParameter(const char* parameter);
        int getFlagParam(char *param_name);
        int toFile(const char* filename);

};

