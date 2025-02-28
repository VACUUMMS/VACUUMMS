/* vacuumms/parameters.hh */
#pragma once

#ifdef BUILD_BOOST_PYTHON_BINDINGS 
#include <boost/python.hpp>
#endif

#include <vacuumms/exports.hh>

#include <string>
#include <vector>

#include <vacuumms/types.h>

class 
#ifdef PYBIND11_EXPORTS 
PYBIND11_EXPORT 
#endif
Parameters
{
    private:

        int parameter_argc;
        std::vector<std::string> parameter_argv;

    public:

        Parameters();
        Parameters(int argc, char **argv);
        Parameters(std::vector<std::string>);

#ifdef BUILD_PYBIND_BINDINGS 
        Parameters(const pybind11::list&);
#endif

        // Maintain a set of methods for boost::python. 
        // Note that these return the value instead of setting pointed value.

#ifdef BUILD_BOOST_PYTHON_BINDINGS 

        Parameters(const boost::python::list&);
//        int getIntParam(char* param_name);
        long getLongParam(char* param_name);
//        vacuumms_float getFloatParam(char* param_name);
        double getDoubleParam(char* param_name);
//        const char* getStringParam(char* param_name);
        boost::python::list getVectorParam(char* param_name);
        boost::python::list getVectorStringParam(char* param_name);

#endif


#ifdef BUILD_PYBIND_BINDINGS 

        pybind11::str __str__();
        pybind11::str __repr__();

//        Parameters(const boost::python::list&);

//        int getIntParam(char* param_name);
//        vacuumms_float getFloatParam(char* param_name);
//        const char* getStringParam(char* param_name);
        pybind11::list getVectorParam(char* param_name);
        pybind11::list getVectorStringParam(char* param_name);

#endif
// FTW pybind11::list getVectorParam(char* param_name);

// pulling this out of pybind/boost world because it can live without
        int getIntParam(char* param_name);
        vacuumms_float getFloatParam(char* param_name);
        const char* getStringParam(char* param_name);

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
	int getVectorStringParam(char *param_name, 
                                 const char **parameter1, 
                                 const char **parameter2, 
                                 const char **parameter3);

        int addParameter(const char* parameter);
        int getFlagParam(char *param_name);

};

