
#include <boost/python.hpp>
#include <iostream>
#include <string>
#include <vector>

#include <vacuumms/parameters.hh>

namespace bp = boost::python;

BOOST_PYTHON_MODULE(vacuumms) 
{
    // Thin wrappers for overleaded functions
// these two can be changed to char*
    int (Parameters::*p_getIntParam)(boost::python::object obj) = &Parameters::getIntParam;
    long (Parameters::*p_getLongParam)(boost::python::object obj) = &Parameters::getLongParam;
     
    float (Parameters::*p_getFloatParam)(char*) = &Parameters::getFloatParam;
    double (Parameters::*p_getDoubleParam)(char*) = &Parameters::getDoubleParam;
    const char* (Parameters::*p_getStringParam)(char*) = &Parameters::getStringParam;
    //std::vector<double> (Parameters::*p_getVectorParam)(char*) = &Parameters::getVectorParam;
    boost::python::list (Parameters::*p_getVectorParam)(char*) = &Parameters::getVectorParam;
    //std::vector<std::string> (Parameters::*p_getVectorStringParam)(char*) = &Parameters::getVectorStringParam;
    boost::python::list (Parameters::*p_getVectorStringParam)(char*) = &Parameters::getVectorStringParam;

/*
        int getVectorParam(char *param_name,
                                 double *parameter1,
                                 double *parameter2,
                                 double *parameter3);
        int getVectorStringParam(char *param_name,
                                 const char **parameter1,
                                 const char **parameter2,
                                 const char **parameter3);

        int getFlagParam(char *param_name);
*/

      
    bp::class_<Parameters>("Parameters")
        .def(bp::init<bp::list>())
//        .def(bp::init<std::string>())
//        .def("Parameters", init<std::string>())
//        .def("setParameters", &Parameters::setParameters)
        .def("addParameter", &Parameters::addParameter)
        .def("getIntParam", p_getIntParam)
        .def("getLongParam", p_getLongParam)

        .def("getFloatParam", p_getFloatParam)
        .def("getDoubleParam", p_getDoubleParam)
        .def("getStringParam", p_getStringParam)
        .def("getVectorParam", p_getVectorParam)
        .def("getVectorStringParam", p_getVectorStringParam)

        .def("getFlagParam", &Parameters::getFlagParam);

}

