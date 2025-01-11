
#include <boost/python.hpp>
#include <iostream>
#include <string>
#include <vector>

#include <vacuumms/parameters.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>

namespace bp = boost::python;

BOOST_PYTHON_MODULE(vacuumms) 
{
    // Thin wrappers for overleaded functions
     
    int (Parameters::*p_getIntParam)(char*) = &Parameters::getIntParam;
    long (Parameters::*p_getLongParam)(char*) = &Parameters::getLongParam;
    float (Parameters::*p_getFloatParam)(char*) = &Parameters::getFloatParam;
    double (Parameters::*p_getDoubleParam)(char*) = &Parameters::getDoubleParam;
    const char* (Parameters::*p_getStringParam)(char*) = &Parameters::getStringParam;
    boost::python::list (Parameters::*p_getVectorParam)(char*) = &Parameters::getVectorParam;
    boost::python::list (Parameters::*p_getVectorStringParam)(char*) = &Parameters::getVectorStringParam;

    bp::class_<Parameters>("Parameters")
        .def(bp::init<bp::list>())
        .def("getIntParam", p_getIntParam)
        .def("getLongParam", p_getLongParam)
        .def("getFloatParam", p_getFloatParam)
        .def("getDoubleParam", p_getDoubleParam)
        .def("getStringParam", p_getStringParam)
        .def("getVectorParam", p_getVectorParam)
        .def("getVectorStringParam", p_getVectorStringParam)

        // These don't require the wrapper since they aren't overloaded
	// and have same function signature with or without boost
        .def("addParameter", &Parameters::addParameter)
        .def("getFlagParam", &Parameters::getFlagParam);

    bp::class_<Configuration>("Configuration", bp::init<char*>())
    ;

    bp::class_<CavityConfiguration>("CavityConfiguration", bp::init<char*>())
    ;

    bp::class_<CavitySizeDistribution>("CavitySizeDistribution", bp::init<CavityConfiguration>())
        .def("print", &CavitySizeDistribution::print)
    ;
}

