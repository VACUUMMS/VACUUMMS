#include <vacuumms/parameters.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>
#include <vacuumms/operations.hh>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>

/*
BOOST_PYTHON_MODULE(vacuumms)
{
    // Parameters type

    // Thin wrappers for overleaded functions

    int (Parameters::*p_getIntParam)(char*) = &Parameters::getIntParam;
    long (Parameters::*p_getLongParam)(char*) = &Parameters::getLongParam;
    float (Parameters::*p_getFloatParam)(char*) = &Parameters::getFloatParam;
    double (Parameters::*p_getDoubleParam)(char*) = &Parameters::getDoubleParam;
    const char* (Parameters::*p_getStringParam)(char*) = &Parameters::getStringParam;
    boost::python::list (Parameters::*p_getVectorParam)(char*) = &Parameters::getVectorParam;
    boost::python::list (Parameters::*p_getVectorStringParam)(char*) = &Parameters::getVectorStringParam;
}
*/

namespace py = pybind11;

// thin wrappers, maybe more readable than the lambdas?
//  int (Parameters::*p_getIntParam)(char*) = &Parameters::getIntParam;
//  py::list (Parameters::*p_getVectorParam)(char*) = &Parameters::getVectorParam;

PYBIND11_MODULE(vacuumms, m)
{
    // Declare a python wrapper and expose member functions
    py::class_<Parameters>(m, "Parameters")
        .def(py::init<>())
        .def("addParameter", &Parameters::addParameter)
        .def("getIntParam", [](Parameters& self, char* arg) -> int {return self.getIntParam(arg);} )
        .def("getVectorParam", [](Parameters& self, char* arg)-> py::list {return self.getVectorParam(arg); })
        ;

    // Same thing as above, implemented using thin wrappers.
    // .def("getIntParam", p_getIntParam)
    // .def("getVectorParam", &Parameters::getVectorParam)
    // .def("getVectorParam", p_getVectorParam)

    // This one won't work because getIntParam is overloaded
    // .def("getIntParam", &Parameters::getIntParam)

    // Can declare a module function outside of a class like this
    // m.def("getlist", &getlist);
}
