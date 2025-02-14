#include <vacuumms/parameters.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>
#include <vacuumms/operations.hh>

#include <pybind11/pybind11.h>
#include <iostream>

//namespace bp = boost::python;

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

/*
namespace lib
{
    int f();
    int g(int a);
}
*/

namespace lib
{
    int f()
    {
        std::cout << "f()\n";
        return 0;
    }

int g(int a)
    {
        std::cout << "g(int)\n";
        std::cout << a << "\n";
        return 1;
    }
}

namespace py = pybind11;

PYBIND11_MODULE(vacuumms, m)
{
    m.def("f", &lib::f);
    m.def("g", &lib::g);

//    m.def("Parameters", &Parameters::Parameters);


py::class_<Parameters>(m, "Parameters")
    .def(py::init<>());
}
