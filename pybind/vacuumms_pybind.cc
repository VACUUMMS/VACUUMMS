#include <vacuumms/parameters.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>
#include <vacuumms/operations.hh>
#include <vacuumms/types.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>


namespace py = pybind11;

PYBIND11_MODULE(vacuumms, m)
{
    // Declare a python wrapper and expose member functions for Parameters class
    py::class_<Parameters>(m, "Parameters")
        .def(py::init<>())
        .def("addParameter", &Parameters::addParameter)
        .def("getFlagParam", &Parameters::getFlagParam)
        .def("getIntParam", [](Parameters& self, char* arg) -> int {return self.getIntParam(arg);} )
        .def("getFloatParam", [](Parameters& self, char* arg) -> vacuumms_float {return self.getFloatParam(arg);})
        .def("getStringParam", [](Parameters& self, char* arg) -> py::str {return self.getStringParam(arg);})
        .def("getVectorParam", [](Parameters& self, char* arg)-> py::list {return self.getVectorParam(arg); })
        .def("getVectorStringParam", [](Parameters& self, char* arg)-> py::list {return self.getVectorStringParam(arg); })
        .def("__repr__", &Parameters::__repr__)
        .def("__str__", &Parameters::__str__)
        ;

    // Can declare a module function outside of a class like this
    // m.def("getlist", &getlist);


    // Configuration type

    py::class_<Configuration>(m, "Configuration")
        .def(py::init<char*>())
        .def("__repr__", &Configuration::__repr__)
        ;


    // CavityConfiguration type

    py::class_<CavityConfiguration>(m, "CavityConfiguration")
        .def(py::init<char*>())
        ;


    // CavitySizeDistribution type

    py::class_<CavitySizeDistribution>(m, "CavitySizeDistribution")
        .def(py::init<CavityConfiguration>())
        .def("print", &CavitySizeDistribution::print)
    ;


    // Operations classes

    // Interface to DDX operation subclass

    py::class_<DDX>(m, "ddx")
        .def(py::init<Configuration, Parameters>())
        .def("getOutput", &DDX::getOutput)
    ;

}
