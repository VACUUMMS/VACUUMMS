#include <vacuumms/exports.hh>

#include <vacuumms/parameters.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>
#include <vacuumms/operations.hh>
#include <vacuumms/ddx.hh>
#include <vacuumms/pddx.hh>
#include <vacuumms/lammps.hh>
#include <vacuumms/voronoi.hh>

#include <vacuumms/types.h>

//#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
#include <iostream>


namespace py = pybind11;

PYBIND11_MODULE(vacuumms, m)
{
    // Declare a python wrapper and expose member functions for Parameters class
    py::class_<Parameters>(m, "Parameters")
        .def(py::init<>())
        .def(py::init<py::list>())
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

    // Configuration type

    py::class_<Configuration>(m, "Configuration")
        .def(py::init<char*>())
        .def("__repr__", &Configuration::__repr__)
        .def("setBoxDimensions", &Configuration::setBoxDimensions)
        .def("cram", &Configuration::cram)
        .def("isCrammed", &Configuration::isCrammed)
        ;

    py::class_<LAMMPSConfiguration>(m, "LAMMPSConfiguration")
        .def(py::init<std::string>())
        .def("getSize", &LAMMPSConfiguration::getSize)
        .def("__repr__", &LAMMPSConfiguration::__repr__)
        ;


    // CavityConfiguration type

    py::class_<CavityConfiguration>(m, "CavityConfiguration")
        .def(py::init<char*>())
        .def("__repr__", &CavityConfiguration::__repr__)
        ;


    // Operations classes

    // Interface to DDX (Operation subclass)

    py::class_<DDX>(m, "DDX")
        .def(py::init<>())
        .def(py::init<Configuration, Parameters>())
        .def("printUsage", &DDX::printUsage)
        .def("setParameters", &DDX::setParameters)
        .def("setConfiguration", &DDX::setConfiguration)
        .def("execute", &DDX::execute)
        .def("getResult", &DDX::getResult)
        .def("__repr__", &DDX::__repr__)
    ;

    // Interface to PDDX (Operation subclass)
    
    py::class_<PDDX>(m, "PDDX")
        .def(py::init<>())
        .def(py::init<Configuration, Parameters>())
        .def("printUsage", &PDDX::printUsage)
        .def("setParameters", &PDDX::setParameters)
        .def("setConfiguration", &PDDX::setConfiguration)
        .def("execute", &PDDX::execute)
        .def("getResult", &PDDX::getResult)
        .def("__repr__", &PDDX::__repr__)
    ;

    // Interface to Voronoi (Operation subclass)
    
    py::class_<Voronoi>(m, "Voronoi")
        .def(py::init<>())
        .def(py::init<Configuration, Parameters>())
        .def("getVertices", &Voronoi::getVertices)
        .def("getEdges", &Voronoi::getEdges)
//        .def("printUsage", &Voronoi::printUsage)
//        .def("setParameters", &Voronoi::setParameters)
//        .def("setConfiguration", &Voronoi::setConfiguration)
//        .def("execute", &Voronoi::execute)
//        .def("getResult", &Voronoi::getResult)
//        .def("__repr__", &Voronoi::__repr__)
    ;

    // Declare the Vertex and Edge classes so they can be mapped in python

    py::class_<VoronoiVertex>(m, "VoronoiVertex")
        .def(py::init<>())
        .def_readwrite("x", &VoronoiVertex::x)
        .def_readwrite("y", &VoronoiVertex::y)
        .def_readwrite("z", &VoronoiVertex::z)
        .def("__repr__", &VoronoiVertex::__repr__)
    ;

    py::class_<VoronoiEdge>(m, "VoronoiEdge")
        .def(py::init<>())
        .def_readwrite("v1", &VoronoiEdge::v1)
        .def_readwrite("v2", &VoronoiEdge::v2)
        .def("__repr__", &VoronoiEdge::__repr__)
    ;

    // Other classes
    
    // Interface to CSD (Histogram subclass)

    py::class_<CavitySizeDistribution>(m, "CavitySizeDistribution")
        .def(py::init<CavityConfiguration, Parameters>())
        .def("setWeightingExponent", &Histogram::setWeightingExponent)
        .def("print", &Histogram::print)
        .def("normalize", &Histogram::normalize)
        .def("smooth", &Histogram::smooth)
        .def("writeToFile", &Histogram::writeToFile)
        .def("__repr__", &CavitySizeDistribution::__repr__)
    ;


    // Histogram type
    
    py::class_<Histogram>(m, "Histogram")
        .def(py::init<>())
        .def(py::init<int, vacuumms_float>())
        .def("bin", &Histogram::bin)
        .def("getMisses", &Histogram::getMisses)        
        .def("writeToFile", &Histogram::writeToFile)
        .def("setWeightingExponent", &Histogram::setWeightingExponent)
        .def("__repr__", &Histogram::__repr__)
        ;


} // end of bindings 
