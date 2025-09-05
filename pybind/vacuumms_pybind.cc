/* pybind/vacuumms_pybind.cc */

//#include <iostream>

#include <vacuumms/types.h>
#include <vacuumms/exports.hh>

#include <vacuumms/parameters.hh>
#include <vacuumms/configuration.hh>
#include <vacuumms/cavity.hh>
#include <vacuumms/operations.hh>
#include <vacuumms/ddx.hh>
#include <vacuumms/pddx.hh>
#include <vacuumms/lammps.hh>
#include <vacuumms/scene.hh>

#ifdef BUILD_VORONOI_UTILS
#include <vacuumms/voronoi.hh>
#endif

#ifdef BUILD_CUDA_COMPONENTS
#include <cuda_runtime.h>
#endif

#include <vacuumms/fvi.hh>

#ifdef BUILD_CUDA_COMPONENTS
void finalize_cuda()
{
    cudaDeviceReset();
}
#endif


namespace py = pybind11;

PYBIND11_MODULE(vacuumms, m)
{

    // Declare a python wrapper and expose member functions for Parameters class
    py::class_<Parameters>(m, "Parameters")
        .def(py::init<>())
        .def(py::init<py::list>())
        .def(py::init<const char*>())
        .def("addParameter", &Parameters::addParameter)
        .def("removeParameter", &Parameters::removeParameter)
        .def("getFlagParam", &Parameters::getFlagParam)
        .def("toFile", &Parameters::toFile)
        .def("getIntParam", [](Parameters& self, char* arg) -> int {return self.getIntParam(arg);} )
        .def("getFloatParam", [](Parameters& self, char* arg) -> vacuumms_float {return self.getFloatParam(arg);})
        .def("getStringParam", [](Parameters& self, char* arg) -> py::str {return self.getStringParam(arg);})
        .def("getVectorParam", [](Parameters& self, char* arg)-> std::vector<vacuumms_float> {return self.getVectorParam(arg); })
        .def("getVectorStringParam", [](Parameters& self, char* arg)-> py::list {return self.getVectorStringParam(arg); })
        .def("__repr__", &Parameters::__repr__)
        .def("__str__", &Parameters::__str__)
        ;


    // Configuration type(s)

    py::class_<Configuration>(m, "Configuration")
        .def(py::init<char*>())
        .def(py::init<char*, std::vector<vacuumms_float>>())
        .def(py::init<char*, std::vector<vacuumms_float>, vacuumms_float>())
        .def("__repr__", &Configuration::__repr__)
        .def("setBoxDimensions", &Configuration::setBoxDimensions)
        .def("getBoxDimensions", &Configuration::getBoxDimensions)
        .def("setTemperature", &Configuration::setTemperature)
        .def("getTemperature", &Configuration::getTemperature)
        .def("cram", &Configuration::cram)
        .def("isCrammed", &Configuration::isCrammed)
        .def("replicate", &Configuration::replicate)
        .def("getSize", &Configuration::getSize)
        ;

    py::class_<LAMMPSConfiguration, Configuration>(m, "LAMMPSConfiguration")
        .def(py::init<std::string>())
        ;


    // CavityConfiguration type

    py::class_<CavityConfiguration>(m, "CavityConfiguration")
        .def(py::init<char*>())
        .def("scrubDuplicates", &CavityConfiguration::scrubDuplicates)
        .def("__repr__", &CavityConfiguration::__repr__)
        ;

    /* Operations subclass notes: 
 
       Declare any operation that will be called on all operations in 
       base class. It will use the method from derived class, even
       though it is not listed explicitly for derived class. 

       There is no Trampoline class because subclasses aren't overriding 
       base class methods in a way that requires it.
    */

    // Operations classes
    
    py::class_<Operation>(m, "Operation")
        .def("getParameters", &Operation::getParameters)
        .def("setParameters", &Operation::setParameters)
        .def("execute", &Operation::execute)
        ;

    // Interface to DDX (Operation subclass)

    py::class_<DDX, Operation>(m, "DDX")
        .def(py::init<>())
        .def(py::init<Configuration>())
        .def(py::init<Configuration, Parameters>())
        .def("printUsage", &DDX::printUsage)
        .def("getConfiguration", &DDX::getConfiguration)
        .def("setConfiguration", &DDX::setConfiguration)
        .def("setNumberOfSamples", &DDX::setNumberOfSamples)
        .def("setSeed", &DDX::setSeed)
        .def("setVerletCutoff", &DDX::setVerletCutoff)
        .def("setNumberOfSteps", &DDX::setNumberOfSteps)
        .def("setMinDiameter", &DDX::setMinDiameter)
        .def("randomize", &DDX::Randomize)
        .def("getResult", &DDX::getResult)
        .def("__repr__", &DDX::__repr__)
    ;

    // Interface to PDDX (Operation subclass)
    
    py::class_<PDDX, Operation>(m, "PDDX")
        .def(py::init<>())
        .def(py::init<Configuration, Parameters>())
        .def("printUsage", &PDDX::printUsage)
        .def("getConfiguration", &PDDX::getConfiguration)
        .def("setConfiguration", &PDDX::setConfiguration)
        .def("execute", &PDDX::execute)
        .def("getResult", &PDDX::getResult)
        .def("__repr__", &PDDX::__repr__)
    ;


    // Scene interface, for generating and rendering POVRay SDL
    
    py::class_<Scene>(m, "Scene")
        .def(py::init<>())
        .def("createSceneFile", &Scene::createSceneFile)  // POV file
        .def("renderScene", &Scene::renderScene)          // PNG file
        .def("setRenderDimensions", &Scene::setRenderDimensions)
        .def("setCameraLocation", &Scene::setCameraLocation)
        .def("setCameraLookAt", &Scene::setCameraLookAt)
        .def("applyStandardLight", &Scene::applyStandardLight)
        .def("applyAmbientLight", &Scene::applyAmbientLight)
        .def("addLightSource", &Scene::addLightSource)
        .def("clearLightSources", &Scene::clearLightSources)
        .def("addSceneComponent", &Scene::addSceneComponent)
        .def("generateContainerSDL", &Scene::generateContainerSDL)
        .def("dumpSDL", &Scene::dumpSDL)

//        SceneComponent componentAt(int i);
//        size_t deleteComponentAt(int i);
//        size_t getNumberOfComponents();
    ;
    
    py::class_<SceneComponent>(m, "SceneComponent")
        .def("getComponentSDL", &SceneComponent::getComponentSDL)
        .def("setTransmit", &SceneComponent::setTransmit)
        .def("setPhong", &SceneComponent::setPhong)
        .def(py::init<>())
    ;

    py::class_<ConfigurationComponent, SceneComponent>(m, "ConfigurationComponent")
        .def(py::init<>())
        .def(py::init<Configuration>())
        .def("getComponentSDL", &ConfigurationComponent::getComponentSDL)
    ;

    py::class_<CavityComponent, SceneComponent>(m, "CavityComponent")
        .def(py::init<>())
        .def(py::init<CavityConfiguration>())
        .def("getComponentSDL", &CavityComponent::getComponentSDL)
    ;

#ifdef BUILD_VORONOI_UTILS

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

#endif


    py::class_<FVIX>(m, "FVIX")
        .def(py::init<>())
        .def(py::init<Configuration>())
        .def(py::init<Configuration, Parameters>())
        .def("printUsage", &FVIX::printUsage)
        .def("getParameters", &FVIX::getParameters)
        .def("setParameters", &FVIX::setParameters)
        .def("getDimensions", &FVIX::getDimensions)
        .def("setDimensions", &FVIX::setDimensions)
        .def("getConfiguration", &FVIX::getConfiguration)
        .def("setConfiguration", &FVIX::setConfiguration)
        .def("execute", &FVIX::execute)
        .def("calculateAttraction", &FVIX::calculateAttraction)
        .def("calculateRepulsion", &FVIX::calculateRepulsion)
        .def("calculateEnergy", &FVIX::calculateEnergy)
        .def("calculateFVI", &FVIX::calculateFVI)
        .def("getRepulsion", &FVIX::getRepulsion)
        .def("getAttraction", &FVIX::getAttraction)
        .def("getEnergy", &FVIX::getEnergy)
        .def("getFVI", &FVIX::getFVI)
#ifdef BUILD_TIFF_UTILS
        .def("generateTIFF", &FVIX::generateTIFF)
#endif
        .def("__repr__", &FVIX::__repr__)
    ;

#ifdef BUILD_CUDA_COMPONENTS
    m.def("finalize_cuda", &finalize_cuda)
#endif
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
