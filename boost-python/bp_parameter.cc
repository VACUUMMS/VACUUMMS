
#include <boost/python.hpp>
#include <iostream>

#include <vacuumms/parameters.hh>

//using namespace boost::python;
namespace bp = boost::python;

BOOST_PYTHON_MODULE(vacuumms) 
{

/* This is correct to make argc, argv type constructor, but I don't know how to get the right types to pass
    class_<Parameters>("Parameters", init<int, char**>())
*/

//    class_<Parameters>("Parameters", init<std::string>())

    bp::class_<Parameters>("Parameters")
        .def(bp::init<bp::list>())
//        .def(bp::init<std::string>())
//        .def("Parameters", init<std::string>())
//        .def("setParameters", &Parameters::setParameters)
        .def("getIntParam", &Parameters::getIntParam) ;

}

