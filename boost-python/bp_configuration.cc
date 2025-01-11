
#include <boost/python.hpp>
#include <iostream>
#include <string>
#include <vector>

#include <vacuumms/configuration.hh>

namespace bp = boost::python;

BOOST_PYTHON_MODULE(vacuumms) 
{
    bp::class_<Configuration>("Configuration", bp::init<char*>())
//    bp::class_<Configuration>("Configuration")
//        .def(bp::init<std::string>())
    ;
}

