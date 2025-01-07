
#include <boost/python.hpp>
#include <iostream>

#include <vacuumms/parameters.hh>

using namespace boost::python;

BOOST_PYTHON_MODULE(vacuumms) 
{

    class_<Parameters>("Parameters", init<int, char**>())
        .def("getIntParam", &Parameters::getIntParam) ;

}

