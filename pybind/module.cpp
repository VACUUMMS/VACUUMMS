#include <pybind11/pybind11.h>
#include "lib.h"


//PYBIND11_MODULE(binding_module, m)
PYBIND11_MODULE(vacuumms, m)
{
    m.def("f", &lib::f);
    m.def("g", &lib::g);
}
