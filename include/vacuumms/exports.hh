/* vacuumms/exports.hh */

#pragma once

/* Declare macros so that symbols are exported/external 
 * when building PYBIND interface
 */

#ifdef BUILD_PYBIND_BINDINGS 
    #include <pybind11/pybind11.h>
    #include <pybind11/stl.h>
    #include <pybind11/numpy.h>
    #define PYBIND11_EXPORT __attribute__((visibility("default")))
#endif
