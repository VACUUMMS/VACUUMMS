/* vacuumms/exports.hh */

/* Declare macros so that symbols are exported/external 
 * when building PYBIND interface
 */

#pragma once

#ifdef BUILD_PYBIND_BINDINGS 
    #include <pybind11/pybind11.h>
    #ifdef PYBIND11_EXPORTS 
        #define PYBIND11_EXPORT __attribute__((visibility("default")))
    #else 
        #define PYBIND11_EXPORT
    #endif
#endif


