# RELEASE NOTES

---

## 1.2.0:

- Introduction of optional modules and First release of the variational module
- The variational module introduces new utilities for the variational calculation of paths of least energy between points on the insertion energy landscape. 
- introduces support for voronoi graph generation (using voro++ backend)
- VACUUMMS internals are moving toward factoring out types, constants, etc. (e.g. include/vacuumms/types.h, include/vacuumms/constants.h)
- Introduces the VACUUMMS C++ library (libvacuumms_cpp) and C/C++ header APIs. The examples in the variational module are written to use this API, and are compiled in user space using Makefiles, rather than as part of the VACUUMMS build system. 
- Refactor of core internal code to four libraries:
    - vacuumms_rt.so: The VACUUMMS runtime library includes most of the original C code used in developing VACUUMMS. 
    - vacuumms_cpp.so: The VACUUMMS C++ library contains components used by many of the more recent extensions which are written in C++, and go here.
    - vacuumms_cuda.so: Optional, contains CUDA implementations used by some applications and utilities. Built when BUILD_CUDA_COMPONENTS is set (+cuda in spack)
    - vacuumms_tiff.so: Optional, built when ENABLE_TIFF_UTILS is set (+tiff in spack)
- The first optional module 'variational' is now available. It is built when BUILD_VARIATIONAL_MODULE is set (+variational in spack). It uses the newly exposed header APIs and libraries and provides additional headers \<vacuuumms/variational/*\> and library vacuumms_variational.so. It also contains developer examples with Makefiles, so that users can see how to build code against the functionality that the module provides, without having to touch the VACUUMMS build system. 
- Basic testing via CMake/CTest is introduced and enabled by default, with unit tests placed in relevant source directories.

Because of the available API and libraries, VACUUMMS is now no longer purely a top-level code and developers are welcome to use or extend any of the available functionality.

## 1.1.4: 

- Increase system size limit from 65536 atoms to 131072 atoms (largest data structure is still <~ 16MB, so no memory issues anticipated)
- add FCC example
- Simplify some of the examples (polystyrene)

## 1.1.3: Fix TIFF issues

- Large TIFF issue addressed (integer type was overflowing for a tiff sized at 4 x 1024 x 1024 x 1024, resulting in a malloc of size 0, which was successful, but useless)
- Fixed non-working FVI2TIFF. 
- Added ability to map to more than one channel in the TIFF
- Change hard-coded limits to allow for larger configurations in gfg2fvi. Extracted hard-coded values to new limits.h file.

## 1.1.2

- Now builds with apple-clang, both x86_64 and arm64
- Add examples
- Remove redundant definitions/includes that were confounding builds

## 1.1.1

- Now builds with linux, gcc for both x86_64 and arm64
- Add explicit dependency on X11
- Remove pointer-to-stack-variable bug that would *sometimes* compile and/or run OK
- Purge the malloc-wrapped-with-an-assertion that was causing release builds to segfault, while debug ran fine

## 1.1.0

- Add first pass at documentation and user guides
- Add first few working utilities

## 1.0.0 

- Introduce first pass at CMake build system generator
