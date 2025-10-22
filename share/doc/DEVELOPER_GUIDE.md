# VACUUMMS DEVELOPER GUIDE

## OVERVIEW

VACUUMMS development is hosted on github at: https://github.com/VACUUMMS/VACUUMMS

Contributions including PR's and feedback are welcome!

VACUUMMS is now composed of:

 - An API, expressed as C++ header files. To develop against the API, include the headers, e.g. #include \<vacuumms/types.h\> and link against the below libraries.

 - libraries:

    - vacuumms_rt.so: The VACUUMMS runtime library includes most of the original C code used in developing VACUUMMS and the original CLI. 
    - vacuumms_cpp.so: The VACUUMMS C++ library contains components used by many of the more recent extensions which are written in C++, and go here.
    - vacuumms_cuda.so: Optional, contains CUDA implementations used by some applications and utilities. Built when BUILD_CUDA_COMPONENTS is set (+cuda in spack) and depends on CUDA libraries.
    - vacuumms_tiff.so: Optional, built when ENABLE_TIFF_UTILS is set (+tiff in spack) and depends on TIFF libraries.
    - vacuumms_variational.so: Supports the variational module, now an integral part of VACUUMMS.

 - Internal applications and utilities. These are always installed when VACUUMMS is built. They use the VACUUMMS headers and libraries above. 
 - Conda package: Alternative way to install, and to make available to python/jupyter.

 - Spack package: Driver for the CMake build system.

 - CMake build system: When in doubt, look here to see how it all fits together. 

## STYLE GUIDE:

    CLI code is now deprecated and should be considered frozen. It should not be modified, even for style. Don't mess!
    When in doubt, try to adhere to standard conventions.
    Case:
    - C++ classes are in PascalCase.
    - C++ methods are in camelCase.
    - C++ data members are in snake_case.
    - C++ constants are in SCREAMING_SNAKE_CASE.
    Classes should be separated by double space in declarations. 
    Class methods should be separated by double space in implementations.
    Avoid wrapping of code. No hard rule on width, just strive for readability.
    Comment freely and frequently and feel free to leave contextual notes for yourself.
    Be consistent in use of iostream vs stdio.h. Favor iostream, and avoid mixing them, as they buffer separately. 
 
## APOCRYPHA

"Dead" code for old projects that were not completed, have not yet been ported or documented, or which have been replaced by a refactored version, but have not (yet) been abandoned. 

## EXAMPLES

If you add a feature, add an example to show how to use the feature.  Also, add a test. 

