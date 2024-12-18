# This is the Developer Guide.

## OVERVIEW

VACUUMMS is now composed of:

 - An API, consisting of C and C++ header files. To develop against the API, include the headers, e.g. #include \<vacuumms/types.h\> and link against the below libraries.

 - Four libraries:

    - vacuumms_rt.so: The VACUUMMS runtime library includes most of the original C code used in developing VACUUMMS. 
    - vacuumms_cpp.so: The VACUUMMS C++ library contains components used by many of the more recent extensions which are written in C++, and go here.
    - vacuumms_cuda.so: Optional, contains CUDA implementations used by some applications and utilities. Built when BUILD_CUDA_COMPONENTS is set (+cuda in spack)
    - vacuumms_tiff.so: Optional, built when ENABLE_TIFF_UTILS is set (+tiff in spack)

 - Internal applications and Utilities. These are always installed when VACUUMMS is built. They use the VACUUMMS headers and libraries above. 

 - Optional Modules: Currently the variational module is the only module available, though contributions are welcome. It is built when BUILD_VARIATIONAL_MODULE is set (+variational in spack). It uses the above API (headers) and libraries and provides additional headers \<vacuuumms/variational/*\> and provides the additional library vacuumms_variational.so. It also contains developer examples with Makefiles, so that users can see how to build code against the functionality that the module provides, without having to touch the VACUUMMS build system. 

Because of the available API and libraries, VACUUMMS is now no longer solely a top-level code and developers are welcome to use or extend any of the available functionality.

## TO DO: (and PR's are welcome)

- re-work command line options for POSIX style
- re-name uniq to something without a name collision with the Linux utility of the same name, e.g. vniq
- display build options as output when running vacuumms executable. Make the vacuumms executable do more than display version info.

## Developer stuff

**apocrypha** - projects that are not yet complete, have not yet been ported or documented, or which have been replaced by a refactored version

**howto_examples** - sample codes for developers showcasing some of the design patterns employed in VACUUMMS


