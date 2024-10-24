cmake_minimum_required(VERSION 3.12)

# set the project name
project(VACUUMMS LANGUAGES C CXX CUDA)

option(BUILD_SHARED_LIBS "shared libs " ON)
option(BUILD_STATIC_LIBS "static libs " OFF)
option(BUILD_CUDA_COMPONENTS "FVI and PMDGPU require cuda" OFF)
option(BUILD_TIFF_UTILS "some utils require libtiff" OFF)
option(BUILD_VORONOI_UTILS "some utils require voropp" OFF)
option(BUILD_VARIATIONAL_MODULE "variational calculation of energy paths" OFF)
option(ENABLE_TESTING "run the unit tests when building" ON)


message("got value for BUILD_CUDA_COMPONENTS=${BUILD_CUDA_COMPONENTS}")
message("got value for BUILD_TIFF_UTILS=${BUILD_TIFF_UTILS}")
message("got value for BUILD_VARIATIONAL_MODULE=${BUILD_VARIATIONAL_MODULE}")
message("got value for BUILD_VORONOI_UTILS=${BUILD_VORONOI_UTILS}")
message("got value for ENABLE_TESTING=${ENABLE_TESTING}")

# Find X11. Usually installed to system area, but when it's not, need to find it.
find_package(X11)
message("X11_INCLUDE_DIR set to ${X11_INCLUDE_DIR}")
message("X11_FOUND set to ${X11_FOUND}")
message("X11_X11_LIB set to ${X11_X11_LIB}")
message("X11_SM_FOUND set to ${X11_SM_FOUND}")
message("X11_SM_LIB set to ${X11_SM_LIB}")
message("X11_Xext_LIB set to ${X11_Xext_LIB}")
message("X11_ICE_LIB set to ${X11_ICE_LIB}")
message("X11_LIBRARIES set to ${X11_LIBRARIES}")

# FTW: Below is an attempt to make the package buildable from CMake with a non-system install
#      of X11, but this effort has been abandoned for now. CMake standalone build will work
#      on a system with an X11 install to /usr, and a spack build will work with non-system X11.
#      10/22/2024 note: I was able to get CMake to use the spack build with:
#      LIBRARY_PATH=/path/to/spack/build/of/X11/lib cmake ..

#list(PREPEND X11_LIBRARIES ${X11_ICE_LIB})
#message("X11_LIBRARIES now set to ${X11_LIBRARIES}")
#list(PREPEND X11_INCLUDE_DIR ${X11_ICE_INCLUDE_PATH})
#message("X11_INCLUDE_DIR now set to ${X11_INCLUDE_DIR}")
#add_library(X11 SHARED IMPORTED)
#set_target_properties(X11 PROPERTIES IMPORTED_LOCATION "${X11_LIBRARIES}")
#add_library(SM SHARED IMPORTED)
#set_target_properties(SM PROPERTIES IMPORTED_LOCATION "${X11_SM_LIB}")
#add_library(ICE SHARED IMPORTED)
#set_target_properties(ICE PROPERTIES IMPORTED_LOCATION "${X11_ICE_INCLUDE_PATH}")
#add_library(Xext SHARED IMPORTED)
#set_target_properties(Xext PROPERTIES IMPORTED_LOCATION "${X11_Xext_LIB}")

# enable/disable testing

if (ENABLE_TESTING)
    message("adding tests subdirectory")
    include(CTest)
    enable_testing()
    add_subdirectory(test)
endif()

if (BUILD_CUDA_COMPONENTS)
    find_package(CUDA)
    find_package(CUDAToolkit)
endif()

# Import the tiff library. This will make it available for linking. 
if (BUILD_TIFF_UTILS)
    find_package(TIFF)
    add_library(tiff SHARED IMPORTED)
    set_target_properties(tiff PROPERTIES IMPORTED_LOCATION "${TIFF_LIBRARIES}")
endif()

# The standard VACUUMMS stuff

add_subdirectory(libraries)
add_subdirectory(utils)
add_subdirectory(applications)

# There's no find_package for voropp, so just add a path
if (BUILD_VORONOI_UTILS)
# how do we set VOROPP_HOME? and check it is set?
    # set paths here
    add_library(voro++ SHARED IMPORTED)
    message("Setting voropp paths... VOROPP_HOME=${VOROPP_HOME}")
    #    set_target_properties(voro++ PROPERTIES IMPORTED_LOCATION "${VOROPP_HOME}/lib")
    set_target_properties(voro++ PROPERTIES IMPORTED_LOCATION "${VOROPP_HOME}/lib/libvoro++.a")
endif()

# optional modules

if (BUILD_VARIATIONAL_MODULE)
    message("adding variational module")
    list(APPEND EXTRA_lIBS vacuumms_variational)
    add_subdirectory(modules/variational)
endif()

# install headers

install(DIRECTORY "${CMAKE_SOURCE_DIR}/" # source directory
        DESTINATION "include" # target directory
        FILES_MATCHING # install only matched files
        PATTERN "*.h" # select C header files
        PATTERN "*.hh" # select C++ header files
)

