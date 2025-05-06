cmake_minimum_required(VERSION 3.12)

option(BUILD_SHARED_LIBS "shared libs " ON)
option(BUILD_STATIC_LIBS "static libs " OFF)
option(BUILD_CUDA_COMPONENTS "FVI and PMDGPU require cuda" OFF)
option(BUILD_TIFF_UTILS "some utils require libtiff" OFF)
option(BUILD_VORONOI_UTILS "some utils require voropp" OFF)
option(BUILD_PYBIND_BINDINGS "Build the pybind binding library (vacuumms.so)" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # Optional: Disable compiler-specific extensions

option(ENABLE_TESTING "run the unit tests when building" ON)

# set the project name
if (BUILD_CUDA_COMPONENTS) 
    project(VACUUMMS VERSION 1.2.1 LANGUAGES C CXX CUDA)
else()
    project(VACUUMMS VERSION 1.2.1 LANGUAGES C CXX)
endif()

message("got value for BUILD_CUDA_COMPONENTS=${BUILD_CUDA_COMPONENTS}")
message("got value for BUILD_TIFF_UTILS=${BUILD_TIFF_UTILS}")
message("got value for BUILD_VORONOI_UTILS=${BUILD_VORONOI_UTILS}")
message("got value for ENABLE_TESTING=${ENABLE_TESTING}")
message("got value for BUILD_PYBIND_BINDINGS=${BUILD_PYBIND_BINDINGS}")

#   Find X11. Usually installed to system area, but when it's not, need to find it.
#   It's possible for CMake to use the spack build of X11 with:
#       LIBRARY_PATH=/path/to/spack/build/of/X11/lib cmake ..
#   or by spack load-ing the necessary libs.

find_package(X11 REQUIRED)
message("X11_INCLUDE_DIR set to ${X11_INCLUDE_DIR}")
message("X11_FOUND set to ${X11_FOUND}")
message("X11_X11_LIB set to ${X11_X11_LIB}")
message("X11_SM_FOUND set to ${X11_SM_FOUND}")
message("X11_SM_LIB set to ${X11_SM_LIB}")
message("X11_Xext_LIB set to ${X11_Xext_LIB}")
message("X11_ICE_LIB set to ${X11_ICE_LIB}")
message("X11_LIBRARIES set to ${X11_LIBRARIES}")

if (ENABLE_TESTING)
    message("adding tests subdirectory")
    include(CTest)
#    enable_testing()
    add_subdirectory(test)
endif()

if (BUILD_VORONOI_UTILS)
    add_definitions(-DBUILD_VORONOI_UTILS)
    # There's no find_package for voropp, so just add a path
    # Use the environment variable if it exists
    set(VOROPP_HOME $ENV{VOROPP_HOME})

    if(NOT VOROPP_HOME)
        # Optionally provide a default if the env var isn't set
        # set(MY_LIBRARY_PATH "/usr/lib")
        message("voro++ support requested but voro++ not found.")
        message("Set the VOROPP_HOME environment variable to point to it.")
        message(FATAL_ERROR)
    else()
        message("Setting voropp paths... VOROPP_HOME=${VOROPP_HOME}")
    endif()
endif()
    
if (BUILD_CUDA_COMPONENTS)
    add_definitions(-DBUILD_CUDA_COMPONENTS)
    find_package(CUDA REQUIRED)
    find_package(CUDAToolkit)
endif()

# Import the tiff library. This will make it available for linking. 
# Note that cmake uses FindTIFF.cmake, which defines TIFF::TIFF not TIFF::tiff
if (BUILD_TIFF_UTILS)
    add_definitions(-DBUILD_TIFF_UTILS)
    find_package(TIFF REQUIRED)
#    add_library(tiff SHARED IMPORTED)
#    set_target_properties(tiff PROPERTIES IMPORTED_LOCATION "${TIFF_LIBRARIES}")
endif()

#  the cmake_policy command to set the policy and suppress this warning, or
#  preferably upgrade to using FindPython, either by calling it explicitly
#  before pybind11, or by setting PYBIND11_FINDPYTHON ON before pybind11.

if (BUILD_PYBIND_BINDINGS)
    add_definitions(-DBUILD_PYBIND_BINDINGS)
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
    message(STATUS "Python3_EXECUTABLE: ${Python3_EXECUTABLE}")
    message(STATUS "Python3_LIBRARIES: ${Python3_LIBRARIES}")
    find_package(pybind11 REQUIRED)
    include_directories(${pybind11_INCLUDE_DIRS})
    add_subdirectory(pybind)
endif()

# The standard VACUUMMS stuff

add_subdirectory(include/vacuumms)
add_subdirectory(conda)
add_subdirectory(libraries)
add_subdirectory(utils)
add_subdirectory(applications)

add_subdirectory(examples/cpp_interface)

# This may be erroneous, maybe needs removed
file(COPY conda DESTINATION conda)

# install headers

install(DIRECTORY "${CMAKE_SOURCE_DIR}/" # source directory
        DESTINATION "include" # target directory
        FILES_MATCHING # install only matched files
        PATTERN "*.h" # select C header files
        PATTERN "*.hh" # select C++ header files
)
