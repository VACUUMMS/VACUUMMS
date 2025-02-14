#cmake_minimum_required(VERSION 3.14)

#project(cmake_example)

add_library(binding_lib SHARED lib.h lib.cpp)

#find_package(pybind11 REQUIRED)

# add /apps/spack/2406/apps/linux-rocky8-x86_64_v3/gcc-13.3.0/python-3.11.9-x74mtjf/include/python3.11/
# to include path

#pybind11_add_module(binding_module module.cpp)
#pybind11_add_module(vacuumms module.cpp)
pybind11_add_module(vacuumms v.cpp)


#target_link_libraries(binding_module PUBLIC binding_lib)
#target_link_libraries(vacuumms PUBLIC binding_lib)
target_link_libraries(vacuumms PUBLIC binding_lib vacuumms_rt vacuumms_cpp)

target_include_directories(vacuumms PUBLIC
    "${PROJECT_SOURCE_DIR}/include"
    )

install(TARGETS vacuumms binding_lib
        DESTINATION lib)


