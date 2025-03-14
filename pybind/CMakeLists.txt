# pybind/CMakeLists.txt

# This creates the 'vacuumms' target based on the bindings decared in the source file here
pybind11_add_module(vacuumms vacuumms_pybind.cc)

# The bindings will need defs of vacuumms types
target_include_directories(vacuumms PUBLIC
    "${PROJECT_SOURCE_DIR}/include"
    )

# vacuumms is the name of the binding module, what gets imported to python
target_link_libraries(vacuumms PUBLIC vacuumms_rt vacuumms_cpp vacuumms_cuda)

install(TARGETS vacuumms
        DESTINATION lib)


