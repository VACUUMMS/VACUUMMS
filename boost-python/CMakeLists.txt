# CMakeLists.txt for boost::python connector lib

set(CMAKE_SHARED_LIBRARY_PREFIX "")

#find_package(Boost REQUIRED COMPONENTS system filesystem)
#find_package(Python COMPONENTS Interpreter Development)

add_library(vacuumms SHARED
      ${PROJECT_SOURCE_DIR}/boost-python/bp_def.cc
)

set_target_properties(vacuumms PROPERTIES
    POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
    )

target_include_directories(vacuumms PUBLIC
    "${PROJECT_SOURCE_DIR}/include"
    )

target_link_libraries(vacuumms PUBLIC 
    vacuumms_rt 
    vacuumms_cpp 
    boost_python312 
    python3.12
    m)

install(TARGETS vacuumms
        DESTINATION lib)

