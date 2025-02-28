
add_executable(cpp_eg "${PROJECT_SOURCE_DIR}/examples/cpp_interface/cpp_eg.cc")
target_include_directories(cav2pov PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(cpp_eg PUBLIC vacuumms_cpp vacuumms_rt m pthread)

install(TARGETS 
        cpp_eg 
        DESTINATION bin
        )

