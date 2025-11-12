# test/CMakeLists.txt

enable_testing()
include(CTest)

# Define global env vars in a variable, if we don't want to type for every test
set(MY_TEST_ENV
    "BinaryPath=${CMAKE_BINARY_DIR}"
    "SourcePath=${CMAKE_SOURCE_DIR}"
    "TestDataPath=${CMAKE_SOURCE_DIR}/test/data"
    "PYTHONPATH=${CMAKE_BINARY_DIR}/pybind:$ENV{PYTHONPATH}"
    "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/libraries:${CMAKE_BINARY_DIR}/libraries/vacuumms_variational"
)

# Testing the tester section.

# example binary executable test
add_executable(binary ${CMAKE_SOURCE_DIR}/test/src/binary.cc)
add_test(test_binary "${CMAKE_BINARY_DIR}/test/binary" 1)

# example with assertion
add_executable(assert ${CMAKE_SOURCE_DIR}/test/src/assert.cc)
add_test(test_assert "${CMAKE_BINARY_DIR}/test/assert" 0)

# example script test
add_test(test_script "${CMAKE_SOURCE_DIR}/test/scripts/test_script.sh" 1)

# example compare output test
add_test(
    NAME test_compare_output 
    COMMAND "${CMAKE_SOURCE_DIR}/test/scripts/test_compare_output.sh" 1
#    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test
)


# PyBind tests:

# Ensure the module is in PYTHONPATH and loadable
add_test(
    NAME test_pybind_import
    COMMAND python3 ${PROJECT_SOURCE_DIR}/test/python/test_import.py 
)

add_test(
    NAME test_lammps_import
    COMMAND ${PROJECT_SOURCE_DIR}/test/scripts/test_lammps_import.sh 
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test
)

add_test(
    NAME test_ddx_fcc_py
    COMMAND ${PROJECT_SOURCE_DIR}/test/scripts/test_ddx_fcc_py.sh
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test
)

#add_test(
#    NAME test_pddx_fcc_py
#    COMMAND ${PROJECT_SOURCE_DIR}/test/scripts/test_pddx_fcc_py.sh
#    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test
#)

add_test(
    NAME test_fvi_fcc_py
    COMMAND ${PROJECT_SOURCE_DIR}/test/scripts/test_fvi_fcc_py.sh
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test
)


# utils testing


add_test(
    NAME test_uniq
    COMMAND "${PROJECT_SOURCE_DIR}/test/scripts/test_uniq.sh"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)


# applications testing

add_test(
    NAME test_ddx_fcc
    COMMAND "${PROJECT_SOURCE_DIR}/test/scripts/test_ddx_fcc.sh"
)

add_test(
    NAME test_pddx_fcc
    COMMAND "${PROJECT_SOURCE_DIR}/test/scripts/test_pddx_fcc.sh"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test
)

add_test(
    NAME test_ljx
    COMMAND "${PROJECT_SOURCE_DIR}/test/scripts/test_ljx.sh"
)


# Variational module tests

add_test(
    NAME test_variational_lj
    COMMAND "${PROJECT_SOURCE_DIR}/test/scripts/test_variational_lj.sh"
)

add_test(
    NAME test_variational_2D
    COMMAND "${PROJECT_SOURCE_DIR}/test/scripts/test_variational_2D.sh"
)

add_test(
    NAME test_variational_3D
    COMMAND "${PROJECT_SOURCE_DIR}/test/scripts/test_variational_3D.sh"
)

add_test(
    NAME test_quaternion
    COMMAND "${PROJECT_SOURCE_DIR}/test/scripts/test_quaternion.sh"
)

#get_test_list(all_tests)
get_property(test_names DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY TESTS)
foreach(test_name ${test_names})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "${MY_TEST_ENV}")
endforeach()

