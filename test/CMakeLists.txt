# test/CMakeLists.txt

# Define global env vars in a variable, if we don't want to type for every test
set(MY_TEST_ENV
    "BinaryPath=${CMAKE_BINARY_DIR}"
    "TestDataPath=${CMAKE_SOURCE_DIR}/test_data"
    "PYTHONPATH=${CMAKE_BINARY_DIR}:$ENV{PYTHONPATH}"
)

enable_testing()

# Testing the tester section.

# example binary executable test
add_executable(binary ${CMAKE_SOURCE_DIR}/test/src/binary.cc)
add_test(test_binary, "${CMAKE_BINARY_DIR}/test/binary" 1)

# example with assertion
add_executable(assert ${CMAKE_SOURCE_DIR}/test/src/assert.cc)
add_test(test_assert, "${CMAKE_BINARY_DIR}/test/assert" 0)

# example script test
add_test(test_script, "${CMAKE_SOURCE_DIR}/test/scripts/test_script.sh" 1)

# example compare output test
add_test(
    NAME test_compare_output 
    COMMAND "${CMAKE_SOURCE_DIR}/test/scripts/test_compare_output.sh" 1
#    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test
)
set_tests_properties(test_compare_output PROPERTIES
    ENVIRONMENT "TestDataPath=${CMAKE_SOURCE_DIR}/test/data"
)



# PyBind tests:

# Ensure the module is in PYTHONPATH and loadable
add_test(
    NAME test_pybind_import
    COMMAND python3 ${PROJECT_SOURCE_DIR}/test/python/test_import.py 
)
set_tests_properties(test_pybind_import 
    PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/pybind:$ENV{PYTHONPATH};"
)

add_test(
    NAME test_lammps_import
    COMMAND ${PROJECT_SOURCE_DIR}/test/scripts/test_lammps_import.sh 
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test
)
set_tests_properties(test_lammps_import 
    PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/pybind:$ENV{PYTHONPATH};SourcePath=${CMAKE_SOURCE_DIR};TestDataPath=${CMAKE_SOURCE_DIR}/test/data"
)

add_test(
    NAME test_ddx_fcc_py
    COMMAND ${PROJECT_SOURCE_DIR}/test/scripts/test_ddx_fcc_py.sh
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test
)
set_tests_properties(test_ddx_fcc_py
    PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/pybind:$ENV{PYTHONPATH};SourcePath=${CMAKE_SOURCE_DIR};TestDataPath=${CMAKE_SOURCE_DIR}/test/data"
)

add_test(
    NAME test_pddx_fcc_py
    COMMAND ${PROJECT_SOURCE_DIR}/test/scripts/test_pddx_fcc_py.sh
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test
)
set_tests_properties(test_pddx_fcc_py
    PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/pybind:$ENV{PYTHONPATH};SourcePath=${CMAKE_SOURCE_DIR};TestDataPath=${CMAKE_SOURCE_DIR}/test/data"
)

add_test(
    NAME test_fvi_fcc_py
    COMMAND ${PROJECT_SOURCE_DIR}/test/scripts/test_fvi_fcc_py.sh
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test
)
set_tests_properties(test_fvi_fcc_py
    PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/pybind:$ENV{PYTHONPATH};SourcePath=${CMAKE_SOURCE_DIR};TestDataPath=${CMAKE_SOURCE_DIR}/test/data"
)




# utils testing


add_test(
    NAME test_uniq
    COMMAND "${PROJECT_SOURCE_DIR}/test/scripts/test_uniq.sh"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(test_uniq PROPERTIES
    ENVIRONMENT "BinaryPath=${CMAKE_BINARY_DIR};SourcePath=${CMAKE_SOURCE_DIR};PYTHONPATH=${CMAKE_BINARY_DIR}:$ENV{PYTHONPATH}"
)


# applications testing

add_test(
    NAME test_ddx_fcc
    COMMAND "${PROJECT_SOURCE_DIR}/test/scripts/test_ddx_fcc.sh"
)
set_tests_properties(test_ddx_fcc PROPERTIES
    ENVIRONMENT "BinaryPath=${CMAKE_BINARY_DIR};TestDataPath=${PROJECT_SOURCE_DIR}/test/data"
)

add_test(
    NAME test_pddx_fcc
    COMMAND "${PROJECT_SOURCE_DIR}/test/scripts/test_pddx_fcc.sh"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test
)
set_tests_properties(test_pddx_fcc PROPERTIES
    ENVIRONMENT "BinaryPath=${CMAKE_BINARY_DIR};TestDataPath=${PROJECT_SOURCE_DIR}/test/data"
)

add_test(
    NAME test_ljx
    COMMAND "${PROJECT_SOURCE_DIR}/test/scripts/test_ljx.sh"
)
set_tests_properties(test_ljx PROPERTIES
    ENVIRONMENT "BinaryPath=${CMAKE_BINARY_DIR};TestDataPath=${PROJECT_SOURCE_DIR}/test/data"
)


