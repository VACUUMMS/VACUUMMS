# test/CMakeLists.txt

# example binary executable test
add_executable(binary ${CMAKE_SOURCE_DIR}/test/src/binary.cc)
#target_include_directories(binary PUBLIC "${PROJECT_SOURCE_DIR}/include")
add_test(test_binary, "${CMAKE_BINARY_DIR}/test/binary" 1)

add_executable(assert ${CMAKE_SOURCE_DIR}/test/src/assert.cc)
add_test(test_assert_0, "${CMAKE_BINARY_DIR}/test/assert" 0)

# example script test
add_test(test_script, "${CMAKE_SOURCE_DIR}/test/scripts/test_script.sh" 1)

# example compare output test
add_test(NAME test_compare_output, 
    COMMAND "${CMAKE_SOURCE_DIR}/test/scripts/test_compare_output.sh" 1
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test)

# Ensure the module is in PYTHONPATH
add_test(
    NAME test_pybind_import
    COMMAND python3 ${PROJECT_SOURCE_DIR}/test/python/test_import.py 
)

set_tests_properties(test_pybind_import 
    PROPERTIES
    ENVIRONMENT "PYTHONPATH=/home1/fwillmor/opt/vacuumms/lib:$ENV{PYTHONPATH}"
)
