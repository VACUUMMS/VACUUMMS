# test/CMakeLists.txt

# testing the tests:

    # example binary executable test
    add_executable(binary ${PROJECT_SOURCE_DIR}/test/binary.cc)
    target_include_directories(binary PUBLIC "${PROJECT_SOURCE_DIR}/include")
    add_test(test_binary, "${CMAKE_SOURCE_DIR}/cmake/test/binary" 1)

    # example script test
    add_test(test_script, "${PROJECT_SOURCE_DIR}/test/test_script.sh" 1)

    # example compare output test
    add_test(test_compare_output, "${PROJECT_SOURCE_DIR}/test/test_compare_output.sh" 1)

