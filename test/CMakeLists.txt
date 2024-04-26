# test/CMakeLists.txt


#add_executable(expr_add "${PROJECT_SOURCE_DIR}/utils/expr_add.c")
#target_include_directories(expr_add PUBLIC "${PROJECT_SOURCE_DIR}/include")
#target_link_libraries(expr_add PUBLIC vacuumms_rt m ) 




#include(CTest)
#enable_testing()

add_executable(poo poo.cc)
target_include_directories(poo PUBLIC "${PROJECT_SOURCE_DIR}/include")

add_test(test_poo, "${CMAKE_SOURCE_DIR}/cmake/test/poo" 1)

#/home/frankwillmore/VACUUMMS/cmake/test/poo

