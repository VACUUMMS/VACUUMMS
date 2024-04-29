# test/CMakeLists.txt


#add_executable(expr_add "${PROJECT_SOURCE_DIR}/utils/expr_add.c")
#target_include_directories(expr_add PUBLIC "${PROJECT_SOURCE_DIR}/include")
#target_link_libraries(expr_add PUBLIC vacuumms_rt m ) 




#include(CTest)
#enable_testing()


add_executable(poo ${PROJECT_SOURCE_DIR}/test/poo.cc)
target_include_directories(poo PUBLIC "${PROJECT_SOURCE_DIR}/include")
add_test(test_poo, "${CMAKE_SOURCE_DIR}/cmake/test/poo" 1)

#add_test(ddx_fcc "${CMAKE_SOURCE_DIR}/cmake/applications/pddx < ${PROJECT_SOURCE_DIR}/applications/ddx/fcc.gfg" -box 1 1 1 )
add_test(cat "/usr/bin/cat" "${PROJECT_SOURCE_DIR}/applications/ddx/fcc.gfg"  )
add_test(test_script, "${PROJECT_SOURCE_DIR}/test/script.sh" 1)



#/home/frankwillmore/VACUUMMS/cmake/test/poo

