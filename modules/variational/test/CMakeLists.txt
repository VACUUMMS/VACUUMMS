# CMakeLists.txt for module/variational/test

#option(BUILD_SHARED_LIBS "shared libs " ON)
#option(BUILD_STATIC_LIBS "static libs " OFF)

#add_library(vacuumms_variational SHARED
#      ${PROJECT_SOURCE_DIR}/modules/variational/configuration.cc
##      ${PROJECT_SOURCE_DIR}/modules/variational/constants.cc
#      ${PROJECT_SOURCE_DIR}/modules/variational/quaternion.cc
##      ${PROJECT_SOURCE_DIR}/modules/variational/variational_2D.cc
#      ${PROJECT_SOURCE_DIR}/modules/variational/variational_3D.cc
#)

## What libraries this module depends on
#target_link_libraries(vacuumms_variational PUBLIC ftw_general m)
#
#set_target_properties(vacuumms_variational PROPERTIES
#    POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
#    )
#
#target_include_directories(vacuumms_variational PUBLIC
#    "${PROJECT_SOURCE_DIR}/include"
#    )

add_executable(test_variational_lj "${PROJECT_SOURCE_DIR}/modules/variational/test/test_variational_lj.cc")
target_include_directories(test_variational_lj PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(test_variational_lj PUBLIC vacuumms_variational vacuumms_rt m)

add_executable(test_variational_2D "${PROJECT_SOURCE_DIR}/modules/variational/test/test_variational_2D.cc")
target_include_directories(test_variational_2D PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(test_variational_2D PUBLIC vacuumms_variational vacuumms_rt m)

add_executable(test_variational_3D "${PROJECT_SOURCE_DIR}/modules/variational/test/test_variational_3D.cc")
target_include_directories(test_variational_3D PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(test_variational_3D PUBLIC vacuumms_variational vacuumms_rt m)

add_executable(test_quaternion "${PROJECT_SOURCE_DIR}/modules/variational/test/test_quaternion.cc")
target_include_directories(test_quaternion PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(test_quaternion PUBLIC vacuumms_variational vacuumms_rt m)

add_executable(test_configuration "${PROJECT_SOURCE_DIR}/modules/variational/test/test_configuration.cc")
target_include_directories(test_configuration PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(test_configuration PUBLIC vacuumms_variational vacuumms_rt m)

install(TARGETS test_variational_lj test_variational_2D test_variational_3D test_quaternion test_configuration
        DESTINATION test)

