# CMakeLists.txt for module/variational

#cmake_minimum_required(VERSION 3.10)

option(BUILD_SHARED_LIBS "shared libs " ON)
option(BUILD_STATIC_LIBS "static libs " OFF)

add_library(vacuumms_variational SHARED
      ${PROJECT_SOURCE_DIR}/modules/variational/configuration.cc
#      ${PROJECT_SOURCE_DIR}/modules/variational/constants.cc
      ${PROJECT_SOURCE_DIR}/modules/variational/quaternion.cc
      ${PROJECT_SOURCE_DIR}/modules/variational/variational_2D.cc
      ${PROJECT_SOURCE_DIR}/modules/variational/variational_3D.cc
)

# What libraries this module depends on
target_link_libraries(vacuumms_variational PUBLIC ftw_general m)

set_target_properties(vacuumms_variational PROPERTIES
    POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
    )

target_include_directories(vacuumms_variational PUBLIC
    "${PROJECT_SOURCE_DIR}/include"
    )

#add_executable(test_variational_lj "${PROJECT_SOURCE_DIR}/modules/variational/test/test_variational_lj.cc")
#target_include_directories(test_variational_lj PUBLIC "${PROJECT_SOURCE_DIR}/include")
#target_link_libraries(test_variational_lj PUBLIC vacuumms_variational ftw_general m)

install(TARGETS vacuumms_variational
        DESTINATION lib)

add_subdirectory(test)

