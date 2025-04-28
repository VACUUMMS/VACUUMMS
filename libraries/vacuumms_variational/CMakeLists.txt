# CMakeLists.txt for module/variational

option(BUILD_SHARED_LIBS "shared libs " ON)
option(BUILD_STATIC_LIBS "static libs " OFF)

add_library(vacuumms_variational SHARED
      ${PROJECT_SOURCE_DIR}/libraries/vacuumms_variational/quaternion.cc
      ${PROJECT_SOURCE_DIR}/libraries/vacuumms_variational/variational_2D.cc
      ${PROJECT_SOURCE_DIR}/libraries/vacuumms_variational/variational_3D.cc
)

target_link_libraries(vacuumms_variational 
    PUBLIC vacuumms_rt vacuumms_cpp m
    )

set_target_properties(vacuumms_variational PROPERTIES
    POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
    )

target_include_directories(vacuumms_variational PUBLIC
    "${PROJECT_SOURCE_DIR}/include"
    )

install(TARGETS vacuumms_variational
        DESTINATION lib
    )

add_subdirectory(test)

