
option(BUILD_SHARED_LIBS "shared libs " ON)
option(BUILD_STATIC_LIBS "static libs " OFF)

add_library(vacuumms_rt SHARED
      ${PROJECT_SOURCE_DIR}/libraries/vacuumms_rt/prng.c
      ${PROJECT_SOURCE_DIR}/libraries/vacuumms_rt/cav_parser.c
      ${PROJECT_SOURCE_DIR}/libraries/vacuumms_rt/rng2.c
      ${PROJECT_SOURCE_DIR}/libraries/vacuumms_rt/command_line_parser.c  
      ${PROJECT_SOURCE_DIR}/libraries/vacuumms_rt/rng.c
      ${PROJECT_SOURCE_DIR}/libraries/vacuumms_rt/config_parser.c
      ${PROJECT_SOURCE_DIR}/libraries/vacuumms_rt/std.c
      ${PROJECT_SOURCE_DIR}/libraries/vacuumms_rt/param.c
      ${PROJECT_SOURCE_DIR}/libraries/vacuumms_rt/constants.c
      ${PROJECT_SOURCE_DIR}/libraries/vacuumms_rt/povheader.c
)

# This is removed in favor of explicitly adding with target_link_libraries
# As apple-clang and newer gcc builds will fail with this

set_target_properties(vacuumms_rt PROPERTIES
    POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
    )

target_include_directories(vacuumms_rt PUBLIC
    "${PROJECT_SOURCE_DIR}/include"
    )

# vacuumms_cpp

add_library(vacuumms_cpp SHARED
    ${PROJECT_SOURCE_DIR}/libraries/vacuumms_cpp/cavity.cc
    ${PROJECT_SOURCE_DIR}/libraries/vacuumms_cpp/vertex.cc
    ${PROJECT_SOURCE_DIR}/libraries/vacuumms_cpp/edge.cc
    ${PROJECT_SOURCE_DIR}/libraries/vacuumms_cpp/param.cc
    ${PROJECT_SOURCE_DIR}/libraries/vacuumms_cpp/pair.cc
)

target_link_libraries(vacuumms_cpp PUBLIC vacuumms_rt m)

set_target_properties(vacuumms_cpp PROPERTIES
    POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
    )

target_include_directories(vacuumms_cpp PUBLIC
    "${PROJECT_SOURCE_DIR}/include"
)

# optional libs

if (BUILD_TIFF_UTILS)
    add_library(vacuumms_tiff SHARED
        ${PROJECT_SOURCE_DIR}/libraries/vacuumms_tiff/vacuumms_makeTIFF.c
    )
    list(APPEND EXTRA_LIBS vacuumms_tiff)

    target_link_libraries(vacuumms_tiff PUBLIC TIFF::tiff)

    set_target_properties(vacuumms_tiff PROPERTIES
        POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
        )
    target_include_directories(vacuumms_tiff PUBLIC
        "${TIFF_INCLUDE_DIRS}"
        )
    target_include_directories(vacuumms_tiff PUBLIC
        "${PROJECT_SOURCE_DIR}/include"
        )
endif() # BUILD_TIFF_UTILS


if (BUILD_CUDA_COMPONENTS)
    add_library(vacuumms_cuda SHARED
        ${PROJECT_SOURCE_DIR}/libraries/vacuumms_cuda/vacuumms_cuda_utils.cu
        ${PROJECT_SOURCE_DIR}/libraries/vacuumms_cuda/vacuumms_gfg2fvi.cu
    )
    target_include_directories(vacuumms_cuda PUBLIC
       "${PROJECT_SOURCE_DIR}/include"
    )

    target_link_libraries(vacuumms_cuda PUBLIC cuda CUDA::cudart m)

    set_property(TARGET vacuumms_cuda PROPERTY CUDA_ARCHITECTURES OFF)

target_include_directories(vacuumms_cpp PUBLIC
    "${PROJECT_SOURCE_DIR}/include"
)

endif()

install(TARGETS vacuumms_rt vacuumms_cpp
        DESTINATION lib)

if (BUILD_TIFF_UTILS)
    install(TARGETS vacuumms_tiff 
        DESTINATION lib)
endif()

if (BUILD_CUDA_COMPONENTS)
    install(TARGETS vacuumms_cuda
        DESTINATION lib)
endif()

