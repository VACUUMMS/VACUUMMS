cmake_minimum_required(VERSION 3.10)

option(BUILD_SHARED_LIBS "shared libs " ON)
option(BUILD_STATIC_LIBS "static libs " OFF)

## Now the executables...

add_executable(rattle "${PROJECT_SOURCE_DIR}/applications/rattle/rattle.c")
target_include_directories(rattle PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(rattle PUBLIC vacuumms_rt m pthread)

add_executable(ljx 
    "${PROJECT_SOURCE_DIR}/applications/ljx/ljx_main.c"
    "${PROJECT_SOURCE_DIR}/applications/ljx/command_line_parser.c"
    "${PROJECT_SOURCE_DIR}/applications/ljx/io_setup.c"
    "${PROJECT_SOURCE_DIR}/applications/ljx/energy.c"
    "${PROJECT_SOURCE_DIR}/applications/ljx/graphics.c"
)
target_include_directories(ljx PUBLIC 
    "${PROJECT_SOURCE_DIR}/include"
    "${X11_INCLUDE_DIR}"
)
target_link_libraries(ljx PUBLIC vacuumms_rt m pthread X11)

add_executable(hs 
    "${PROJECT_SOURCE_DIR}/applications/hs/hs_main.c"
    "${PROJECT_SOURCE_DIR}/applications/hs/io_setup.c"
    "${PROJECT_SOURCE_DIR}/applications/hs/graphics.c"
)
target_include_directories(hs PUBLIC 
    "${PROJECT_SOURCE_DIR}/include"
    "${X11_INCLUDE_DIR}"
)
target_link_libraries(hs PUBLIC vacuumms_rt m pthread X11)

add_executable(center "${PROJECT_SOURCE_DIR}/applications/center/center.c")
target_include_directories(center PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(center PUBLIC vacuumms_rt m pthread)

add_executable(mfp
    "${PROJECT_SOURCE_DIR}/applications/mfp/mfp.c"
    "${PROJECT_SOURCE_DIR}/applications/mfp/io_setup.c"
)
target_include_directories(mfp PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(mfp PUBLIC vacuumms_rt m pthread)

add_executable(ddx "${PROJECT_SOURCE_DIR}/applications/ddx/ddx.c" "${PROJECT_SOURCE_DIR}/applications/ddx/io_setup.c" )
target_include_directories(ddx PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(ddx PUBLIC vacuumms_rt m)

add_executable(pddx "${PROJECT_SOURCE_DIR}/applications/ddx/pddx.c" "${PROJECT_SOURCE_DIR}/applications/ddx/io_setup.c" )
#add_executable(pddx "${PROJECT_SOURCE_DIR}/applications/ddx/pddx.c" )
target_include_directories(pddx PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(pddx PUBLIC vacuumms_rt m pthread)

add_executable(vddx "${PROJECT_SOURCE_DIR}/applications/ddx/vddx.c" 
                    "${PROJECT_SOURCE_DIR}/applications/ddx/vddx_io_setup.c" 
                    "${PROJECT_SOURCE_DIR}/applications/ddx/io_setup.c" )
target_include_directories(vddx PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(vddx PUBLIC vacuumms_rt m pthread)

add_executable(vuniq "${PROJECT_SOURCE_DIR}/applications/ddx/vuniq.cc") 
target_include_directories(vuniq PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(vuniq PUBLIC vacuumms_rt m pthread vacuumms_cpp)

add_executable(edges2var "${PROJECT_SOURCE_DIR}/applications/ddx/edges2var.cc")
target_include_directories(edges2var PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(edges2var PUBLIC vacuumms_rt vacuumms_cpp m)

add_executable(2pc 
    "${PROJECT_SOURCE_DIR}/applications/2pc/2pc.c" 
)
target_include_directories(2pc PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(2pc PUBLIC vacuumms_rt m )

add_executable(csd 
    "${PROJECT_SOURCE_DIR}/applications/csd/csd.c" 
    "${PROJECT_SOURCE_DIR}/applications/csd/command_line_parser.c" 
)
target_include_directories(csd PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(csd PUBLIC vacuumms_rt m pthread)

add_executable(rog 
    "${PROJECT_SOURCE_DIR}/applications/rog/rog.c" 
)
target_include_directories(rog PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(rog PUBLIC vacuumms_rt m pthread)

add_executable(cv 
    "${PROJECT_SOURCE_DIR}/applications/sv/cv.c" 
)
target_include_directories(cv PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(cv PUBLIC vacuumms_rt m )

add_executable(csa 
    "${PROJECT_SOURCE_DIR}/applications/sv/csa.c" 
)
target_include_directories(csa PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(csa PUBLIC vacuumms_rt m )

add_executable(cav2cluster 
    "${PROJECT_SOURCE_DIR}/applications/cav2cluster/cav2cluster.c" 
)
target_include_directories(cav2cluster PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(cav2cluster PUBLIC vacuumms_rt m )

add_executable(sdhs 
    "${PROJECT_SOURCE_DIR}/applications/sdhs/sdhs_main.c" 
    "${PROJECT_SOURCE_DIR}/applications/sdhs/io_setup.c" 
    "${PROJECT_SOURCE_DIR}/applications/sdhs/graphics.c" 
)
target_include_directories(sdhs PUBLIC 
    "${PROJECT_SOURCE_DIR}/include"
    "${X11_INCLUDE_DIR}"
)
target_link_libraries(sdhs PUBLIC vacuumms_rt m X11 )

add_executable(end2end 
    "${PROJECT_SOURCE_DIR}/applications/end2end/end2end.c" 
)
target_include_directories(end2end PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(end2end PUBLIC vacuumms_rt m pthread)

add_executable(vis 
    "${PROJECT_SOURCE_DIR}/applications/vis/vis.c" 
    "${PROJECT_SOURCE_DIR}/applications/vis/command_line_parser.c"
    "${PROJECT_SOURCE_DIR}/applications/vis/io_setup.c" 
    "${PROJECT_SOURCE_DIR}/applications/vis/graphics.c" 
)
target_include_directories(vis PUBLIC 
    "${PROJECT_SOURCE_DIR}/include"
    "${X11_INCLUDE_DIR}"
)
target_link_libraries(vis PUBLIC vacuumms_rt m pthread X11)

add_executable(tcv 
    "${PROJECT_SOURCE_DIR}/applications/sv/tcv.cc" 
)
target_include_directories(tcv PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(tcv PUBLIC vacuumms_cpp vacuumms_rt m )

install(TARGETS 
        rattle 
        ljx 
        center 
        mfp 
        ddx 
        pddx 
        vddx 
        vuniq
        edges2var
        hs 
        2pc
        csd    
        rog  
        cav2cluster
        sdhs
        end2end
        cv
        csa
        vis
        DESTINATION bin)

# testing

add_test(
    NAME test_ddx_fcc
    COMMAND "${PROJECT_SOURCE_DIR}/applications/test/test_ddx_fcc.sh"
)
set_property(
    TEST test_ddx_fcc
    PROPERTY ENVIRONMENT BinaryPath=${CMAKE_SOURCE_DIR} TestDataPath=${PROJECT_SOURCE_DIR}
)

add_test(
    NAME test_pddx_fcc
    COMMAND "${PROJECT_SOURCE_DIR}/applications/test/test_pddx_fcc.sh"
)
set_property(
    TEST test_pddx_fcc
    PROPERTY ENVIRONMENT BinaryPath=${CMAKE_SOURCE_DIR} TestDataPath=${PROJECT_SOURCE_DIR}
)

add_test(
    NAME test_ljx
    COMMAND "${PROJECT_SOURCE_DIR}/applications/test/test_ljx.sh"
)
set_property(
    TEST test_ljx
    PROPERTY ENVIRONMENT BinaryPath=${CMAKE_SOURCE_DIR} TestDataPath=${PROJECT_SOURCE_DIR}
)

