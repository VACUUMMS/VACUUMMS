message("Configuring version.h")
message("CMAKE_BINARY_DIR=${CMAKE_BINARYDIR}")
configure_file(version.h.in version.h)
