add_library(stduse)

add_library(stduse::stduse ALIAS stduse)

target_sources(
    stduse
    PRIVATE
        stduse.f90
        stduse_utilities.f90
        stduse_random.f90
        stduse_print.f90
        stduse_character.f90
)

set_target_properties(
    stduse
    PROPERTIES
        Fortran_MODULE_DIRECTORY
            "${CMAKE_CURRENT_BINARY_DIR}/modules"
)

target_include_directories(
    stduse
    PUBLIC
        "${CMAKE_CURRENT_BINARY_DIR}/modules"
)