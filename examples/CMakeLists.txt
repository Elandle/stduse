add_executable(basic)

target_sources(
    basic
    PRIVATE
        basic.f90
)

target_link_libraries(
    basic
    PRIVATE
        stduse)