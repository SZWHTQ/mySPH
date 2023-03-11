set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
    APPEND srcs
    "${dir}/tools.f90"
)

add_subdirectory("Initialize")
add_subdirectory("NNPS")
add_subdirectory("Numeric")
add_subdirectory("Output")
add_subdirectory("Physic")
add_subdirectory("TimeIntegration")

list(
    APPEND srcs
    # "${dir}/Initialize/input.f90"
    # "${dir}/NNPS/nnps.f90"
    "${dir}/main.f90"
)

set(srcs "${srcs}" PARENT_SCOPE)