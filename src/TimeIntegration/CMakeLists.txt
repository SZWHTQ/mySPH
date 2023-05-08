set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
    APPEND srcs
    "${dir}/singleStep.f90"
    "${dir}/timeIntegration.f90"
)

set(srcs "${srcs}" PARENT_SCOPE)