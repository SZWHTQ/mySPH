set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
    APPEND srcs
    "${dir}/controlDict.f90"
    "${dir}/getValueFromToml.f90"
    "${dir}/allocateMemory.f90"
    "${dir}/input.f90"
    "${dir}/SPH.f90"
)

set(srcs "${srcs}" PARENT_SCOPE)