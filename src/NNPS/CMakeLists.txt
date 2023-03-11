set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
    APPEND srcs
    "${dir}/graphics.f90"
    "${dir}/linkList.f90"
    "${dir}/nnps.f90"
    "${dir}/tree.f90"
)

set(srcs "${srcs}" PARENT_SCOPE)