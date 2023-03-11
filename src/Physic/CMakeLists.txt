set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
    APPEND srcs
    "${dir}/boundaryCondition.f90"
    "${dir}/density.f90"
    "${dir}/eos.f90"
    "${dir}/externalForce.f90"
    "${dir}/HE.f90"
    "${dir}/internalForce.f90"
    "${dir}/viscosity.f90"
)

set(srcs "${srcs}" PARENT_SCOPE)