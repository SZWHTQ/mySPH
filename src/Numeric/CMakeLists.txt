set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
    APPEND srcs
    "${dir}/area.f90"
    "${dir}/artificialHeat.f90"
    "${dir}/artificialTerm.f90"
    "${dir}/artificialViscosity.f90"
    "${dir}/correctVelocity.f90"
    "${dir}/courantNumber.f90"
    "${dir}/decomposePar.f90"
    "${dir}/divergence.f90"
    "${dir}/hsml.f90"
    "${dir}/kernel.f90"
    "${dir}/kernelGradientCorrection.f90"
    "${dir}/shifting.f90"
)

set(srcs "${srcs}" PARENT_SCOPE)