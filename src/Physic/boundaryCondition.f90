module boundary_condition_m
    use sph
    use geometry_m
    implicit none

    private

contains
    !!! Open Boundary Non-Reflecting Condition
    subroutine none_reflecting_bc(Particles, FluidDomain, BufferDomain)
        class(geometry_t), intent(in) :: FluidDomain(:), BufferDomain(:)
        type(Particle), intent(inout) :: Particles(:)
        type(Particle), allocatable :: GhostParticles(:)
    
        integer i
        
    end subroutine none_reflecting_bc

end module boundary_condition_m