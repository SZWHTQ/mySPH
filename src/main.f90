program main
    use ctrl_dict, only: maxn, dim, max_interaction
    use parse_toml_m
    use sph
    use initial_m
    use input_m
    use time_integration_m
    implicit none
    integer :: ntotal
    type(Particle), allocatable :: Particles(:)
    integer startT, endT, rate
    character(len=16) :: buffer

    call system_clock(startT)

    call fetch_control_value()

    ! call initialize()
    call allocateParticleList(Particles, maxn, dim, max_interaction/maxn)

    call input(ntotal, Particles)

    call time_integration(ntotal, Particles)

    call system_clock(endT, rate)
    write(buffer, "(F15.3)") dble((endT - startT))/rate
    write(*,*) "Time used: ", trim(adjustl(buffer)), "s"


end program main
