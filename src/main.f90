program main
    use ctrl_dict,          only: Field
    use parse_toml_m,       only: fetch_control_value
    use sph,                only: Particle, allocateParticleList
    use input_m,            only: input
    use time_integration_m, only: time_integration
    use tools_m,            only: to_string, round
    implicit none
    integer :: ntotal
    type(Particle), allocatable :: Particles(:)
    integer startT, endT, rate
    character(len=16) :: buffer

    call system_clock(startT)

    call fetch_control_value()

    call allocateParticleList(Particles, Field%Maxn, Field%Dim, Field%pairNum)

    call input(ntotal, Particles)

    call time_integration(ntotal, Particles)

    call system_clock(endT, rate)
    write(buffer, "(F15.3)") dble((endT - startT))/rate
    write(*,*) "Time used: ", trim(adjustl(buffer)), "s"


end program main
