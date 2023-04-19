#include "../macro.h"
subroutine single_step(ntotal, ndummy, nbuffer, Particles, Delta, aver_v, Shear, dSdt)
    use, intrinsic :: iso_fortran_env, only: err => error_unit
    !$ use omp_lib
    use ctrl_dict,          only: Config, Field
    use time_integration_m, only: Update
    use tools_m,            only: to_string
    use sph,                only: Particle, allocateNeighborList, allocateParticleList
    use nnps_m,             only: search_particles, print_statistics
    use APS_M,              only: BGGS !! Asymmetric Particle Search
    use density_m,          only: sum_density, con_density, con_density_riemann, sum_density_dsph
    use visc_m,             only: viscosity
    use divergence_m,       only: divergence
    use in_force_m,         only: in_force
    use arti_visc_m,        only: arti_visc
    use ex_force_m,         only: ex_force
    use hsml_m,             only: h_upgrade
    use arti_heat_m,        only: arti_heat
    use corr_velo_m,        only: aver_velo
    use dummy_part_m,       only: gen_dummy_particle
    use he_m,               only: detonation_wave
    use decompose_m,        only: decompose
    use area_m,             only: calculate_area
    use bc_m,               only: gen_non_reflecting_bc
    implicit none
    integer,        intent(inout) :: ntotal
    integer,        intent(inout) :: ndummy, nbuffer
    type(Particle), intent(inout) :: Particles(:)
    type(Update),   intent(inout) :: Delta(:)
    real(8),        intent(inout) :: aver_v(:, :)
    integer :: kpair, N
    real(8), dimension(Field%Dim, Field%Maxn) :: indvdt, exdvdt, avdvdt
    real(8), dimension(Field%Maxn) :: indedt, avdedt, ahdedt
    ! real(8) :: area
    ! integer :: number, index(maxn)
#if SOLID
    real(8), intent(in),    optional :: Shear(:, :, :)
    real(8), intent(inout), optional :: dSdt(:, :, :)
#endif
    type(Particle), allocatable :: Buffers(:)

    integer i

#ifdef _OPENMP
    Config%chunkSize = N / Config%nthreads
#endif

    do i = 1, Field%Maxn
        indvdt(:, i) = 0
        indedt(i)    = 0
        exdvdt(:, i) = 0
        avdvdt(:, i) = 0
        avdedt(i)    = 0
        ahdedt(i)    = 0
    end do

    !!! Positions of dummy (boundary) particles
    if ( Config%dummy_parti_w ) call gen_dummy_particle(ndummy, Particles(1:ntotal))

    if ( Config%open_boundary_w) then
        call allocateParticleList(Buffers, nbuffer, Field%Dim, 1)
        call gen_non_reflecting_bc(ntotal, Particles, nbuffer, Buffers)
    else
        nbuffer = 0
    end if

    N = ntotal + ndummy + nbuffer

    !!! Interactions parameters, calculating neighboring particles
    !!! and optimizing smoothing length
    if ( Config%open_boundary_w) then
        call BGGS(Particles(1:ntotal+ndummy), Particles(1:N), skipItsSelf=.true.)
    else
        call search_particles(Config%nnps, Particles(1:N))
    end if

    !!! Density approximation or change rate
    select case (Config%pa_sph)
    case (1, 2)
        if ( Config%sum_density_w ) then
            call sum_density(Particles)
        else
            call con_density(ntotal+ndummy, Particles, Delta%Density)
        end if
    case (3)
        call con_density_riemann(Particles(1:N), Delta%Density)
    case (4)
        call sum_density_dsph(Particles(1:N))
    case default
        write(err, "(1X, A, I0)") "Error density scheme ", Config%pa_sph
        error stop
    end select

    call detonation_wave(Config%i_time_step, Config%delta_t, Particles(1:ntotal))

    call divergence(ntotal+ndummy, Particles)

    !!! Internal forces
#if SOLID
    call in_force(ntotal+ndummy, Particles, indvdt, indedt, Shear, dSdt)
#else
    call in_force(ntotal+ndummy, Particles, indvdt, indedt)
#endif

    !!! Artificial viscosity
    if ( Config%arti_visc_w ) call arti_visc(ntotal+ndummy, Particles, avdvdt, avdedt)

    !!! External force
    if ( Config%ex_force_w ) call ex_force(ntotal+ndummy, Particles, exdvdt)

    if ( Config%arti_heat_w ) call arti_heat(ntotal+ndummy, Particles, ahdedt)

    !!! Calculating average velocity of each particle for avoiding penetration
    if ( Config%aver_velocity_w ) call aver_velo(ntotal, Particles, aver_v)

    !!! Calculating the neighboring particles and updating HSML
    call h_upgrade(ntotal, Particles(1:ntotal))

    !!! Convert velocity, force and energy to f and dfdt
    do i = 1, ntotal
        Delta(i)%Velocity = indvdt(:, i) + avdvdt(:, i) + exdvdt(:, i)
        Delta(i)%Energy   = indedt(i)    + avdedt(i)    + ahdedt(i)
    end do


    if ( mod(Config%i_time_step, Config%print_interval) == 0 ) then
        !!! Statistics for the interaction
        if (Config%print_statistics_w) then
            if ( Config%dummy_parti_w ) then
                write(*,*) ">> Statistics: Dummy Particles:"
                write(*,*) "   Number of dummy particles: ", to_string(ndummy)
            end if
            if ( Config%open_boundary_w ) then
                write(*,*) ">> Statistics: Particles Number:"
                write(*,*) "   Number of particles: ", to_string(ntotal)
            end if
            call print_statistics(Particles(1:ntotal+ndummy))
        end if
            
        ! index = 0
        ! number = 0
        ! do i = 1, ntotal+ndummy
        !     if ( P(i)%Type == 2 ) then
        !         number = number + 1
        !         index(number) = i
        !     end if
        ! end do
        ! call calculate_area(number, mass(index(1:number)), rho(index(1:number)), area)
        ! write(*,"(A, I0, ES15.7)") " >> Statistics: Target Area: ", number, area

        write(*,*) ">> Information for particle ", to_string(Config%monitor_particle)
        write(*,1001) "Internal a ", "Artificial a", "External a", "Total a"
        do i = 1, Field%Dim
            write(*,1002) indvdt(i, Config%monitor_particle), avdvdt(i, Config%monitor_particle), &
                          exdvdt(i, Config%monitor_particle), Delta(Config%monitor_particle)%Velocity(i)
        end do
    end if

    kpair = maxval(Particles%neighborNum)
    if ( kpair > Field%pairNum ) then
        call allocateNeighborList(Particles, Field%Dim, kpair + 3)
    else
        do i = 1, size(Particles)
            Particles(i)%neighborNum = 0
            Particles(i)%neighborList = 0
            Particles(i)%w = 0
            Particles(i)%dwdx = 0
        end do
    end if

    1001 format(A17, A14, 2(A15))
    1002 format(1X, 4(2X, ES13.6))

end subroutine single_step