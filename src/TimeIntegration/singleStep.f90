#include "../macro.h"
subroutine single_step(ntotal, ndummy, Particles, Delta, aver_v, Shear, dSdt)
    use, intrinsic :: iso_fortran_env, only: err => error_unit
    !$ use omp_lib
    use ctrl_dict,          only: Config, Field
    use time_integration_m, only: Update
    use tools_m,            only: to_string
    use sph,                only: Particle, allocateNeighborList
    use nnps_m,             only: search_particles, print_statistics
    use BGGS_m,             only: BGGS
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
    use bc_m,               only: gen_non_reflecting_bc, calc_nrbc_property
    implicit none
    integer,        intent(inout) :: ntotal
    integer,        intent(inout) :: ndummy
    type(Particle), intent(inout) :: Particles(:)
    type(Update),   intent(inout) :: Delta(:)
    real(8),        intent(inout) :: aver_v(:, :)
    real(8), dimension(Field%dim, Field%maxn) :: indvdt, exdvdt, avdvdt
    real(8), dimension(Field%maxn) :: indedt, avdedt, ahdedt
    ! real(8) :: area
    ! integer :: number, index(maxn)
#if SOLID
    real(8), intent(in),    optional :: Shear(:, :, :)
    real(8), intent(inout), optional :: dSdt(:, :, :)
#endif

    integer i

    do i = 1, Field%maxn
        indvdt(:, i) = 0
        indedt(i)    = 0
        exdvdt(:, i) = 0
        avdvdt(:, i) = 0
        avdedt(i)    = 0
        ahdedt(i)    = 0
    end do

    !!! Positions of dummy (boundary) particles
    if ( Config%dummy_parti_w ) call gen_dummy_particle(ndummy, Particles(1:ntotal))

    if ( Config%open_boundary_w) call gen_non_reflecting_bc(ntotal, Particles)

#ifdef _OPENMP
    Config%chunkSize = (ntotal+ndummy) / Config%nthreads
#endif
    !!! Interactions parameters, calculating neighboring particles
    !!! and optimizing smoothing length
    call search_particles(Config%nnps, Particles(1:ntotal+ndummy))

    if ( Config%open_boundary_w) call calc_nrbc_property(Particles(1:ntotal))

    !!! Density approximation or change rate
    select case (Config%pa_sph)
    case (1, 2)
        if ( Config%sum_density_w ) then
            call sum_density(Particles(1:ntotal+ndummy))
        else
            call con_density(Particles(1:ntotal+ndummy), Delta%Density)
        end if
    case (3)
        call con_density_riemann(Particles(1:ntotal+ndummy), Delta%Density)
    case (4)
        call sum_density_dsph(Particles(1:ntotal+ndummy))
    case default
        write(err, "(1X, A, I0)") "Error density scheme ", Config%pa_sph
        error stop
    end select

    call detonation_wave(Config%i_time_step, Config%delta_t, Particles(1:ntotal))

    call divergence(Particles(1:ntotal+ndummy))

    !!! Internal forces
#if SOLID
    call in_force(Particles(1:ntotal+ndummy), indvdt, indedt, Shear, dSdt)
#else
    call in_force(Particles(1:ntotal+ndummy), indvdt, indedt)
#endif

    !!! Artificial viscosity
    if ( Config%arti_visc_w ) call arti_visc(Particles(1:ntotal+ndummy), avdvdt, avdedt)

    !!! External force
    if ( Config%ex_force_w ) call ex_force(Particles(1:ntotal+ndummy), exdvdt)

    if ( Config%arti_heat_w ) call arti_heat(Particles(1:ntotal+ndummy), ahdedt)

    !!! Calculating average velocity of each particle for avoiding penetration
    if ( Config%aver_velocity_w ) call aver_velo(Particles(1:ntotal), aver_v)

    !!! Calculating the neighboring particles and updating HSML
    call h_upgrade(Config%sle, Config%delta_t, Particles(1:ntotal))

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
        do i = 1, Field%dim
            write(*,1002) indvdt(i, Config%monitor_particle), avdvdt(i, Config%monitor_particle), &
                          exdvdt(i, Config%monitor_particle), Delta(Config%monitor_particle)%Velocity(i)
        end do
    end if

    call allocateNeighborList(Particles, Field%dim, maxval(Particles%neighborNum)+5)

    1001 format(A17, A14, 2(A15))
    1002 format(1X, 4(2X, ES13.6))

end subroutine single_step