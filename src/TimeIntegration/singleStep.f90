subroutine single_step(ntotal, ndummy, nbuffer, Particles, Delta, aver_v, Shear, dSdt)
#include "../macro.h"
    use, intrinsic :: iso_fortran_env, only: err => error_unit
    !$ use omp_lib
    use ctrl_dict,          only: Config, Field
    use tools_m,            only: to_string
    use SPH,                only: Particle, Update, allocateNeighborList, allocateParticleList
    use nnps_m,             only: search_particles, print_statistics
    use APS_M,              only: BGGS !! Asymmetric Particle Search
    use KGC_m,              only: kernelGradientCorrection
    use shifting_m,         only: shifting
    use density_m,          only: sum_density, con_density, con_density_riemann, &
                                  sum_density_dsph, norm_density
    use visc_m,             only: viscosity
    use divergence_m,       only: divergence
    use in_force_m,         only: internal_force
    use arti_visc_m,        only: arti_visc
    use arti_heat_m,        only: arti_heat
    use AT_m,               only: artificialStress
    use ex_force_m,         only: external_force
    use hsml_m,             only: h_upgrade
    use corr_velo_m,        only: aver_velo
    use dummy_part_m,       only: gen_dummy_particle
    use he_m,               only: detonation_wave
    use decompose_m,        only: decompose
    use area_m,             only: calculate_area
    use bc_m,               only: gen_non_reflecting_bc, boundary
    implicit none
    integer,        intent(inout) :: ntotal
    integer,        intent(inout) :: ndummy, nbuffer
    type(Particle), intent(inout) :: Particles(:)
    type(Update),   intent(inout) :: Delta(:)
    real(8),        intent(inout) :: aver_v(:, :)
    integer :: N!, kpair
    real(8), dimension(Field%Dim, Field%Maxn) :: indvdt, exdvdt, avdvdt, asdvdt
    real(8), dimension(Field%Maxn) :: indedt, avdedt, ahdedt
    ! real(8) :: area
    ! integer :: number, index(maxn)
#if SOLID
    real(8), intent(in),    optional :: Shear(:, :, :)
    real(8), intent(inout), optional :: dSdt(:, :, :)
#endif
    type(Particle), allocatable :: Buffers(:)

    integer i

    do i = 1, Field%Maxn
        indvdt(:, i) = 0
        indedt(i)    = 0
        exdvdt(:, i) = 0
        avdvdt(:, i) = 0
        asdvdt(:, i) = 0
        avdedt(i)    = 0
        ahdedt(i)    = 0
    end do

    !!! Positions of dummy (boundary) particles
    if ( Config%dummy_parti_w ) call gen_dummy_particle(ntotal, ndummy, Particles)

    if ( Config%open_boundary_w) then
        if (Config%i_time_step == 1) then
            call allocateParticleList(Buffers, nbuffer, Field%Dim, 1)
        else
            allocate(Buffers, source=Particles(ntotal+1:))
        end if
        call gen_non_reflecting_bc(ntotal, Particles, nbuffer, Buffers)
    else
        nbuffer = 0
    end if

    N = ntotal + ndummy + nbuffer

#ifdef _OPENMP
    Config%chunkSize = N / Config%nthreads
#endif

    !!! Interactions parameters, calculating neighboring particles
    !!! and optimizing smoothing length
    ! if ( Config%open_boundary_w) then
    !     call BGGS(Particles(1:ntotal), Particles(1:N), skipItsSelf=.true.)
    ! else
    !     call search_particles(ntotal, Particles(1:N))
    ! end if
    call BGGS(Particles(1:ntotal), Particles(1:N), skipItsSelf=.true.)

    if ( Config%kernel_correciton_w ) call kernelGradientCorrection(ntotal, Particles)

    if ( Config%shifting_w ) call shifting(ntotal, Particles)

    !!! Density approximation or change rate
    select case (Config%pa_sph)
    case (1, 2)
        if ( Config%sum_density_w ) then
            call sum_density(ntotal, Particles(1:N))
        else
            call con_density(ntotal, Particles(1:N), Delta%Density)
        end if
    case (3)
        call con_density_riemann(Particles(1:N), Delta%Density)
    case (4)
        call sum_density_dsph(Particles(1:N))
    case default
        write(err, "(1X, A, I0)") "Error density scheme ", Config%pa_sph
        error stop
    end select

    ! call norm_density(Particles(1:ntotal))

    call detonation_wave(Config%i_time_step, Config%delta_t, Particles(1:ntotal))

    call divergence(ntotal, Particles)

    !!! Internal forces
#if SOLID
    call internal_force(ntotal, Particles, indvdt, indedt, Shear, dSdt)
#else
    call internal_force(ntotal, Particles, indvdt, indedt)
#endif

    !!! External force
    if ( Config%ex_force_w ) call external_force(ntotal, Particles, exdvdt)

    !!! Artificial viscosity
    if ( Config%arti_visc_w ) call arti_visc(ntotal, Particles, avdvdt, avdedt)

    !!! Artificial heat
    if ( Config%arti_heat_w ) call arti_heat(ntotal, Particles, ahdedt)
    
    !!! Artificial stress
    if ( Config%arti_stress_w ) call artificialStress(ntotal, Particles, asdvdt)

    !!! Calculating average velocity of each particle for avoiding penetration
    if ( Config%aver_velocity_w ) call aver_velo(ntotal, Particles, aver_v)

    !!! Calculating the neighboring particles and updating HSML
    call h_upgrade(ntotal, Particles(1:ntotal))

    !!! Convert velocity, force and energy to f and dfdt
    do i = 1, ntotal
        Delta(i)%Velocity = indvdt(:, i) + avdvdt(:, i) + exdvdt(:, i) + asdvdt(:, i)
        Delta(i)%Energy   = indedt(i)    + avdedt(i)    + ahdedt(i)
    end do

    call boundary(ntotal, Particles, Delta)

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
            call print_statistics(Particles(1:ntotal))
        end if

        ! index = 0
        ! number = 0
        ! do i = 1, ntotal
        !     if ( P(i)%Type == 2 ) then
        !         number = number + 1
        !         index(number) = i
        !     end if
        ! end do
        ! call calculate_area(number, mass(index(1:number)), rho(index(1:number)), area)
        ! write(*,"(A, I0, ES15.7)") " >> Statistics: Target Area: ", number, area

        write(*,*) ">> Information for Particle ", to_string(Config%monitor_particle), &
                    " Type ", to_string(Particles(Config%monitor_particle)%Type)
        if ( Particles(Config%monitor_particle)%Type > 100 ) then
            write(*,1001) "Internal a ", "Artificial a", "External a", "AS a", "Total a"
            do i = 1, Field%Dim
                write(*,1002) indvdt(i, Config%monitor_particle), avdvdt(i, Config%monitor_particle), &
                              exdvdt(i, Config%monitor_particle), asdvdt(i, Config%monitor_particle), &
                              Delta(Config%monitor_particle)%Velocity(i)
            end do
        else
            write(*,1001) "Internal a ", "Artificial a", "External a", "Total a"
            do i = 1, Field%Dim
                write(*,1002) indvdt(i, Config%monitor_particle), avdvdt(i, Config%monitor_particle), &
                              exdvdt(i, Config%monitor_particle), &
                              Delta(Config%monitor_particle)%Velocity(i)
            end do
        end if
    end if

    ! kpair = maxval(Particles%neighborNum)
    ! if ( kpair > Field%pairNum ) then
    !     call allocateNeighborList(Particles, Field%Dim, kpair + 3)
    ! else
    !     do i = 1, size(Particles)
    !         Particles(i)%neighborNum = 0
    !         Particles(i)%neighborList = 0
    !         Particles(i)%w = 0
    !         Particles(i)%dwdx = 0
    !     end do
    ! end if
    do i = 1, ntotal + ndummy + nbuffer
        Particles(i)%neighborNum  = 0
        Particles(i)%neighborList = 0
        Particles(i)%w            = 0
        Particles(i)%dwdx         = 0
    end do

    1001 format(A17, A14, *(A15))
    1002 format(1X, *(2X, ES13.6))

end subroutine single_step