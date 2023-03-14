#include "../macro.h"
module time_integration_m
    use tools_m

    implicit none

contains
! #if SOLID
    subroutine single_step(ntotal, ndummy, Particles, drhodt, dvdt, dedt, tdsdt, aver_v, Shear, dSdt)
! #else
!     subroutine single_step(ntotal, ndummy, Particles, drhodt, dvdt, dedt, tdsdt, aver_v)
! #endif
        use, intrinsic :: iso_fortran_env, only: err => error_unit
        !$ use omp_lib
#ifndef _OPENMP
        use ctrl_dict, only: dim, maxn, i_time_step, delta_t,nnps,sle, print_interval, &
                             monitor_particle, sum_density_w, arti_visc_w, ex_force_w, &
                             arti_heat_w, aver_velocity_w, print_statistics_w, dummy_parti_w
#else
        use ctrl_dict, only: dim, maxn, i_time_step, delta_t,nnps,sle, print_interval,        &
                             monitor_particle, sum_density_w, arti_visc_w, ex_force_w,        &
                             arti_heat_w, aver_velocity_w, print_statistics_w, dummy_parti_w, &
                             nthreads, chunkSize
#endif
        use sph
        use nnps_m
        use density_m
        use visc_m
        use divergence_m
        use in_force_m
        use arti_visc_m
        use ex_force_m
        use hsml_m
        use arti_heat_m
        use corr_velo_m
        use dummy_part_m
        use he_m
        use decompose_m
        use area_m
        ! use boundary_condition_m
        integer, intent(in)    :: ntotal
        integer, intent(inout) :: ndummy
        type(Particle), intent(inout) :: Particles(:)
        real(8), intent(inout) :: drhodt(:), dvdt(:, :), dedt(:), &
                                  tdsdt(:), aver_v(:, :)
        real(8) :: indvdt(dim, maxn), indedt(maxn), exdvdt(dim, maxn), &
                   avdvdt(dim, maxn), avdedt(maxn), ahdedt(maxn)
        ! real(8) :: area
        ! integer :: number, index(maxn)
#if SOLID
        real(8), intent(in),    optional :: Shear(:, :, :)
        real(8), intent(inout), optional :: dSdt(:, :, :)
#endif

        integer i

        do i = 1, maxn
            indvdt(:, i) = 0
            indedt(i)    = 0
            exdvdt(:, i) = 0
            avdvdt(:, i) = 0
            avdedt(i)    = 0
            ahdedt(i)    = 0
        end do

        !!! Positions of dummy (boundary) particles
        if ( dummy_parti_w ) call gen_dummy_particle(ndummy, Particles(1:ntotal))

#ifdef _OPENMP
        chunkSize = (ntotal+ndummy) / nthreads
#endif
        !!! Interactions parameters, calculating neighboring particles
        !!! and optimizing smoothing length
        ! write(*,*) niac
        call search_particles(nnps, Particles(1:ntotal+ndummy))
        ! write(*,*) niac

        !!! Density approximation or change rate
        select case (pa_sph)
        case (1, 2)
            if ( sum_density_w ) then
                call sum_density(Particles(1:ntotal+ndummy))
            else
                call con_density(Particles(1:ntotal+ndummy), drhodt)
            end if
        case (3)
            call con_density_riemann(Particles(1:ntotal+ndummy), drhodt)
        case (4)
            call sum_density_dsph(Particles(1:ntotal+ndummy))
        case default
            write(err, "(1X, A, I0)") "Error density scheme ", pa_sph
            error stop
        end select

        call detonation_wave(i_time_step, delta_t, Particles(1:ntotal))

        call divergence(Particles(1:ntotal+ndummy))

        !!! Internal forces
#if SOLID
        call in_force(Particles(1:ntotal+ndummy), indvdt, tdsdt, indedt, Shear, dSdt)
#else
        call in_force(Particles(1:ntotal+ndummy), indvdt, tdsdt, indedt)
#endif

        !!! Artificial viscosity
        if ( arti_visc_w ) call arti_visc(Particles(1:ntotal+ndummy), avdvdt, avdedt)

        !!! External force
        if ( ex_force_w ) call ex_force(Particles(1:ntotal+ndummy), exdvdt)

        if ( arti_heat_w ) call arti_heat(Particles(1:ntotal+ndummy), ahdedt)

        !!! Calculating average velocity of each particle for avoiding penetration
        if ( aver_velocity_w ) call aver_velo(Particles(1:ntotal), aver_v)

        !!! Calculating the neighboring particles and updating HSML
        call h_upgrade(sle, delta_t, Particles(1:ntotal))

        !!! Convert velocity, force and energy to f and dfdt
        do i = 1, ntotal
            dvdt(:, i) = indvdt(:, i) + avdvdt(:, i) + exdvdt(:, i)
            dedt(i)    = indedt(i)    + avdedt(i)    + ahdedt(i)
        end do


        if ( mod(i_time_step, print_interval) == 0 ) then
            !!! Statistics for the interaction
            if (print_statistics_w) then
                write(*,*) ">> Statistics: Dummy particles:"
                write(*,*) "   Number of dummy particles: ", to_string(ndummy)
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

            write(*,*) ">> Information for particle ", to_string(monitor_particle)
            write(*,1001) "Internal a ", "Artificial a", "External a", "Total a"
            do i = 1, dim
                write(*,1002) indvdt(i, monitor_particle), avdvdt(i, monitor_particle), &
                              exdvdt(i, monitor_particle),   dvdt(i, monitor_particle)
            end do
        end if

        call allocateNeighborList(Particles, dim, maxval(Particles%neighborNum)+5)

        1001 format(A17, A14, 2(A15))
        1002 format(1X, 4(2X, ES13.6))

    end subroutine single_step

    subroutine time_integration(ntotal, P)
#ifndef _OPENMP
        use ctrl_dict, only: dim, maxn, i_time_step, max_time_step, &
                             save_interval, print_interval,         &
                             monitor_particle, sum_density_w, nsym
#else
        use ctrl_dict, only: dim, maxn, i_time_step, max_time_step, &
                             save_interval, print_interval,         &
                             monitor_particle, sum_density_w, nsym, &
                             nthreads, chunkSize
#endif
        use sph
        use cour_num_m
        use output_m
        implicit none
        integer, intent(in) :: ntotal
        type(Particle), intent(inout) :: P(:)
        integer :: ndummy = 0
        real(8) :: v_prev(dim, maxn), e_prev(maxn), rho_prev(maxn)
        real(8) :: tdsdt(maxn), dvdt(dim, maxn), dedt(maxn), &
                   drhodt(maxn), aver_v(dim, maxn)
        real(8) :: temp_rho, temp_e
        real(8) :: time = 0
        real(8) :: aver_courant = 0, max_courant = 0, cntemp
#if SOLID
        real(8) :: Shear_prev(dim, dim, maxn)
        real(8) :: Shear(dim, dim, maxn), dSdt(dim, dim, maxn)
        real(8) :: J2, SigmaY
#endif

        integer i

        forall (i=1:maxn)
            v_prev(:, i) = 0
            e_prev(i)    = 0
            rho_prev(i)  = 0
            tdsdt(i)     = 0
            dvdt(:, i)   = 0
            dedt(i)      = 0
            drhodt(i)    = 0
            aver_v(:, i) = 0
#if SOLID
            Shear_prev(:, :, i) = 0
            Shear(:, :, i)      = 0
            dSdt(:, :, i)       = 0
#endif
        end forall

#if SOLID
        SigmaY = 6e8
#endif

#ifdef _OPENMP
        call omp_set_num_threads(nthreads)
#endif

        call pbout(0, max_time_step, .true.)
        call pbflush()

        do i_time_step = 1, max_time_step

            if ( mod(i_time_step, print_interval) == 0 ) then
                call pbflush()
#ifndef _WIN32
                write(*, "(A)") repeat("—", 72)
#else
                write(*, "(A)") repeat("-", 72)
#endif
                write(*, "(2(A, G0))") " Courant Number mean: ", aver_courant, &
                                       " max: ", max_courant
                write(*,*) "Time step = ", to_string(i_time_step)
                write(*, "(A, G0, A)") " deltaT = ", delta_t, "s"
                write(*, "(A, G0, A)") " Time   = ", time+delta_t, "s"
            end if

            !!! If not first time step, then update thermal energy, density
            !!! and velocity half a time step
            if ( i_time_step /= 1 ) then
                !$OMP PARALLEL DO PRIVATE(i)
                do i = 1, ntotal
                    e_prev(i) = P(i)%InternalEnergy
                    temp_e = 0
                    if ( dim == 1 ) temp_e = -nsym &
                        * P(i)%Pressure * P(i)%v(1) / P(i)%x(1) / P(i)%Density
                    P(i)%InternalEnergy = P(i)%InternalEnergy + (delta_t/2) * (dedt(i)+temp_e)
                    if ( P(i)%InternalEnergy < 0 ) P(i)%InternalEnergy = 0
                    if ( .not. sum_density_w ) then
                        rho_prev(i) = P(i)%Density
                        temp_rho = 0
                        if ( dim == 1 ) temp_rho = -nsym &
                            * P(i)%Density * P(i)%v(1) / P(i)%x(1)
                        P(i)%Density = P(i)%Density + (delta_t/2) * (drhodt(i) + temp_rho)
                    end if
                    v_prev(:, i) = P(i)%v(:)
                    P(i)%v(:) = P(i)%v(:) + (delta_t/2)*dvdt(:, i)
#if SOLID
                    if ( abs(P(i)%Type) == 8 ) then
                        Shear_prev(:, :, i) = Shear(:, :, i)
                        Shear(:, :, i) = Shear(:, :, i) + (delta_t/2)*dSdt(:, :, i)
                        J2 = sum( Shear(:, :, i)**2 )
                        Shear(:, :, i) = Shear(:, :, i) * min(1., sqrt(((SigmaY**2)/3)/J2))
                    end if
#endif
                end do
                !$OMP END PARALLEL DO
            end if

#if SOLID
            call single_step(ntotal, ndummy, P, drhodt, dvdt, dedt, tdsdt, aver_v, Shear, dSdt)
#else
            call single_step(ntotal, ndummy, P, drhodt, dvdt, dedt, tdsdt, aver_v)
#endif

            if ( i_time_step == 1 ) then
                !$OMP PARALLEL DO PRIVATE(i)
                do i = 1, ntotal
                    temp_e = 0
                    if ( dim == 1 ) temp_e = -nsym*P(i)%Pressure*P(i)%v(1)/P(i)%x(1)/P(i)%Density

                    P(i)%InternalEnergy = P(i)%InternalEnergy + (delta_t/2) * (dedt(i) + temp_e)
                    if ( P(i)%InternalEnergy < 0 ) P(i)%InternalEnergy = 0

                    if ( .not. sum_density_w ) then
                        temp_rho = 0
                         if ( dim == 1 ) temp_rho = -nsym*P(i)%Density*P(i)%v(1)/P(i)%x(1)

                        P(i)%Density = P(i)%Density + (delta_t/2) * (drhodt(i)+temp_rho)
                    end if

                    P(i)%v(:) = P(i)%v(:) + (delta_t/2) * dvdt(:, i) + aver_v(:, i)
                    P(i)%x(:) = P(i)%x(:) + delta_t * P(i)%v(:)
#if SOLID
                    if ( abs(P(i)%Type) == 8 ) then
                        Shear(:, :, i) = Shear(:, :, i) + (delta_t/2) * dSdt(:, :, i)
                        J2 = sum( Shear(:, :, i)**2 )
                        Shear(:, :, i) = Shear(:, :, i) * min(1., sqrt(((SigmaY**2)/3)/J2))
                    end if
#endif
                end do
                !$OMP END PARALLEL DO
            else
                max_courant = 0
                !$OMP PARALLEL DO PRIVATE(i) REDUCTION(max:max_courant) REDUCTION(+:aver_courant)
                do i = 1, ntotal
                    temp_e = 0._8
                    if ( dim == 1 ) temp_e = -nsym*P(i)%Pressure*P(i)%v(1)/P(i)%x(1)/P(i)%Density

                    P(i)%InternalEnergy = e_prev(i) + delta_t * (dedt(i)+temp_e)
                    if ( P(i)%InternalEnergy < 0 ) P(i)%InternalEnergy = 0

                    if ( .not. sum_density_w ) then
                        temp_rho = 0
                        if ( dim == 1 ) temp_rho = -nsym*P(i)%Density*P(i)%v(1)/P(i)%x(1)

                        P(i)%Density = rho_prev(i) + delta_t * (drhodt(i)+temp_rho)
                    end if

                    P(i)%v(:) = v_prev(:, i) + delta_t * dvdt(:, i) + aver_v(:, i)
                    P(i)%x(:) = P(i)%x(:) + delta_t * P(i)%v(:)
#if SOLID
                    if ( abs(P(i)%Type) == 8 ) then
                        Shear(:, :, i) = Shear_prev(:, :, i) + delta_t * dSdt(:, :, i)
                        J2 = sum( Shear(:, :, i)**2 )
                        Shear(:, :, i) = Shear(:, :, i) * min(1., sqrt(((SigmaY**2)/3)/J2))
                    end if
#endif
                    cntemp = courant_num(P(i)%SmoothingLength, P(i)%divergenceVelocity, P(i)%SoundSpeed)
                    aver_courant = aver_courant + cntemp
                    if ( cntemp > max_courant ) max_courant = cntemp
                end do
                !$OMP END PARALLEL DO
                aver_courant = aver_courant / ntotal
            end if

            time = time + delta_t

            if (mod(i_time_step, save_interval) == 0) then
                call output((i_time_step/save_interval), P(1:ntotal+ndummy))
            end if

            if ( mod(i_time_step, print_interval) == 0 ) then
                write(*,1000) "Location", "Velocity", "Acceleration"
                do i = 1, dim
                    write(*,1001) P(monitor_particle)%x(i), &
                                  P(monitor_particle)%v(i), &
                                  dvdt(i, monitor_particle)
                end do
#ifndef _WIN32
                write(*, "(A)") repeat("—", 72)
#else
                write(*, "(A)") repeat("-", 72)
#endif
                write(*,*)
                call pbout(i_time_step, max_time_step, .true.)
            end if
        1000    format(1X, 3(A15))
        1001    format(1X, 3(2X, ES13.6))

        end do

    end subroutine time_integration

end module time_integration_m