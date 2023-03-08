#include "../macro.h"
module time_integration_m
    use tools_m

    implicit none

contains
    subroutine single_step(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml, c, &
                           tdsdt, dvdt, dedt, drhodt, aver_v, div_v, div_r, Stress, dSdt)
        use, intrinsic :: iso_fortran_env, only: err => error_unit
        !$ use omp_lib
        use ctrl_dict, only: dim, maxn, i_time_step, delta_t,nnps,sle, print_interval, &
                             monitor_particle, sum_density_w, arti_visc_w, ex_force_w, &
                             arti_heat_w, aver_velocity_w, print_statistics_w, dummy_parti_w
        use initial_m, only: neighborList, neighborNum, w, dwdx
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
        integer, intent(inout) :: itype(:)
        real(8), intent(inout) :: x(:, :)
        real(8), intent(inout) :: v(:, :)
        real(8), intent(inout) :: mass(:)
        real(8), intent(inout) :: rho(:)
        real(8), intent(inout) :: p(:)
        real(8), intent(inout) :: e(:)
        real(8), intent(inout) :: hsml(:)
        real(8), intent(inout) :: c(:)
        real(8), intent(inout) :: tdsdt(:), dvdt(:, :), dedt(:), &
                                  drhodt(:), aver_v(:, :), &
                                  div_v(:), div_r(:)
        real(8) :: indvdt(dim, maxn), indedt(maxn), exdvdt(dim, maxn), &
                   avdvdt(dim, maxn), avdedt(maxn), ahdedt(maxn)
        real(8) :: area
        integer :: number, index(maxn)
#if SOLID
        real(8), intent(in),    optional :: Stress(:, :, :)
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
            neighborNum(i)     = 0
            neighborList(i, :) = 0
            w(i, :)       = 0
            dwdx(:, i, :) = 0
        end do

        !!! Positions of dummy (boundary) particles
        if ( dummy_parti_w ) call gen_dummy_particle(ntotal, ndummy, &
                                                     itype, x, v, mass, rho, p, e, c, hsml)

        !!! Interactions parameters, calculating neighboring particles
        !!! and optimizing smoothing length
        ! write(*,*) niac
        call search_particles(nnps, x(:, 1:ntotal+ndummy), hsml(1:ntotal+ndummy), &
                              neighborNum, neighborList, w, dwdx)
        ! write(*,*) niac

        !!! Density approximation or change rate
        select case (pa_sph)
        case (1, 2)
            if ( sum_density_w ) then
                call sum_density(ntotal+ndummy, mass, rho, hsml, &
                                 neighborNum, neighborList, w)
            else
                call con_density(ntotal+ndummy, v, mass, &
                                 neighborNum, neighborList, dwdx, drhodt)
            end if
        case (3)
            call con_density_riemann(ntotal+ndummy, &
                                     x, v, mass, rho, p, c,  &
                                     neighborNum, neighborList, dwdx, drhodt)
        case (4)
            call sum_density_dsph(ntotal+ndummy, mass, rho, hsml, &
                                  neighborNum, neighborList, w)
        case default
            write(err, "(1X, A, I0)") "Error density scheme ", pa_sph
            error stop
        end select

        call detonation_wave(ntotal, i_time_step, delta_t, x, itype)

        ! if ( nick == "undex_cylinder" ) call free_surface()

        call divergence(ntotal+ndummy, v, mass, rho, neighborNum, neighborList, dwdx, div_v)

        !!! Internal forces
#if SOLID
        call in_force(ntotal+ndummy, itype, x, v, mass, rho, p, e, c, &
                      neighborNum, neighborList, dwdx, indvdt, tdsdt, indedt, Stress, dSdt)
#else
        call in_force(ntotal+ndummy, itype, x, v, mass, rho, p, e, c, &
                      neighborNum, neighborList, dwdx, indvdt, tdsdt, indedt)
#endif

        !!! Artificial viscosity
        if ( arti_visc_w ) call arti_visc(ntotal+ndummy, x, v, mass, rho, c, hsml, &
                                          neighborNum, neighborList, dwdx, avdvdt, avdedt)

        !!! External force
        if ( ex_force_w ) call ex_force(ntotal+ndummy, itype, x, hsml, &
                                        neighborNum, neighborList, exdvdt)

        if ( arti_heat_w ) call arti_heat(ntotal+ndummy, x, mass, rho, e, c, hsml, &
                                          neighborNum, neighborList, dwdx, div_v, ahdedt)

        !!! Calculating average velocity of each particle for avoiding penetration
        if ( aver_velocity_w ) call aver_velo(ntotal, v, mass, rho, &
                                              neighborNum, neighborList, w, aver_v)

        !!! Calculating the neighboring particles and updating HSML
        call h_upgrade(ntotal, sle, delta_t, mass, rho, div_v, hsml)

        call divergence(ntotal+ndummy, x, mass, rho, neighborNum, neighborList, dwdx, div_r)

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
                call print_statistics(ntotal+ndummy, neighborNum)
            end if
            index = 0
            number = 0
            do i = 1, ntotal+ndummy
                if ( itype(i) == 2 ) then
                    number = number + 1
                    index(number) = i
                end if
            end do
            call calculate_area(number, mass(index(1:number)), rho(index(1:number)), area)
            write(*,"(A, I0, ES15.7)") " >> Statistics: Target Area: ", number, area
            write(*,*) ">> Information for particle ", to_string(monitor_particle)
            write(*,1001) "Internal a ", "Artificial a", "External a", "Total a"
            do i = 1, dim
                write(*,1002) indvdt(i, monitor_particle), avdvdt(i, monitor_particle), &
                              exdvdt(i, monitor_particle),   dvdt(i, monitor_particle)
            end do
        end if

        1001 format(A17, A14, 2(A15))
        1002 format(1X, 4(2X, ES13.6))

    end subroutine single_step

    subroutine time_integration()
#ifndef _OPENMP
        use ctrl_dict, only: dim, maxn, i_time_step, max_time_step, ntotal,   &
                             save_interval, print_interval, monitor_particle, &
                             sum_density_w, nsym
#else
        use ctrl_dict, only: dim, maxn, i_time_step, max_time_step, ntotal,   &
                             save_interval, print_interval, monitor_particle, &
                             sum_density_w, nsym, nthreads, chunkSize
#endif
        use initial_m
        use cour_num_m
        use output_m
        implicit none
        integer :: ndummy
        real(8) :: v_prev(dim, maxn), e_prev(maxn), rho_prev(maxn)
        real(8) :: tdsdt(maxn), dvdt(dim, maxn), dedt(maxn), &
                   drhodt(maxn), aver_v(dim, maxn), &
                   div_v(maxn), div_r(maxn)
        real(8) :: time = 0, temp_rho, temp_e
        real(8) :: aver_courant = 0, max_courant = 0, cntemp
#if SOLID
        ! integer :: solid_num
        real(8) :: Stress_prev(dim, dim, maxn)
        real(8) :: Stress(dim, dim, maxn), dSdt(dim, dim, maxn)
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
            Stress_prev(:, :, i) = 0
            Stress(:, :, i)      = 0
            dSdt(:, :, i)        = 0
        end forall

! #if SOLID
!         ! solid_num = 0
!         ! do i = 1, ntotal
!         !     if ( itype(i) == 8 ) solid_num = solid_num + 1
!         ! end do

!         allocate(Stress(dim, dim, maxn), dSdt(dim, dim, maxn), source=0._8)
!         SigmaY = 5e8
! #endif
    
#ifdef _OPENMP
        call omp_set_num_threads(nthreads)
        chunkSize = maxn / nthreads
#endif

        do i_time_step = 1, max_time_step

            if ( mod(i_time_step, print_interval) == 0 ) then
                call pbflush()
                write(*, "(A)") repeat("—", 68)
                write(*, "(2(A, G0))") " Courant Number mean: ", aver_courant, &
                                       " max: ", max_courant
                write(*,*) "Time step = ", to_string(i_time_step)
                write(*, "(A, G0, A)") " deltaT = ", delta_t, "s"
                ! write(buffer, "(F8.3)") round(time+delta_t, 3)
                ! write(*,*) "Time   = ", trim(adjustl(buffer)), "s"
                write(*, "(A, G0, A)") " Time   = ", time+delta_t, "s"
            end if

            !!! If not first time step, then update thermal energy, density
            !!! and velocity half a time step
            if ( i_time_step /= 1 ) then
                !$OMP PARALLEL DO PRIVATE(i)
                do i = 1, ntotal
                    e_prev(i) = e(i)
                    temp_e = 0
                    if ( dim == 1 ) temp_e = -nsym*p(i)*v(1, i)/x(1, i)/rho(i)
                    e(i) = e(i) + (delta_t/2) * (dedt(i)+temp_e)
                    if ( e(i) < 0 ) e(i) = 0
                    if ( .not. sum_density_w ) then
                        rho_prev(i) = rho(i)
                        temp_rho = 0
                        if ( dim == 1 ) temp_rho = -nsym*rho(i)*v(1, i)/x(1, i)
                        rho(i) = rho(i) + (delta_t/2) * (drhodt(i) + temp_rho)
                    end if
                    v_prev(:, i) = v(:, i)
                    v(:, i) = v(:, i) + (delta_t/2)*dvdt(:, i)
#if SOLID
                    Stress_prev(:, :, i) = Stress(:, :, i)
                    Stress(:, :, i) = Stress(:, :, i) + (delta_t/2)*dSdt(:, :, i)
                    J2 = sqrt( 1.5 * sum( Stress(:, :, i)**2 ))
                    if ( J2 > SigmaY ) Stress(:, :, i) = Stress(:, :, i) * SigmaY / J2
#endif
                end do
                !$OMP END PARALLEL DO
            end if

#if SOLID
            call single_step(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml, c, &
                             tdsdt, dvdt, dedt, drhodt, aver_v, div_v, div_r, Stress, dSdt)
#else
            call single_step(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml, c, &
                             tdsdt, dvdt, dedt, drhodt, aver_v, div_v, div_r)
#endif

            if ( i_time_step == 1 ) then
                !$OMP PARALLEL DO PRIVATE(i)
                do i = 1, ntotal
                    temp_e = 0
                    if ( dim == 1 ) temp_e = -nsym*p(i)*v(1, i)/x(1, i)/rho(i)

                    e(i) = e(i) + (delta_t/2) * (dedt(i) + temp_e)
                    if ( e(i) < 0 ) e(i) = 0

                    if ( .not. sum_density_w ) then
                        temp_rho = 0
                         if ( dim == 1 ) temp_rho = -nsym*rho(i)*v(1, i)/x(1, i)

                        rho(i) = rho(i) + (delta_t/2) * (drhodt(i)+temp_rho)
                    end if

                    v(:, i) = v(:, i) + (delta_t/2) * dvdt(:, i) + aver_v(:, i)
                    x(:, i) = x(:, i) + delta_t * v(:, i)
#if SOLID
                    Stress(:, :, i) = Stress(:, :, i) + (delta_t/2) * dSdt(:, :, i)
                    J2 = sqrt( 1.5 * sum( Stress(:, :, i)**2 ))
                    if ( J2 > SigmaY ) Stress(:, :, i) = Stress(:, :, i) * SigmaY / J2
#endif
                end do
                !$OMP END PARALLEL DO
            else
                max_courant = 0
                !$OMP PARALLEL DO PRIVATE(i) REDUCTION(max:max_courant) REDUCTION(+:aver_courant)
                do i = 1, ntotal
                    temp_e = 0._8
                    if ( dim == 1 ) temp_e = -nsym*p(i)*v(1, i)/x(1, i)/rho(i)

                    e(i) = e_prev(i) + delta_t * (dedt(i)+temp_e)
                    if ( e(i) < 0 ) e(i) = 0

                    if ( .not. sum_density_w ) then
                        temp_rho = 0
                        if ( dim == 1 ) temp_rho = -nsym*rho(i)*v(1, i)/x(1, i)

                        rho(i) = rho_prev(i) + delta_t * (drhodt(i)+temp_rho)
                    end if

                    v(:, i) = v_prev(:, i) + delta_t * dvdt(:, i) + aver_v(:, i)
                    x(:, i) = x(:, i) + delta_t * v(:, i)
#if SOLID
                    Stress(:, :, i) = Stress_prev(:, :, i) + delta_t * dSdt(:, :, i)
                    J2 = sqrt( 1.5 * sum( Stress(:, :, i)**2 ))
                    if ( J2 > SigmaY ) Stress(:, :, i) = Stress(:, :, i) * SigmaY / J2
#endif
                    cntemp = courant_num(hsml(i), div_v(i), c(i))
                    aver_courant = aver_courant + cntemp
                    if ( cntemp > max_courant ) max_courant = cntemp
                end do
                !$OMP END PARALLEL DO

                ! do i = 1, ntotal
                !     cntemp = courant_num(hsml(i), div_v(i), c(i))
                !     aver_courant = aver_courant + cntemp
                !     if ( cntemp > max_courant ) max_courant = cntemp
                ! end do
                aver_courant = aver_courant / ntotal
            end if

            time = time + delta_t

            if (mod(i_time_step, save_interval) == 0) then
                call output((i_time_step/save_interval), ntotal+ndummy, itype, x, v, mass, rho, p, e, c, hsml, div_r)
            end if

            if ( mod(i_time_step, print_interval) == 0 ) then
                write(*,1000) "Location", "Velocity", "Acceleration"
                do i = 1, dim
                    write(*,1001) x(i, monitor_particle), &
                                  v(i, monitor_particle), &
                                  dvdt(i, monitor_particle)
                end do
                write(*, "(A)") repeat("—", 68)
                write(*,*)
                call pbout(i_time_step, max_time_step, .true.)
            end if
        1000    format(1X, 3(A15))
        1001    format(1X, 3(2X, ES13.6))

        end do

    end subroutine time_integration

end module time_integration_m