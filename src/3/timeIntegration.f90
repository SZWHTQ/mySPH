module time_integration_m
    use tools_m

    implicit none

contains
    subroutine single_step(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml, c, &
                           tdsdt, dvdt, dedt, drhodt, aver_v, div_v, div_r)
        use, intrinsic :: iso_fortran_env, only: err => error_unit
        !$ use omp_lib
        use ctrl_dict, only: dim, maxn, max_interaction, i_time_step,delta_t,nnps,skf,sle, &
                             print_interval, monitor_particle, sum_density_w, arti_visc_w, &
                             ex_force_w, arti_heat_w, aver_velocity_w, print_statistics_w, &
                             dummy_parti_w, nthreads
        use initial_m, only: pair, neighborNum, w, dwdx
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
        integer :: niac
        real(8) :: area
        integer :: number, index(maxn)
! #if defined _OPENMP
!         integer, dimension(nthreads) :: sub_ntotal, sub_ndummy
!         integer :: sub_index(nthreads, maxn)
!         integer :: omp_ntotal, omp_ndummy, omp_niac
!         integer :: omp_max_interaction
!         integer, allocatable :: omp_index(:)
!         integer, allocatable :: omp_itype(:)
!         real(8), allocatable :: omp_x(:, :)
!         real(8), allocatable :: omp_v(:, :)
!         real(8), allocatable :: omp_mass(:)
!         real(8), allocatable :: omp_rho(:)
!         real(8), allocatable :: omp_p(:)
!         real(8), allocatable :: omp_e(:)
!         real(8), allocatable :: omp_hsml(:)
!         real(8), allocatable :: omp_c(:)
!         integer, allocatable :: omp_pair(:, :)
!         real(8), allocatable :: omp_w(:), omp_dwdx(:, :)
!         real(8), allocatable :: omp_tdsdt(:), omp_drhodt(:), &
!                                 omp_aver_v(:, :), omp_div_v(:), omp_div_r(:)
!         real(8), allocatable :: omp_indvdt(:, :), omp_indedt(:), omp_exdvdt(:, :), &
!                                 omp_avdvdt(:, :), omp_avdedt(:), omp_ahdedt(:)
!         integer, allocatable :: omp_countiac(:)
!         integer :: scale_k
!         integer :: tid
! #endif
        integer i

        do i = 1, maxn
            indvdt(:, i) = 0
            indedt(i)    = 0
            exdvdt(:, i) = 0
            avdvdt(:, i) = 0
            avdedt(i)    = 0
            ahdedt(i)    = 0
            pair(i, :)     = 0
            neighborNum(i) = 0
            w(i, :)        = 0
            dwdx(:, i, :)  = 0
        end do

        !!! Positions of dummy (boundary) particles
        if ( dummy_parti_w ) call gen_dummy_particle(ntotal, ndummy, &
                                                     itype, x, v, mass, rho, p, e, c, hsml)

#if defined _OPENMP
        call omp_set_num_threads(nthreads)
#endif

! #if defined _OPENMP

!         select case (skf)
!         case (1)
!             scale_k = 2
!         case (2, 3)
!             scale_k = 3
!         end select

!         call omp_set_num_threads(nthreads)
!         call decompose(itype(1:ntotal+ndummy), x(:, 1:ntotal+ndummy), 1.2*maxval(hsml)*scale_k, nthreads, &
!                        sub_ntotal, sub_ndummy, sub_index)

!         !$OMP PARALLEL PRIVATE(tid, i) &
!         !$OMP PRIVATE(omp_ntotal, omp_ndummy) &
!         !$OMP PRIVATE(omp_niac, omp_max_interaction) &
!         !$OMP PRIVATE(omp_index) &
!         !$OMP PRIVATE(omp_itype) &
!         !$OMP PRIVATE(omp_x, omp_v, omp_mass, omp_rho) &
!         !$OMP PRIVATE(omp_p, omp_e, omp_hsml, omp_c) &
!         !$OMP PRIVATE(omp_pair, omp_w, omp_dwdx) &
!         !$OMP PRIVATE(omp_tdsdt, omp_drhodt, omp_aver_v, omp_div_v, omp_div_r) &
!         !$OMP PRIVATE(omp_indvdt, omp_indedt, omp_exdvdt) &
!         !$OMP PRIVATE(omp_avdvdt, omp_avdedt, omp_ahdedt) &
!         !$OMP PRIVATE(omp_countiac)
!         tid = omp_get_thread_num()
!         omp_ntotal = sub_ntotal(tid+1)
!         omp_ndummy = sub_ndummy(tid+1)
!         omp_max_interaction = max_interaction / maxn * (omp_ntotal+omp_ndummy)
!         allocate(omp_index, source=sub_index(tid+1, 1:omp_ntotal+omp_ndummy))
!         allocate(omp_itype, source=itype(omp_index))
!         allocate(omp_x,     source=x(:,omp_index))
!         allocate(omp_v,     source=v(:,omp_index))
!         allocate(omp_mass,  source=mass(omp_index))
!         allocate(omp_rho,   source=rho(omp_index))
!         allocate(omp_p,     source=p(omp_index))
!         allocate(omp_e,     source=e(omp_index))
!         allocate(omp_hsml,  source=hsml(omp_index))
!         allocate(omp_c,     source=c(omp_index))
!         allocate(omp_pair(omp_max_interaction, 2),   source=0)
!         allocate(omp_w(omp_max_interaction),         source=0._8)
!         allocate(omp_dwdx(dim, omp_max_interaction), source=0._8)
!         allocate(omp_tdsdt(omp_ntotal+omp_ndummy),       source=0._8)
!         allocate(omp_drhodt(omp_ntotal+omp_ndummy),      source=0._8)
!         allocate(omp_aver_v(dim, omp_ntotal+omp_ndummy), source=0._8)
!         allocate(omp_div_v(omp_ntotal+omp_ndummy),       source=0._8)
!         allocate(omp_div_r(omp_ntotal+omp_ndummy),       source=0._8)
!         allocate(omp_indvdt(dim, omp_ntotal+omp_ndummy), source=0._8)
!         allocate(omp_indedt(omp_ntotal+omp_ndummy),      source=0._8)
!         allocate(omp_exdvdt(dim, omp_ntotal+omp_ndummy), source=0._8)
!         allocate(omp_avdvdt(dim, omp_ntotal+omp_ndummy), source=0._8)
!         allocate(omp_avdedt(omp_ntotal+omp_ndummy),      source=0._8)
!         allocate(omp_ahdedt(omp_ntotal+omp_ndummy),      source=0._8)

!         !!! Interactions parameters, calculating neighboring particles
!         !!! and optimizing smoothing length
!         allocate(omp_countiac(omp_ntotal+omp_ndummy), source=0)
!         call search_particles(nnps, omp_ntotal+omp_ndummy, &
!                               x(:, omp_index), hsml(omp_index), &
!                               omp_max_interaction, omp_niac, &
!                               omp_pair, omp_w, omp_dwdx, omp_countiac)

!         !!! Density approximation or change rate
!         select case (pa_sph)
!         case (1, 2)
!             if ( sum_density_w ) then
!                 call sum_density(omp_ntotal+omp_ndummy, omp_niac, &
!                                  omp_mass, omp_rho, omp_hsml, &
!                                  omp_pair, omp_w)
!             else
!                 call con_density(omp_ntotal+omp_ndummy, omp_niac, &
!                                  omp_v, omp_mass, &
!                                  omp_pair, omp_dwdx, omp_drhodt)
!             end if
!         case (3)
!             call con_density_riemann(omp_ntotal+omp_ndummy, omp_niac, &
!                                      omp_x, omp_v, omp_mass, omp_rho, omp_p, omp_c,  &
!                                      omp_pair, omp_dwdx, omp_drhodt)
!         case (4)
!             call sum_density_dsph(omp_ntotal+omp_ndummy, omp_niac, &
!                                   omp_mass, omp_rho, omp_hsml, &
!                                   omp_pair, omp_w)
!         case default
!             write(err, "(1X, A, I0)") "Error density scheme ", pa_sph
!             error stop
!         end select

!         call detonation_wave(omp_ntotal, i_time_step, delta_t, omp_x, omp_itype)

!         ! if ( nick == "undex_cylinder" ) call free_surface()

!         call divergence(omp_ntotal, omp_niac, &
!                       omp_v, omp_mass, omp_rho, &
!                       omp_pair, omp_dwdx, &
!                       omp_div_v)

!         !!! Internal forces
!         call in_force(omp_ntotal+omp_ndummy, omp_niac, &
!                       omp_itype, omp_x, omp_v, omp_mass, &
!                       omp_rho, omp_p, omp_e, omp_c, omp_pair, omp_dwdx, &
!                       omp_indvdt, omp_tdsdt, omp_indedt)

!         !!! Artificial viscosity
!         if ( arti_visc_w ) call arti_visc(omp_ntotal+omp_ndummy, omp_niac, &
!                                           omp_x, omp_v, &
!                                           omp_mass, omp_rho, omp_c, omp_hsml, &
!                                           omp_pair, omp_dwdx, &
!                                           omp_avdvdt, omp_avdedt)

!         !!! External force
!         if ( ex_force_w ) call ex_force(omp_ntotal+omp_ndummy, omp_niac, &
!                                         omp_itype, omp_x, &
!                                         omp_hsml, &
!                                         omp_pair, &
!                                         omp_exdvdt)

!         if ( arti_heat_w ) call arti_heat(omp_ntotal+omp_ndummy, omp_niac, &
!                                           omp_x, omp_mass, &
!                                           omp_rho, omp_e, omp_c, omp_hsml, &
!                                           omp_pair, omp_dwdx, &
!                                           omp_div_v, omp_ahdedt)

!         !!! Calculating average velocity of each particle for avoiding penetration
!         if ( aver_velocity_w ) call aver_velo(omp_ntotal, omp_niac, &
!                                               omp_v, omp_mass, omp_rho, &
!                                               omp_pair, omp_w, &
!                                               omp_aver_v)

!         !!! Calculating the neighboring particles and updating HSML
!         call h_upgrade(omp_ntotal, sle, delta_t, omp_mass, omp_rho, omp_div_v, omp_hsml)

!         call divergence(omp_ntotal, omp_niac, &
!                         omp_x, omp_mass, omp_rho, &
!                         omp_pair, omp_dwdx, &
!                         omp_div_r)

!         ! do i = 1, omp_ntotal
!         !     if ( omp_index(i) == 70 ) write(*,*) omp_x(:, i)
!         ! end do

!         do i = 1, omp_ntotal
!             rho(omp_index(i))  = omp_rho(i)
!             p(omp_index(i))    = omp_p(i)
!             e(omp_index(i))    = omp_e(i)
!             hsml(omp_index(i)) = omp_hsml(i)
!             c(omp_index(i))    = omp_c(i)
!             tdsdt(omp_index(i))     = omp_tdsdt(i)
!             drhodt(omp_index(i))    = omp_drhodt(i)
!             aver_v(:, omp_index(i)) = omp_aver_v(:,i)
!             div_v(omp_index(i))     = omp_div_v(i)
!             div_r(omp_index(i))     = omp_div_r(i)
!             indvdt(:, omp_index(i)) = omp_indvdt(:,i)
!             indedt(omp_index(i))    = omp_indedt(i)
!             exdvdt(:, omp_index(i)) = omp_exdvdt(:,i)
!             avdvdt(:, omp_index(i)) = omp_avdvdt(:,i)
!             avdedt(omp_index(i))    = omp_avdedt(i)
!             ahdedt(omp_index(i))    = omp_ahdedt(i)
!             neighborNum(omp_index(i))  = omp_countiac(i)
!         end do
!         !$OMP END PARALLEL
!         ! write(*,*) "x:", x(:,70)
! #else
        !!! Interactions parameters, calculating neighboring particles
        !!! and optimizing smoothing length
        ! write(*,*) niac
        call search_particles(nnps, x(:, 1:ntotal+ndummy), hsml(1:ntotal+ndummy), &
                              pair, neighborNum, w, dwdx)
        ! write(*,*) niac

        !!! Density approximation or change rate
        select case (pa_sph)
        case (1, 2)
            if ( sum_density_w ) then
                call sum_density(ntotal+ndummy, mass, rho, hsml, &
                                 pair, neighborNum, w)
            else
                call con_density(ntotal+ndummy, v, mass, &
                                 pair, neighborNum, dwdx, drhodt)
            end if
        case (3)
            call con_density_riemann(ntotal+ndummy, &
                                     x, v, mass, rho, p, c,  &
                                     pair, neighborNum, dwdx, drhodt)
        case (4)
            call sum_density_dsph(ntotal+ndummy, mass, rho, hsml, &
                                  pair, neighborNum, w)
        case default
            write(err, "(1X, A, I0)") "Error density scheme ", pa_sph
            error stop
        end select

        call detonation_wave(ntotal, i_time_step, delta_t, x, itype)

        ! if ( nick == "undex_cylinder" ) call free_surface()

        call divergence(ntotal+ndummy, v, mass, rho, pair, neighborNum, dwdx, div_v)

        !!! Internal forces
        call in_force(ntotal+ndummy, itype, x, v, mass, rho, &
                      p, e, c, pair, neighborNum, dwdx, &
                      indvdt, tdsdt, indedt)

        !!! Artificial viscosity
        if ( arti_visc_w ) call arti_visc(ntotal+ndummy, x, v, &
                                          mass, rho, c, hsml,  &
                                          pair, dwdx, neighborNum, avdvdt, avdedt)

        !!! External force
        if ( ex_force_w ) call ex_force(ntotal+ndummy, itype, x, &
                                        hsml, pair, neighborNum, exdvdt)

        if ( arti_heat_w ) call arti_heat(ntotal+ndummy, x, mass, rho, e, c, hsml,  &
                                          pair, neighborNum, dwdx, div_v, ahdedt)

        !!! Calculating average velocity of each particle for avoiding penetration
        if ( aver_velocity_w ) call aver_velo(ntotal, v, mass, rho, &
                                              pair, neighborNum, w, aver_v)

        !!! Calculating the neighboring particles and updating HSML
        call h_upgrade(ntotal, sle, delta_t, mass, rho, div_v, hsml)

        call divergence(ntotal+ndummy, x, mass, rho, pair, neighborNum, dwdx, div_r)

! #endif

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
                call print_statistics(ntotal, neighborNum)
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
        use ctrl_dict, only: dim, maxn, &
                             i_time_step, max_time_step, &
                             ntotal, ndummy, &
                             save_interval, print_interval, monitor_particle, &
                             sum_density_w, nsym
        use initial_m
        use cour_num_m
        use output_m
        implicit none
        real(8) :: x_min(dim, maxn), v_min(dim, maxn), &
                   e_min(maxn), rho_min(maxn)
        real(8) :: tdsdt(maxn), dvdt(dim, maxn), dedt(maxn), &
                   drhodt(maxn), aver_v(dim, maxn), &
                   div_v(maxn), div_r(maxn)
        real(8) :: time = 0, temp_rho, temp_e
        real(8) :: aver_courant = 0, max_courant = 0, cntemp

        integer i

        forall (i=1:maxn)
            x_min(:, i)  = 0
            v_min(:, i)  = 0
            e_min(i)     = 0
            rho_min(i)   = 0
            tdsdt(i)     = 0
            dvdt(:, i)   = 0
            dedt(i)      = 0
            drhodt(i)    = 0
            aver_v(:, i) = 0
        end forall

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
                !!$OMP PARALLEL DO PRIVATE(i) NUM_THREADS(nthreads)
                do i = 1, ntotal
                    e_min(i) = e(i)
                    temp_e = 0
                    if ( dim == 1 ) temp_e = -nsym*p(i)*v(1, i)/x(1, i)/rho(i)
                    e(i) = e(i) + (delta_t/2) * (dedt(i)+temp_e)
                    if ( e(i) < 0 ) e(i) = 0
                    if ( .not. sum_density_w ) then
                        rho_min(i) = rho(i)
                        temp_rho = 0
                        if ( dim == 1 ) temp_rho = -nsym*rho(i)*v(1, i)/x(1, i)
                        rho(i) = rho(i) + (delta_t/2) * (drhodt(i) + temp_rho)
                    end if
                    v_min(:, i) = v(:, i)
                    v(:, i) = v(:, i) + (delta_t/2)*dvdt(:, i)
                end do
                !!$OMP END PARALLEL DO
            end if

            !!! Definition of variables out of the function vector:
            call single_step(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml, c, &
                             tdsdt, dvdt, dedt, drhodt, aver_v, div_v, div_r)

            if ( i_time_step == 1 ) then
                !!$OMP PARALLEL DO PRIVATE(i) NUM_THREADS(nthreads)
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
                end do
                !!$OMP END PARALLEL DO
            else
                max_courant = 0
                !!$OMP PARALLEL DO PRIVATE(i) NUM_THREADS(nthreads)
                do i = 1, ntotal

                    temp_e = 0._8
                    if ( dim == 1 ) temp_e = -nsym*p(i)*v(1, i)/x(1, i)/rho(i)

                    e(i) = e_min(i) + delta_t * (dedt(i)+temp_e)
                    if ( e(i) < 0 ) e(i) = 0

                    if ( .not. sum_density_w ) then
                        temp_rho = 0
                        if ( dim == 1 ) temp_rho = -nsym*rho(i)*v(1, i)/x(1, i)

                        rho(i) = rho_min(i) + delta_t * (drhodt(i)+temp_rho)
                    end if

                    v(:, i) = v_min(:, i) + delta_t * dvdt(:, i) + aver_v(:, i)
                    x(:, i) = x(:, i) + delta_t * v(:, i)
                    cntemp = courant_num(hsml(i), div_v(i), c(i))
                    aver_courant = aver_courant + cntemp
                    if ( cntemp > max_courant ) max_courant = cntemp
                end do
                !!$OMP END PARALLEL DO
                aver_courant = aver_courant / ntotal
            end if

            time = time + delta_t

            if (mod(i_time_step, save_interval) == 0) then
                call output(i_time_step, ntotal, ndummy, itype, x, v, mass, rho, p, e, c, hsml, div_r)
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