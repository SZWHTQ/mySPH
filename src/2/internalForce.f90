module in_force_m
    use ctrl_dict, only: dim, maxn, &
                         viscosity_w, pa_sph
    use parse_toml_m, only: nick
    ! use boundary_condition_m

    implicit none

contains
    subroutine in_force(ntotal, itype, x, v, mass, rho, p, e, c, pair, neighborNum, dwdx, dvdt, tdsdt, dedt)
        use, intrinsic :: iso_fortran_env, only: err => error_unit
        use visc_m
        use eos_m
        use initial_m, only: eta
        !$ use omp_lib
        integer, intent(in)  :: ntotal
        integer, intent(in)  :: itype(:)
        real(8), intent(in)  :: x(:, :)
        real(8), intent(in)  :: v(:, :)
        real(8), intent(in)  :: mass(:)
        real(8), intent(in)  :: rho(:)
        real(8), intent(inout) :: p(:)
        real(8), intent(in)  :: e(:)
        real(8), intent(inout) :: c(:)
        integer, intent(in)  :: pair(:, :)
        integer, intent(in)  :: neighborNum(:)
        real(8), intent(in)  :: dwdx(:, :, :)
        real(8), intent(inout) :: dvdt(:, :), tdsdt(:), dedt(:)
        real(8), allocatable :: edot(:, :, :)
        real(8) :: dv(dim)
        real(8) :: Z_l, Z_r, v_l, v_r, v_ij, v_star(dim), e_ij(dim)
        real(8) :: rhoij, aux
        real(8) :: p_star

        integer i, j, k, d, dd, ddd

        !!! Initialization of shear tensor, velocity divergence,
        !!! viscous energy, internal energy, acceleration
        allocate(edot(dim, dim, ntotal), source=0._8)

        !!! Dynamic viscosity
        if ( viscosity_w ) call viscosity(itype, eta)


        !!! Calculate SPH sum for shear tensor
        !!! = va,b + vb,a - 2/3*delta_ab*vc,c
        if ( viscosity_w ) then
                !$OMP PARALLEL DO PRIVATE(i, j, k, d, dd, ddd, dv, rhoij, aux)
                do i = 1, ntotal !! All particles
                    if ( abs(itype(i)) < 8 ) then !! Fluid
                    do k = 1, neighborNum(i) !! All neighbors of each particle
                        j = pair(i, k)
                        dv = v(:, j) - v(:, i)
                        do d = 1, dim !! All dimensions For the First Order of Strain Rate Tensor, Loop 1
                            do dd = 1, dim !! All dimensions For the Second Order of Strain Rate Tensor, Loop 2
                                edot(d, dd, i) = edot(d, dd, i) &
                                    + mass(j)/rho(j)            &
                                    * dv(d)                     &
                                    * dwdx(dd, i, k)
                                edot(dd, d, i) = edot(dd, d, i) &
                                    + mass(j)/rho(j)            &
                                    * dv(dd)                    &
                                    * dwdx(d, i, k)
                                if ( d == dd ) then !! δii = 1, δij = 0
                                    do ddd = 1, dim
                                        edot(d, d, i) = edot(d, d, i) &
                                            - 2._8 / 3                &
                                            * mass(j)/rho(j)          &
                                            * dv(ddd)                 &
                                            * dwdx(ddd, i, k)
                                    end do
                                end if !! d == dd
                            end do !! dd
                        end do !! d
                    end do !! k
                else !! Solid
                end if !! Fluid of Solid
            end do !! i
            !$OMP END PARALLEL DO
        end if !! Viscoisty
        

        !$OMP PARALLEL DO PRIVATE(i)
        do i = 1, ntotal !! All particles
            !!! Viscous entropy Tds/dt = 1/2 eta/rho
            if ( viscosity_w ) then
                tdsdt(i) = 0.5_8 &
                         * eta(i) & !! Dynamic viscosity
                         / rho(i) & !! Density
                         * sum(edot(:, :, i)**2)
            end if

            !!! Pressure from equation of state
            select case ( abs(itype(i)) )
            case (1)
                call gas_eos(rho(i), e(i), p(i), c(i))
            case (2)
                call arti_water_eos_1(rho(i), p(i), c(i))
            case (3)
                call arti_water_eos_2(rho(i), p(i), c(i))
            case (4)
                call tnt_eos(rho(i),e(i),p(i))
            case (5)
                call jwl_eos(rho(i), e(i), p(i))
            case (6)
                call mie_gruneisen_eos_of_water(rho(i), e(i), p(i))
            case (7)
                call water_polynomial_eos(rho(i), e(i), p(i))
            case (8)
                call mie_gruneisen_eos_of_solid(rho(i), e(i), p(i))
            end select !! abs(itype(i))
        end do !! i
        !$OMP END PARALLEL DO

        ! if ( nick == "undex_cylinder" ) call free_surface()

        !!! Calculate SPH sum for pressure force -p_a/rho
        !!! and viscous force eta_b/rho
        !!! and the internal energy change de/dt due to -p/rho*div_v
        !$OMP PARALLEL DO PRIVATE(i, j, k, d, rhoij, aux, Z_l, Z_r, e_ij, v_l, v_r, v_ij, p_star, v_star)
        do i = 1, ntotal  !! All particles
            if ( abs(itype(i)) < 8 ) then !! Fluid
                do k = 1, neighborNum(i) !! All neighbors of each particle
                    j = pair(i, k)

                    select case (pa_sph)

                    !!! For SPH algorithm 1
                    case (1)
                        !!! Auxiliary variables
                        rhoij = 1._8 / (rho(i)*rho(j))
                        aux = (p(i) + p(j)) * rhoij

                        !!! Conservation of Momentum
                        dvdt(:, i) = dvdt(: ,i)                  &
                            + mass(j)                            &
                            * ( - aux * dwdx(:, i, k)            &
                                + matmul(  eta(i)*edot(:, :, i)  &
                                         + eta(j)*edot(:, :, j), &
                                         dwdx(:, i, k)) * rhoij )

                        !!! Conservation of Energy
                        dedt(i) = dedt(i) &
                            + mass(j)     &
                            * aux         &
                            * dot_product(v(:, i)-v(:, j), dwdx(:, i, k))

                    case (2)
                        !!! Auxiliary variables
                        aux = (p(i)/rho(i)**2 + p(j)/rho(j)**2)

                        !!! Conservation of Momentum
                        dvdt(:, i) = dvdt(: ,i)                            &
                            + mass(j)                                      &
                            * ( - aux * dwdx(:, i, k)                      &
                                + matmul(  eta(i)*edot(:, :, i)/rho(i)**2  &
                                         + eta(j)*edot(:, :, j)/rho(j)**2, &
                                         dwdx(:, i, k)) )

                        !!! Conservation of Energy
                        dedt(i) = dedt(i) &
                            + mass(j)     &
                            * aux         &
                            * dot_product(v(:, i)-v(:, j), dwdx(:, i, k))
                    case default
                        write(err, "(1X, A, I0)") "Error SPH scheme ", pa_sph
                        error stop
                    end select !! pa_sph
                end do !! k
            end if
        end do !! i
        !$OMP END PARALLEL DO

        !!! Change of specific internal energy de/dt = Tds/dt - p/rho*div_v ?
        do i = 1, ntotal
            dedt(i) = tdsdt(i) + 0.5_8*dedt(i)
        end do

    end subroutine in_force

end module in_force_m