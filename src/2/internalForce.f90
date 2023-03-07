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
        real(8) :: v_aux(dim, dim), e_aux, rhoij, aux
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
        end if !! Viscoisty

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

        ! if ( nick == "undex_cylinder" ) call free_surface()

        !!! Calculate SPH sum for pressure force -p_a/rho
        !!! and viscous force eta_b/rho
        !!! and the internal energy change de/dt due to -p/rho*div_v
        ! do k = 1, niac
        !     i = pair(k, 1)
        !     j = pair(k, 2)
        !     e_aux = 0

        !     !!! For SPH algorithm 1
        !     select case (pa_sph)
        !     case (1)
        !         rhoij = 1._8 / (rho(i)*rho(j))
        !         do d = 1, dim
        !             !!! Pressure part
        !             v_aux = -(p(i) + p(j)) * dwdx(d, k)
        !             e_aux = e_aux + (v(d, j) - v(d, i)) * v_aux
        !             !!! Viscous force
        !             if ( viscosity_w ) then
        !                 select case (d)
        !                 case (1)
        !                     !!! x-coordinate of acceleration
        !                     v_aux = v_aux &
        !                     + (eta(i)*txx(i) + eta(j)*txx(j))*dwdx(1, k)
        !                     if ( dim >= 2 ) then
        !                         v_aux = v_aux &
        !                         + (eta(i)*txy(i) + eta(j)*txy(j))*dwdx(2, k)
        !                         if ( dim == 3 ) then
        !                             v_aux = v_aux &
        !                             + (eta(i)*txz(i) + eta(j)*txz(j))*dwdx(3, k)
        !                         end if
        !                     end if

        !                 case (2)
        !                     !!! y-coordinate of acceleration
        !                     v_aux = v_aux &
        !                     + (eta(i)*txx(i) + eta(j)*txy(j))*dwdx(1, k) &
        !                     + (eta(i)*tyy(i) + eta(j)*tyy(j))*dwdx(2, k)
        !                     if ( dim == 3 ) then
        !                         v_aux = v_aux &
        !                         + (eta(i)*tyz(i) + eta(j)*tyz(j))*dwdx(3, k)
        !                     end if

        !                 case (3)
        !                     !!! z-coordinate of acceleration
        !                     v_aux = v_aux &
        !                     + (eta(i)*txz(i) + eta(j)*txz(j))*dwdx(1, k) &
        !                     + (eta(i)*tyz(i) + eta(j)*tyz(j))*dwdx(2, k) &
        !                     + (eta(i)*tzz(i) + eta(j)*tzz(j))*dwdx(3, k)
        !                 end select
        !             end if

        !             v_aux = v_aux * rhoij
        !             dvdt(d, i) = dvdt(d, i) + mass(j) * v_aux
        !             dvdt(d, j) = dvdt(d, j) - mass(i) * v_aux
        !         end do

        !         e_aux = e_aux * rhoij
        !         dedt(i) = dedt(i) + mass(j) * e_aux
        !         dedt(j) = dedt(j) + mass(i) * e_aux

        !     !!!  For SPH algorithm 2
        !     case (2)
        !         do d = 1, dim
        !             !!! Pressure part
        !             v_aux = -(p(i)/rho(i)**2 + p(j)/rho(j)**2) * dwdx(d, k)
        !             e_aux = e_aux + (v(d, j) - v(d, i)) * v_aux

        !             !!! Viscous force
        !             if ( viscosity_w ) then
        !                 select case (d)
        !                 case (1)
        !                     !!! x-coordinate of acceleration
        !                     v_aux = v_aux &
        !                     + (eta(i)*txx(i)/rho(i)**2 + eta(j)*txx(j)/rho(j)**2)*dwdx(1, k)
        !                     if ( dim >= 2 ) then
        !                         v_aux = v_aux &
        !                         + (eta(i)*txy(i)/rho(i)**2 + eta(j)*txy(j)/rho(j)**2)*dwdx(2, k)
        !                         if ( dim == 3 ) then
        !                             v_aux = v_aux &
        !                             + (eta(i)*txz(i)/rho(i)**2 + eta(j)*txz(j)/rho(j)**2)*dwdx(3, k)
        !                         end if
        !                     end if

        !                 case (2)
        !                     !!! y-coordinate of acceleration
        !                     v_aux = v_aux &
        !                     + (eta(i)*txy(i)/rho(i)**2 + eta(j)*txy(j)/rho(j)**2)*dwdx(1, k) &
        !                     + (eta(i)*tyy(i)/rho(i)**2 + eta(j)*tyy(j)/rho(j)**2)*dwdx(2, k)
        !                     if ( dim == 3 ) then
        !                         v_aux = v_aux &
        !                         + (eta(i)*tyz(i)/rho(i)**2 + eta(j)*tyz(j)/rho(j)**2)*dwdx(3, k)
        !                     end if

        !                 case (3)
        !                     !!! z-coordinate of acceleration
        !                     v_aux = v_aux &
        !                     + (eta(i)*txz(i)/rho(i)**2 + eta(j)*txz(j)/rho(j)**2)*dwdx(1, k) &
        !                     + (eta(i)*tyz(i)/rho(i)**2 + eta(j)*tyz(j)/rho(j)**2)*dwdx(2, k) &
        !                     + (eta(i)*tzz(i)/rho(i)**2 + eta(j)*tzz(j)/rho(j)**2)*dwdx(3, k)

        !                 end select
        !             end if

        !             dvdt(d, i) = dvdt(d, i) + mass(j) * v_aux
        !             dvdt(d, j) = dvdt(d, j) - mass(i) * v_aux
        !         end do

        !         dedt(i) = dedt(i) + mass(j) * e_aux
        !         dedt(j) = dedt(j) + mass(i) * e_aux

        !     !!! Riemann solver
        !     case (3)
        !         rhoij = 1._8 / (rho(i)*rho(j))

        !         Z_l = rho(i) * c(i)
        !         Z_r = rho(j) * c(j)

        !         e_ij = (x(:, j) - x(:, i)) / norm2(x(:, j) - x(:, i))
        !         v_l = dot_product(v(:, i), e_ij)
        !         v_r = dot_product(v(:, j), e_ij)

        !         v_ij = ( Z_l*v_l + Z_r*v_r + (p(i)-p(j)) ) &
        !              / ( Z_l + Z_r )
        !         p_star = ( Z_l*p(j) + Z_r*p(i)   &
        !                +   Z_l*Z_r*(v_l - v_r) ) &
        !                / ( Z_l + Z_r )
        !         v_star = v_ij*e_ij + ((v(:, i)+v(:, j))/2 - ((v_l+v_r)/2)*e_ij)

        !         !!! For particle i
        !         do d = 1, dim
        !             !!! Pressure part
        !             v_aux = -2 * p_star * dwdx(d, k)
        !             e_aux = e_aux - (v(d, i) - v_star(d)) * v_aux
        !             !!! Viscous force
        !             if ( viscosity_w ) then
        !                 select case (d)
        !                 case (1)
        !                     !!! x-coordinate of acceleration
        !                     v_aux = v_aux &
        !                     + (eta(i)*txx(i) + eta(j)*txx(j))*dwdx(1, k)
        !                     if ( dim >= 2 ) then
        !                         v_aux = v_aux &
        !                         + (eta(i)*txy(i) + eta(j)*txy(j))*dwdx(2, k)
        !                         if ( dim == 3 ) then
        !                             v_aux = v_aux &
        !                             + (eta(i)*txz(i) + eta(j)*txz(j))*dwdx(3, k)
        !                         end if
        !                     end if

        !                 case (2)
        !                     !!! y-coordinate of acceleration
        !                     v_aux = v_aux &
        !                     + (eta(i)*txx(i) + eta(j)*txy(j))*dwdx(1, k) &
        !                     + (eta(i)*tyy(i) + eta(j)*tyy(j))*dwdx(2, k)
        !                     if ( dim == 3 ) then
        !                         v_aux = v_aux &
        !                         + (eta(i)*tyz(i) + eta(j)*tyz(j))*dwdx(3, k)
        !                     end if

        !                 case (3)
        !                     !!! z-coordinate of acceleration
        !                     v_aux = v_aux &
        !                     + (eta(i)*txz(i) + eta(j)*txz(j))*dwdx(1, k) &
        !                     + (eta(i)*tyz(i) + eta(j)*tyz(j))*dwdx(2, k) &
        !                     + (eta(i)*tzz(i) + eta(j)*tzz(j))*dwdx(3, k)
        !                 end select
        !             end if

        !             v_aux = v_aux * rhoij
        !             dvdt(d, i) = dvdt(d, i) + mass(j) * v_aux
        !             dvdt(d, j) = dvdt(d, j) - mass(i) * v_aux

        !         end do

        !         e_aux = e_aux * rhoij
        !         dedt(i) = dedt(i) + mass(j) * e_aux

        !         !!! For particle j
        !         e_aux = sum(v(:, j) - v_star(:)) * v_aux
        !         dedt(j) = dedt(j) + mass(i) * e_aux

        !     case default
        !         write(err, "(1X, A, I0)") "Error SPH scheme ", pa_sph
        !         error stop
        !     end select
        ! end do

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
                            * dot_product(v(:, i), dwdx(:, i, k))

                    case (2)
                        !!! Auxiliary variables
                        aux = (p(i)/rho(i)**2 + p(j)/rho(j)**2)

                        !!! Conservation of Momentum
                        dvdt(:, i) = dvdt(: ,i)                              &
                            + mass(j)                                        &
                            * ( - aux * dwdx(:, i, k)                        &
                                + matmul(  eta(i)*edot(:, :, i) / rho(i)**2  &
                                         + eta(j)*edot(:, :, j) / rho(j)**2, &
                                         dwdx(:, i, k)) )

                        !!! Conservation of Energy
                        dedt(i) = dedt(i) &
                            + mass(j)     &
                            * aux         &
                            * dot_product(v(:, i), dwdx(:, i, k))
                    end select !! pa_sph
                end do !! k

            end if
        end do !! i

        !!! Change of specific internal energy de/dt = Tds/dt - p/rho*div_v ?
        do i = 1, ntotal
            dedt(i) = tdsdt(i) + 0.5_8*dedt(i)
        end do

    end subroutine in_force

end module in_force_m