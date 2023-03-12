#include "../macro.h"
module in_force_m
    use ctrl_dict, only: dim, maxn, &
                         viscosity_w, pa_sph
    use parse_toml_m, only: nick
    ! use boundary_condition_m

    implicit none

contains
    subroutine in_force(ntotal, itype, x, v, mass, rho, p, e, c, &
                        neighborNum, neighborList, dwdx, dvdt, tdsdt, dedt, Shear, dSdt, Stress)
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
        integer, intent(in)  :: neighborNum(:)
        integer, intent(in)  :: neighborList(:, :)
        real(8), intent(in)  :: dwdx(:, :, :)
        real(8), intent(inout) :: dvdt(:, :), tdsdt(:), dedt(:)
        real(8), allocatable :: edot(:, :, :)
        real(8) :: dv(dim)
        real(8) :: rhoij, aux
        real(8) :: Z_l, Z_r, v_l, v_r, v_ij, e_ij(dim), v_star(dim), p_star
#if SOLID
        ! integer :: solid_num
        real(8), intent(in),    optional :: Shear(:, :, :)
        real(8), intent(inout), optional :: dSdt(:, :, :)
        real(8), intent(inout), optional :: Stress(:, :, :)
        real(8), allocatable :: rdot(:, :, :), aver_edot(:, :)
#endif
        integer i, j, k, d, dd, ddd

        !!! Initialization of Strain Rate Tensor, velocity divergence,
        !!! viscous energy, internal energy, acceleration
        allocate(edot(dim, dim, ntotal), source=0._8)

#if SOLID
        allocate(rdot(dim, dim, ntotal), source=0._8)
        allocate(aver_edot(dim, dim), source=0._8)
#endif

        !!! Dynamic viscosity
        if ( viscosity_w ) call viscosity(itype(1:ntotal), eta)


        if ( viscosity_w ) then
                !$OMP PARALLEL DO PRIVATE(i, j, k, d, dd, ddd, dv, rhoij, aux)
                do i = 1, ntotal !! All particles
                !!! Calculate SPH sum for Strain Rate Tensor of Fluid
                !!! εab = va,b + vb,a - 2/3*delta_ab*vc,c
#if SOLID
                if ( abs(itype(i)) < 8 ) then !! Fluid
#endif
                    do k = 1, neighborNum(i) !! All neighbors of each particle
                        j = neighborList(i, k)
                        dv = v(:, j) - v(:, i)
                        do d = 1, dim !! All dimensions For the First Order of Strain Rate Tensor, Loop 1
                            do dd = 1, dim !! All dimensions For the Second Order of Strain Rate Tensor, Loop 2
                                edot(d, dd, i) = edot(d, dd, i) &
                                    + mass(j)/rho(j)            &
                                    * dv(d) * dwdx(dd, i, k)
                                edot(dd, d, i) = edot(dd, d, i) &
                                    + mass(j)/rho(j)            &
                                    * dv(dd) * dwdx(d, i, k)
                                if ( d == dd ) then !! δii = 1, δij = 0
                                    do ddd = 1, dim
                                        edot(d, d, i) = edot(d, d, i) &
                                            - 2._8 / 3                &
                                            * mass(j)/rho(j)          &
                                            * dv(ddd) * dwdx(ddd, i, k)
                                    end do
                                end if !! d == dd
                            end do !! dd
                        end do !! d
                    end do !! k
#if SOLID
                !!! Calculate SPH sum for Strain Rate Tensor and Rotation Rate Tensor of Fluid
                !!! εab = 1/2 * (va,b + vb,a)
                !!! Rab = 1/2 * (va,b + vb,a)
                else if ( present(Shear) ) then !! Solid
                    do k = 1, neighborNum(i) !! All neighbors of each particle
                        j = neighborList(i, k)
                            dv = v(:, j) - v(:, i)
                            do d = 1, dim !! All dimensions For the First Order of Strain/Rotation Rate Tensor, Loop 1
                                do dd = 1, dim !! All dimensions For the Second Order of Strain/Rotation Rate Tensor, Loop 2
                                    associate (vab => mass(j)/rho(j) * dv(d) * dwdx(dd, i, k))
                                    associate (vba => mass(j)/rho(j) * dv(dd) * dwdx(d, i, k))
                                        edot(d, dd, i) = edot(d, dd, i) &
                                            + 0.5 * (vab + vba)
                                        rdot(d, dd, i) = rdot(d, dd, i) &
                                            + 0.5 * (vab - vba)
                                    end associate
                                    end associate
                                end do !! dd
                            end do !! d
                    end do !! k
            end if !! Fluid or Solid
#endif
            end do !! i
            !$OMP END PARALLEL DO
        end if !! Viscoisty


        !$OMP PARALLEL DO PRIVATE(i, aux, aver_edot, j, k, d, dd, ddd)
        do i = 1, ntotal !! All particles
#if SOLID
        if ( abs(itype(i)) < 8 ) then !! Fluid
#endif
            !!! Viscous entropy Tds/dt = 1/2 eta/rho εab•εab
            if ( viscosity_w ) then
                tdsdt(i) = 0.5_8 * eta(i) / rho(i) &
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
            end select !! abs(itype(i))

#if SOLID
        else if ( present(Shear) ) then !! Solid

            aux = 0
            do d = 1, dim
                aux = aux + edot(d, d, i)
            end do
            aux = aux / 3

            aver_edot = edot(:, :, i) !! Strain Rate Tensor of Solid
            do d = 1, dim
                aver_edot(d, d) = aver_edot(d, d) - aux
            end do

            tdsdt(i) = 1 / rho(i) * sum(Shear(:, :, i) * aver_edot(:, :))

            select case ( abs(itype(i)) )
            case (8)
                call mie_gruneisen_eos_of_solid(rho(i), e(i), p(i))
            end select

            do k = 1, neighborNum(i)
                j = neighborList(i, k)

                !!! Deviatoric Stress Rate Tensor
                do d = 1, dim !! All dimensions For the First Order of Deviatoric Stress Rate Tensor, Loop 1
                    do dd = 1, dim !! All dimensions For the Second Order of Deviatoric Stress Rate Tensor, Loop 2
                        dSdt(d, dd, i) = 2 * eta(i) * aver_edot(d, dd)
                        do ddd = 1, dim
                            dSdt(d, dd, i) = dSdt(d, dd, i) + &
                                Shear(d, ddd, i) * rdot(dd, ddd, i) &
                              + Shear(dd, ddd, i) * rdot(d, ddd, i)
                        end do !!! ddd
                    end do !! dd
                end do !! d
            end do !! k

            do d = 1, dim
                do dd = 1, dim
                    Stress(d, dd, i) = Shear(d, dd, i)
                end do
                Stress(d, d, i) = Stress(d, d, i) - p(i)
            end do

        end if !! Fluid or Solid
#endif
        end do !! i
        !$OMP END PARALLEL DO

        !!! Calculate SPH sum for pressure force -p_a/rho
        !!! and viscous force eta_b/rho
        !!! and the internal energy change de/dt due to -p/rho*div_v
        select case (pa_sph)

        !!! For SPH algorithm 1
        case (1)
            !$OMP PARALLEL DO PRIVATE(i, j, k, rhoij, aux)
            do i = 1, ntotal  !! All particles
                if ( abs(itype(i)) < 8 ) then !! Fluid
                    do k = 1, neighborNum(i) !! All neighbors of each particle
                        j = neighborList(i, k)

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
                        dedt(i) = dedt(i)   &
                            + mass(j) * aux &
                            * dot_product(v(:, i)-v(:, j), dwdx(:, i, k))

                    end do !! k
                else !! Solid
                end if !! Fluid or Solid
            end do !! i
            !$OMP END PARALLEL DO

        !!! For SPH algorithm 2
        case (2)
            !$OMP PARALLEL DO PRIVATE(i, j, k, aux)
            do i = 1, ntotal  !! All particles
                if ( abs(itype(i)) < 8 ) then !! Fluid
                    do k = 1, neighborNum(i) !! All neighbors of each particle
                        j = neighborList(i, k)

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
                        dedt(i) = dedt(i)   &
                            + mass(j) * aux &
                            * dot_product( v(:, i)-v(:, j), dwdx(:, i, k) )

                    end do !! k
                else if ( present(Shear) ) then!! Solid
                    do k = 1, neighborNum(i)
                        j = neighborList(i, k)

                        do d = 1, dim
                            do dd = 1, dim
                                !!! Conservation of Momentum
                                dvdt(d, i) = dvdt(d, i) &
                                    + mass(j) * ( Stress(d, dd, i) / rho(i)**2   &
                                                + Stress(d, dd, j) / rho(j)**2 ) &
                                              * dwdx(d, i, k)
                            end do
                            !!! Conservation of Energy
                            dedt(i) = dedt(i) &
                                + mass(j) * ( p(i) / rho(i)**2 + p(j) / rho(j)**2 ) &
                                * ( v(d, i)-v(d, j) ) * dwdx(d, i, k)
                        end do

                    end do !! k
                end if !! Fluid or Solid
            end do !! i
            !$OMP END PARALLEL DO

        !!! Riemann Solver
        case (3)
            !$OMP PARALLEL DO PRIVATE(i, j, k, rhoij, aux) &
            !$OMP PRIVATE(Z_l, Z_r, v_l, v_r, v_ij, v_star, p_star, e_ij)
            do i = 1, ntotal !! All particles
                if ( abs(itype(i)) < 8 ) then !! Fluid
                    do k = 1, neighborNum(i) !! All neighbors of each particle
                        j = neighborList(i, k)

                        !!! Auxiliary variables
                        rhoij = 1._8 / (rho(i)*rho(j))

                        Z_l = rho(i) * c(i)
                        Z_r = rho(j) * c(j)

                        e_ij = (x(:, j) - x(:, i)) / norm2(x(:, j) - x(:, i))
                        v_l = dot_product(v(:, i), e_ij)
                        v_r = dot_product(v(:, j), e_ij)

                        v_ij = ( Z_l*v_l + Z_r*v_r + (p(i)-p(j)) ) &
                             / ( Z_l + Z_r )

                        !!! Riemann Solution
                        p_star = ( Z_l*p(j) + Z_r*p(i) &
                            +   Z_l*Z_r*(v_l - v_r) )  &
                            / ( Z_l + Z_r )
                        v_star = v_ij * e_ij        &
                            + ( (v(:, i)+v(:, j))/2 &
                                - ((v_l+v_r)/2)*e_ij )

                        aux = 2 * p_star * rhoij

                        !!! Conservation of Momentum
                        dvdt(:, i) = dvdt(: ,i)                  &
                            + mass(j)                            &
                            * ( - aux * dwdx(:, i, k)            &
                                + matmul(  eta(i)*edot(:, :, i)  &
                                         + eta(j)*edot(:, :, j), &
                                         dwdx(:, i, k)) * rhoij )

                        !!! Conservation of Energy
                        dedt(i) = dedt(i)   &
                            + mass(j) * aux &
                            * dot_product(2 * (v(:, i)-v_star), dwdx(:, i, k))

                    end do
                else !! Solid
                end if !! Fluid or Solid
            end do
            !$OMP END PARALLEL DO
        case default
            write(err, "(1X, A, I0)") "Error SPH Scheme ", pa_sph
            error stop
        end select

        !!! Change of specific internal energy de/dt = Tds/dt - p/rho*div_v ?
        do i = 1, ntotal
            dedt(i) = tdsdt(i) + 0.5_8 * dedt(i)
        end do

        deallocate(edot)

#if SOLID
        deallocate(rdot, aver_edot)
#endif

    end subroutine in_force

end module in_force_m