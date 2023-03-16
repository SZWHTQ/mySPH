#include "../macro.h"
module in_force_m
    use ctrl_dict, only: Config, Field
    use sph,       only: Particle
    ! use boundary_condition_m

    implicit none

contains
    subroutine in_force(P, dvdt, dedt, Shear, dSdt)
        use, intrinsic :: iso_fortran_env, only: err => error_unit
        use visc_m
        use eos_m
        !$ use omp_lib
        type(Particle), intent(inout) :: P(:)
        real(8), intent(inout) :: dvdt(:, :), dedt(:)
        integer :: ntotal
        real(8), allocatable :: tdsdt(:), edot(:, :, :)
        real(8) :: dv(Field%dim)
        real(8) :: rhoij, aux
        real(8) :: Z_l, Z_r, v_l, v_r, v_ij, e_ij(Field%dim), v_star(Field%dim), p_star
#if SOLID
        ! integer :: solid_num
        real(8), intent(in),    optional :: Shear(:, :, :)
        real(8), intent(inout), optional :: dSdt(:, :, :)
        real(8), allocatable :: rdot(:, :, :), aver_edot(:, :)
#endif
        integer i, j, k, d, dd, ddd

        ntotal = size(P)
        !!! Initialization of Strain Rate Tensor, velocity divergence,
        !!! viscous energy, internal energy, acceleration
        allocate(tdsdt(ntotal), source=0._8)
        allocate(edot(Field%dim, Field%dim, ntotal), source=0._8)

#if SOLID
        allocate(rdot(Field%dim, Field%dim, ntotal), source=0._8)
        allocate(aver_edot(Field%dim, Field%dim), source=0._8)
#endif

        !!! Dynamic viscosity
        if ( Config%viscosity_w ) call viscosity(P(1:ntotal)%Type, P(1:ntotal)%KineticViscocity)


        if ( Config%viscosity_w ) then
                !$OMP PARALLEL DO PRIVATE(i, j, k, d, dd, ddd, dv, rhoij, aux)
                do i = 1, ntotal !! All particles
                !!! Calculate SPH sum for Strain Rate Tensor of Fluid
                !!! εab = va,b + vb,a - 2/3*delta_ab*vc,c
#if SOLID
                if ( abs(P(i)%Type) < 8 ) then !! Fluid
#endif
                    do k = 1, P(i)%neighborNum !! All neighbors of each particle
                        j = P(i)%neighborList(k)
                        dv = P(j)%v(:) - P(i)%v(:)
                        do d = 1, Field%dim !! All dimensions For the First Order of Strain Rate Tensor, Loop 1
                            do dd = 1, Field%dim !! All dimensions For the Second Order of Strain Rate Tensor, Loop 2
                                edot(d, dd, i) = edot(d, dd, i) &
                                    + P(j)%Mass/P(j)%Density            &
                                    * dv(d) * P(i)%dwdx(dd, k)
                                edot(dd, d, i) = edot(dd, d, i) &
                                    + P(j)%Mass/P(j)%Density            &
                                    * dv(dd) * P(i)%dwdx(d, k)
                                if ( d == dd ) then !! δii = 1, δij = 0
                                    do ddd = 1, Field%dim
                                        edot(d, d, i) = edot(d, d, i) &
                                            - 2._8 / 3                &
                                            * P(j)%Mass/P(j)%Density          &
                                            * dv(ddd) * P(i)%dwdx(ddd, k)
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
                    do k = 1, P(i)%neighborNum !! All neighbors of each particle
                        j = P(i)%neighborList(k)
                            dv = P(j)%v(:) - P(i)%v(:)
                            do d = 1, Field%dim !! All dimensions For the First Order of Strain/Rotation Rate Tensor, Loop 1
                                do dd = 1, Field%dim !! All dimensions For the Second Order of Strain/Rotation Rate Tensor, Loop 2
                                    associate (vab => &
                                        P(j)%Mass/P(j)%Density * dv(d) * P(i)%dwdx(dd, k))
                                    associate (vba => &
                                        P(j)%Mass/P(j)%Density * dv(dd) * P(i)%dwdx(d, k))
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
        if ( abs(P(i)%Type) < 8 ) then !! Fluid
#endif
            !!! Viscous entropy Tds/dt = 1/2 eta/rho εab•εab
            if ( Config%viscosity_w ) then
                tdsdt(i) = 0.5_8 * P(i)%KineticViscocity / P(i)%Density &
                         * sum(edot(:, :, i)**2)
            end if

            !!! Pressure from equation of state
            select case ( abs(P(i)%Type) )
            case (1)
                call gas_eos(P(i)%Density, P(i)%InternalEnergy, P(i)%Pressure, P(i)%SoundSpeed)
            case (2)
                call arti_water_eos_1(P(i)%Density, P(i)%Pressure, P(i)%SoundSpeed)
            case (3)
                call arti_water_eos_2(P(i)%Density, P(i)%Pressure, P(i)%SoundSpeed)
            case (4)
                call tnt_eos(P(i)%Density,P(i)%InternalEnergy,P(i)%Pressure)
            case (5)
                call jwl_eos(P(i)%Density, P(i)%InternalEnergy, P(i)%Pressure)
            case (6)
                call mie_gruneisen_eos_of_water(P(i)%Density, P(i)%InternalEnergy, P(i)%Pressure)
            case (7)
                call water_polynomial_eos(P(i)%Density, P(i)%InternalEnergy, P(i)%Pressure)
            end select !! abs(P(i)%Type)

#if SOLID
        else if ( present(Shear) ) then !! Solid

            aux = 0
            do d = 1, Field%dim
                aux = aux + edot(d, d, i)
            end do
            aux = aux / 3

            aver_edot = edot(:, :, i) !! Strain Rate Tensor of Solid
            do d = 1, Field%dim
                aver_edot(d, d) = aver_edot(d, d) - aux
            end do

            tdsdt(i) = 1 / P(i)%Density * sum(Shear(:, :, i) * aver_edot(:, :))

            select case ( abs(P(i)%Type) )
            case (8)
                call mie_gruneisen_eos_of_solid(P(i)%Density, P(i)%InternalEnergy, P(i)%Pressure)
            end select

            do k = 1, P(i)%neighborNum
                j = P(i)%neighborList(k)

                !!! Deviatoric Stress Rate Tensor
                do d = 1, Field%dim !! All dimensions For the First Order of Deviatoric Stress Rate Tensor, Loop 1
                    do dd = 1, Field%dim !! All dimensions For the Second Order of Deviatoric Stress Rate Tensor, Loop 2
                        dSdt(d, dd, i) = 2 * P(i)%KineticViscocity * aver_edot(d, dd)
                        do ddd = 1, Field%dim
                            dSdt(d, dd, i) = dSdt(d, dd, i) + &
                                Shear(d, ddd, i) * rdot(dd, ddd, i) &
                              + Shear(dd, ddd, i) * rdot(d, ddd, i)
                        end do !!! ddd
                    end do !! dd
                end do !! d
            end do !! k

            do d = 1, Field%dim
                do dd = 1, Field%dim
                    P(i)%Stress(d, dd) = Shear(d, dd, i)
                end do
                P(i)%Stress(d, d) = P(i)%Stress(d, d) - P(i)%Pressure
            end do

        end if !! Fluid or Solid
#endif
        end do !! i
        !$OMP END PARALLEL DO

        !!! Calculate SPH sum for pressure force -p_a/rho
        !!! and viscous force eta_b/rho
        !!! and the internal energy change de/dt due to -p/rho*div_v
        select case (Config%pa_sph)

        !!! For SPH algorithm 1
        case (1)
            !$OMP PARALLEL DO PRIVATE(i, j, k, rhoij, aux)
            do i = 1, ntotal  !! All particles
                if ( abs(P(i)%Type) < 8 ) then !! Fluid
                    do k = 1, P(i)%neighborNum !! All neighbors of each particle
                        j = P(i)%neighborList(k)

                        !!! Auxiliary variables
                        rhoij = 1._8 / (P(i)%Density*P(j)%Density)
                        aux = (P(i)%Pressure + P(j)%Pressure) * rhoij

                        !!! Conservation of Momentum
                        dvdt(:, i) = dvdt(: ,i)                  &
                            + P(j)%Mass                            &
                            * ( - aux * P(i)%dwdx(:, k)            &
                                + matmul(  P(i)%KineticViscocity*edot(:, :, i)  &
                                         + P(j)%KineticViscocity*edot(:, :, j), &
                                         P(i)%dwdx(:, k)) * rhoij )

                        !!! Conservation of Energy
                        dedt(i) = dedt(i)   &
                            + P(j)%Mass * aux &
                            * dot_product(P(i)%v(:)-P(j)%v(:), P(i)%dwdx(:, k))

                    end do !! k
                else !! Solid
                end if !! Fluid or Solid
            end do !! i
            !$OMP END PARALLEL DO

        !!! For SPH algorithm 2
        case (2)
            !$OMP PARALLEL DO PRIVATE(i, j, k, aux)
            do i = 1, ntotal  !! All particles
                if ( abs(P(i)%Type) < 8 ) then !! Fluid
                    do k = 1, P(i)%neighborNum !! All neighbors of each particle
                        j = P(i)%neighborList(k)

                        !!! Auxiliary variables
                        aux = (P(i)%Pressure/P(i)%Density**2 + P(j)%Pressure/P(j)%Density**2)

                        !!! Conservation of Momentum
                        dvdt(:, i) = dvdt(: ,i)                                                 &
                            + P(j)%Mass                                                         &
                            * ( - aux * P(i)%dwdx(:, k)                                         &
                                + matmul(  P(i)%KineticViscocity*edot(:, :, i)/P(i)%Density**2  &
                                         + P(j)%KineticViscocity*edot(:, :, j)/P(j)%Density**2, &
                                         P(i)%dwdx(:, k)) )

                        !!! Conservation of Energy
                        dedt(i) = dedt(i)   &
                            + P(j)%Mass * aux &
                            * dot_product( P(i)%v(:)-P(j)%v(:), P(i)%dwdx(:, k) )

                    end do !! k
                else if ( present(Shear) ) then!! Solid
                    do k = 1, P(i)%neighborNum
                        j = P(i)%neighborList(k)

                        do d = 1, Field%dim
                            do dd = 1, Field%dim
                                !!! Conservation of Momentum
                                dvdt(d, i) = dvdt(d, i) &
                                    + P(j)%Mass * ( P(i)%Stress(d, dd) / P(i)%Density**2   &
                                                  + P(j)%Stress(d, dd) / P(j)%Density**2 ) &
                                              * P(i)%dwdx(d, k)
                            end do
                            !!! Conservation of Energy
                            dedt(i) = dedt(i) &
                                + P(j)%Mass * ( P(i)%Pressure / P(i)%Density**2   &
                                              + P(j)%Pressure / P(j)%Density**2 ) &
                                * ( P(i)%v(d)-P(j)%v(d) ) * P(i)%dwdx(d, k)
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
                if ( abs(P(i)%Type) < 8 ) then !! Fluid
                    do k = 1, P(i)%neighborNum !! All neighbors of each particle
                        j = P(i)%neighborList(k)

                        !!! Auxiliary variables
                        rhoij = 1._8 / (P(i)%Density*P(j)%Density)

                        Z_l = P(i)%Density * P(i)%SoundSpeed
                        Z_r = P(j)%Density * P(j)%SoundSpeed

                        e_ij = (P(j)%x(:) - P(i)%x(:)) / norm2(P(j)%x(:) - P(i)%x(:))
                        v_l = dot_product(P(i)%v(:), e_ij)
                        v_r = dot_product(P(j)%v(:), e_ij)

                        v_ij = ( Z_l*v_l + Z_r*v_r + (P(i)%Pressure-P(j)%Pressure) ) &
                             / ( Z_l + Z_r )

                        !!! Riemann Solution
                        p_star = ( Z_l*P(j)%Pressure + Z_r*P(i)%Pressure &
                            +   Z_l*Z_r*(v_l - v_r) )  &
                            / ( Z_l + Z_r )
                        v_star = v_ij * e_ij        &
                            + ( (P(i)%v(:)+P(j)%v(:))/2 &
                                - ((v_l+v_r)/2)*e_ij )

                        aux = 2 * p_star * rhoij

                        !!! Conservation of Momentum
                        dvdt(:, i) = dvdt(: ,i)                                 &
                            + P(j)%Mass                                         &
                            * ( - aux * P(i)%dwdx(:, k)                         &
                                + matmul(  P(i)%KineticViscocity*edot(:, :, i)  &
                                         + P(j)%KineticViscocity*edot(:, :, j), &
                                         P(i)%dwdx(:, k)) * rhoij )

                        !!! Conservation of Energy
                        dedt(i) = dedt(i)   &
                            + P(j)%Mass * aux &
                            * dot_product(2 * (P(i)%v(:)-v_star), P(i)%dwdx(:, k))

                    end do
                else !! Solid
                end if !! Fluid or Solid
            end do
            !$OMP END PARALLEL DO
        case default
            write(err, "(1X, A, I0)") "Error SPH Scheme ", Config%pa_sph
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