module arti_visc_m
    use ctrl_dict, only: Field, Project
    use sph,       only: Particle
    implicit none

contains
    subroutine arti_visc(ntotal, P, dvdt, dedt)
        integer, intent(in) :: ntotal
        type(Particle), intent(in) :: P(:)
        real(8), intent(inout) :: dvdt(:,:)
        real(8), intent(inout) :: dedt(:)
        real(8) :: alpha !! Bulk viscosity
        real(8) :: beta  !! Shear viscosity 
                         !! To Avoid non-physical penetration
                         !! Shall be 10 in HE simulation
        real(8) :: psi   !! Parameter to avoid singularities
        real(8) :: dx(Field%Dim), dv(Field%Dim)
        real(8) :: xv, hsml_ij, rho_ij, c_ij, phi_ij, PI_ij

        integer i, j, k

        do i=1, ntotal
            dvdt(:, i) = 0
            dedt(i)    = 0
        end do
        

        select case(Project%nick)
        ! case ("dam_break")
        !     alpha = 1
        !     beta  = 1
        !     psi   = 0.1
        ! case("shear_cavity","shock_tube")
        !     alpha = 1
        !     beta  = 1
        !     psi   = 0.1
        case("tnt_bar", "tnt_cylinder", "undex_cylinder", "undex_chamber")
            alpha = 1
            beta  = 10
            psi   = 0.1
        case("taylor_rod", "can_beam")
            alpha = 0.5
            beta  = 0.5
            psi   = 0.1
        case default
            alpha = 1
            beta  = 1
            psi   = 0.1
        end select
        
        !!! Calculate SPH sum for artificial viscous
        !$OMP PARALLEL DO PRIVATE(i, j, k, dx, dv, xv, hsml_ij, rho_ij, c_ij, phi_ij, PI_ij)
        do i = 1, ntotal
            do k = 1, P(i)%neighborNum
                j = P(i)%neighborList(k)
                
                dx = P(i)%x(:) - P(j)%x(:)
                dv = P(i)%v(:) - P(j)%v(:)
                xv = sum( dx * dv )

                !!! Artificial viscous force only if v_ij * r_ij < 0
                if (xv < 0._8) then
                    hsml_ij = 0.5_8 * ( P(i)%SmoothingLength + P(j)%SmoothingLength )
                    rho_ij  = 0.5_8 * ( P(i)%Density  + P(j)%Density  )
                    c_ij    = 0.5_8 * ( P(i)%SoundSpeed    + P(j)%SoundSpeed    )

                    !!! Calculate phi_ij = hsml*v_ij*r_ij / (r_ij^2+psi^2*hsml_ij^2)
                    !!! varphi = psi * hsml_ij
                    phi_ij = hsml_ij        &
                           * xv             &
                           / ( sum( dx**2 ) + (psi * hsml_ij)**2 )

                    !!! Calculate PI_ij = (-alpha*phi_ij*c_ij + beta*phi_ij^2) / rho_ij
                    PI_ij = (beta * phi_ij - alpha * c_ij) &
                          * phi_ij                         &
                          / rho_ij

                    !!! Calculate SPH sum for artificial viscous force
                    associate (aux => PI_ij*P(i)%dwdx(:, k))
                        dvdt(:, i) = dvdt(:, i) - P(j)%mass * aux
                        dedt(i)    = dedt(i) + 0.5_8 * P(j)%mass * sum(dv*aux)
                    end associate
                end if
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine arti_visc

end module arti_visc_m