module arti_visc_m
    use ctrl_dict, only: dim
    use parse_toml_m, only: nick
    implicit none

contains
    subroutine arti_visc(ntotal, x, v, mass, rho, c, hsml, pair, dwdx, neighborNum, dvdt, dedt)
        integer, intent(in)  :: ntotal
        real(8), intent(in) :: x(:, :)
        real(8), intent(in) :: v(:, :)
        real(8), intent(in) :: mass(:)
        real(8), intent(in) :: rho(:)
        real(8), intent(in) :: c(:)
        real(8), intent(in) :: hsml(:)
        integer, intent(in) :: pair(:, :)
        real(8), intent(in) :: dwdx(:, :, :)
        integer, intent(in) :: neighborNum(:)
        real(8), intent(inout) :: dvdt(:,:)
        real(8), intent(inout) :: dedt(:)
        real(8) :: alpha !! Bulk viscosity
        real(8) :: beta  !! Shear viscosity 
                         !! To Avoid non-physical penetration
                         !! Shall be 10 in HE simulation
        real(8) :: psi   !! Parameter to avoid singularities
        real(8) :: dx(dim), dv(dim)
        real(8) :: xv, hsml_ij, rho_ij, c_ij, phi_ij, PI_ij

        integer i, j, k

        forall (i=1:ntotal)
            dvdt(:, i) = 0
            dedt(i)    = 0
        end forall
        

        select case(nick)
        ! case ("dam_break")
        !     alpha = 1
        !     beta  = 1
        !     psi   = 0.1
        ! case("shear_cavity","shock_tube")
        !     alpha = 1
        !     beta  = 1
        !     psi   = 0.1
        case("tnt_bar", "tnt_cylinder", "undex_cylinder", "undex_chamber", "armco_iron_collide")
            alpha = 1
            beta  = 10
            psi   = 0.1
        case default
            alpha = 1
            beta  = 1
            psi   = 0.1
        end select
        
        !!! Calculate SPH sum for artificial viscous
        do i = 1, ntotal
            do k = 1, neighborNum(i)
                j = pair(i, k)
                
                dx = x(:, i) - x(:, j)
                dv = v(:, i) - v(:, j)
                xv = sum( dx * dv )

                !!! Artificial viscous force only if v_ij * r_ij < 0
                if (xv < 0._8) then
                    hsml_ij = 0.5_8 * ( hsml(i) + hsml(j) )
                    rho_ij  = 0.5_8 * ( rho(i)  + rho(j)  )
                    c_ij    = 0.5_8 * ( c(i)    + c(j)    )

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
                    associate (aux => -PI_ij*dwdx(:, i, k))
                        dvdt(:, i) = dvdt(:, i) + mass(j) * aux
                        dedt(i)    = dedt(i) - 0.5_8 * mass(j) * sum(dv*aux)
                    end associate
                end if
            end do
        end do

    end subroutine arti_visc

end module arti_visc_m