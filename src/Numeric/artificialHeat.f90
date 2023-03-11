module arti_heat_m
    implicit none

contains
    subroutine arti_heat(ntotal, x, mass, rho, e, c, hsml, neighborNum, neighborList, dwdx, div_v, dedt)
        integer, intent(in) :: ntotal
        real(8), intent(in) :: x(:, :)
        real(8), intent(in) :: mass(:)
        real(8), intent(in) :: rho(:)
        real(8), intent(in) :: e(:)
        real(8), intent(in) :: c(:)
        real(8), intent(in) :: hsml(:)
        integer, intent(in) :: neighborNum(:)
        integer, intent(in) :: neighborList(:, :)
        real(8), intent(in) :: dwdx(:, :, :)
        real(8), intent(in)  :: div_v(:)
        real(8), intent(inout) :: dedt(:)
        real(8), parameter :: alpha = 0.1_8
        real(8), parameter :: beta  = 1.0_8
        real(8), parameter :: psi   = 0.1_8
        real(8) :: q_i, q_j, q_ij, hsml_ij, rho_ij, aux

        integer i, j, k

        !!! Calculate SPH sum for artificial heat conduction

        forall (i=1:ntotal)
            dedt(i)  = 0
        end forall

        !$OMP PARALLEL DO PRIVATE(i, j, k, q_i, q_j, q_ij, hsml_ij, rho_ij, aux)
        do i = 1, ntotal
            do k = 1, neighborNum(i)
                j = neighborList(i, k)

                q_i = alpha*hsml(i)*c(i) + beta*(hsml(i)**2)*(abs(div_v(i))-div_v(i))
                q_j = alpha*hsml(j)*c(j) + beta*(hsml(j)**2)*(abs(div_v(j))-div_v(j))  !! @todo: P375 ?
                ! q_i = (alpha*hsml(i)*rho(i)*c(i) + beta*(hsml(i)**2)*rho(i)*abs(div_v(i))) &  !! P123
                !       * abs(div_v(i))
                ! q_j = (alpha*hsml(j)*rho(j)*c(j) + beta*(hsml(j)**2)*rho(j)*abs(div_v(j))) &
                !       * abs(div_v(j))
                q_ij    = 0.5_8 * (q_i + q_j)
                hsml_ij = 0.5_8 * (hsml(i) + hsml(j))
                rho_ij  = 0.5_8 * (rho(i)  + rho(j) )

                associate (dx => x(:, i) - x(:, j))
                    aux = q_ij * sum(dx*dwdx(:, i, k)) &
                        / (rho_ij * (sum(dx**2) + 0.01_8*hsml_ij**2))
                end associate

                dedt(i) = dedt(i) + 2._8 * mass(j) * aux * (e(i)-e(j))

            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine arti_heat

end module arti_heat_m