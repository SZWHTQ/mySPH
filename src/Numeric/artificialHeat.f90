module arti_heat_m
    use sph
    implicit none

contains
    subroutine arti_heat(P, dedt)
        type(Particle), intent(in) :: P(:)
        real(8), intent(inout) :: dedt(:)
        integer :: ntotal
        real(8), parameter :: alpha = 0.1_8
        real(8), parameter :: beta  = 1.0_8
        real(8), parameter :: psi   = 0.1_8
        real(8) :: q_i, q_j, q_ij, hsml_ij, rho_ij, aux

        integer i, j, k

        !!! Calculate SPH sum for artificial heat conduction
        ntotal = size(P)
        forall (i=1:ntotal) dedt(i)  = 0

        !$OMP PARALLEL DO PRIVATE(i, j, k, q_i, q_j, q_ij, hsml_ij, rho_ij, aux)
        do i = 1, ntotal
            do k = 1, P(i)%neighborNum
                j = P(i)%neighborList(k)

                q_i = alpha * P(i)%SmoothingLength*P(i)%SoundSpeed &
                    + beta * (P(i)%SmoothingLength**2)             &
                           * (abs(P(i)%divergenceVelocity)-P(i)%divergenceVelocity)
                q_j = alpha * P(j)%SmoothingLength*P(j)%SoundSpeed &
                    + beta * (P(j)%SmoothingLength**2)             &
                           * (abs(P(j)%divergenceVelocity)-P(j)%divergenceVelocity)  !! @todo: P375 ?
                ! q_i = (alpha*P(i)%SmoothingLength*P(i)%Density*P(i)%SoundSpeed + beta*(P(i)%SmoothingLength**2)*P(i)%Density*abs(P(i)%divergenceVelocity)) &  !! P123
                !       * abs(P(i)%divergenceVelocity)
                ! q_j = (alpha*P(j)%SmoothingLength*P(j)%Density*P(j)%SoundSpeed + beta*(P(j)%SmoothingLength**2)*P(j)%Density*abs(P(j)%divergenceVelocity)) &
                !       * abs(P(j)%divergenceVelocity)
                q_ij    = 0.5_8 * (q_i + q_j)
                hsml_ij = 0.5_8 * (P(i)%SmoothingLength + P(j)%SmoothingLength)
                rho_ij  = 0.5_8 * (P(i)%Density  + P(j)%Density )

                associate (dx => P(i)%x(:) - P(j)%x(:))
                    aux = q_ij * sum(dx*P(i)%dwdx(:, k)) &
                        / (rho_ij * (sum(dx**2) + 0.01_8*hsml_ij**2))
                end associate

                dedt(i) = dedt(i)            &
                    + 2._8 * P(j)%Mass * aux &
                        * (P(i)%InternalEnergy-P(j)%InternalEnergy)

            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine arti_heat

end module arti_heat_m