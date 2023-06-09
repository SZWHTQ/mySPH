module divergence_m
    use sph, only: Particle
    implicit none

contains
    subroutine divergenceVelocity(ntotal, P)
        integer, intent(in) :: ntotal
        type(Particle), intent(inout) :: P(:)

        integer i, j, k

        do i=1, ntotal
            P(i)%divergenceVelocity = 0
        end do

        !$OMP PARALLEL DO PRIVATE(i, j, k)
        do i = 1, ntotal
            do k = 1, P(i)%neighborNum
                j = P(i)%neighborList(k)
                P(i)%divergenceVelocity = P(i)%divergenceVelocity &
                    + P(j)%Mass/P(j)%Density * sum((P(j)%v(:)-P(i)%v(:)) * P(i)%dwdx(:, k))
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine divergenceVelocity

    subroutine divergencePosition(ntotal, P)
        integer, intent(in) :: ntotal
        type(Particle), intent(inout) :: P(:)

        integer i, j, k

        do i=1, ntotal
            P(i)%divergencePosition = 0
        end do

        !$OMP PARALLEL DO PRIVATE(i, j, k)
        do i = 1, ntotal
            do k = 1, P(i)%neighborNum
                j = P(i)%neighborList(k)
                P(i)%divergencePosition = P(i)%divergencePosition &
                    + P(j)%Mass/P(j)%Density * sum((P(j)%x(:)-P(i)%x(:)) * P(i)%dwdx(:, k))
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine divergencePosition

end module divergence_m