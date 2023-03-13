module corr_velo_m
    use ctrl_dict, only: dim
    use sph
    use parse_toml_m, only: nick
    implicit none

contains
    subroutine aver_velo(ntotal, P, aver_v)
        integer, intent(in) :: ntotal
        type(Particle), intent(in) :: P(:)
        real(8), intent(inout) :: aver_v(:, :)
        real(8) :: epsilon
        real(8) :: dv(dim)

        integer i, j, k

        !!! ε (epsilon), a small constant chosen by experience,
        !!! may lead to instability.
        !!! for example, for the 1 dimensional shock tube problem,
        !!! the ε ≤ 0.3
        parameter(epsilon = 0.3)

        forall(i=1:dim, j=1:ntotal) aver_v(i, j) = 0._8

        do i = 1, ntotal
            do k = 1, P(i)%neighborNum
                j = P(i)%neighborList(k)
                
                dv = P(i)%v(:) - P(j)%v(:)
                aver_v(:, i) = aver_v(:, i) &
                    - epsilon               &
                    * 2 * P(j)%mass * dv      &
                    / (P(i)%Density + P(j)%Density)     &
                    * P(i)%w(k)
            end do
        end do

    end subroutine aver_velo

end module corr_velo_m