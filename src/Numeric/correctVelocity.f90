module corr_velo_m
    use ctrl_dict, only: dim
    use parse_toml_m, only: nick
    implicit none

contains
    subroutine aver_velo(ntotal, v, mass, rho, neighborNum, neighborList, w, aver_v)
        integer, intent(in) :: ntotal
        real(8), intent(in) :: v(:, :)
        real(8), intent(in) :: mass(:)
        real(8), intent(in) :: rho(:)
        integer, intent(in) :: neighborNum(:)
        integer, intent(in) :: neighborList(:, :)
        real(8), intent(in) :: w(:, :)
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
            do k = 1, neighborNum(i)
                j = neighborList(i, k)
                
                dv = v(:, i) - v(:, j)
                aver_v(:, i) = aver_v(:, i) &
                    - epsilon               &
                    * 2 * mass(j) * dv      &
                    / (rho(i) + rho(j))     &
                    * w(i, k)
            end do
        end do

    end subroutine aver_velo

end module corr_velo_m