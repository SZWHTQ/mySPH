module divergence_m
    implicit none

contains
    subroutine divergence(ntotal, F, mass, rho, neighborNum, neighborList, dwdx, div_F)
        integer, intent(in)  :: ntotal
        real(8), intent(in)  :: F(:, :)
        real(8), intent(in)  :: mass(:)
        real(8), intent(in)  :: rho(:)
        integer, intent(in)  :: neighborNum(:)
        integer, intent(in)  :: neighborList(:, :)
        real(8), intent(in)  :: dwdx(:, :, :)
        real(8), intent(inout) :: div_F(:)

        integer :: i, j, k

        forall (i=1:ntotal) div_F(i) = 0

        !$OMP PARALLEL DO PRIVATE(i, j, k)
        do i = 1, ntotal
            do k = 1, neighborNum(i)
                j = neighborList(i, k)
                div_F(i) = div_F(i) &
                    + mass(j)/rho(j) * sum((F(:, j)-F(:, i)) * dwdx(:, i, k))
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine divergence

end module divergence_m