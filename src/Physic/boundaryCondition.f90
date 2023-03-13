module boundary_condition_m
    implicit none

    private

contains
    !!! Open Boundary Non-Reflecting Condition
    subroutine NRBC(fluidDomain, boundaryDomain, &
                    ntotal, nbuffer, Type, status, x, v, mass, rho, p, e, c, hsml)
        use geometry_m
        class(geometry_t), intent(in) :: fluidDomain(:), boundaryDomain(:)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: nbuffer
        integer, intent(inout) :: Type(:)
        integer, intent(inout) :: status(:)
        real(8), intent(inout) :: x(:, :)
        real(8), intent(inout) :: v(:, :)
        real(8), intent(inout) :: mass(:)
        real(8), intent(inout) :: rho(:)
        real(8), intent(inout) :: p(:)
        real(8), intent(inout) :: e(:)
        real(8), intent(inout) :: c(:)
        real(8), intent(inout) :: hsml(:)

        integer i, j, k, d

        do i = 1, ntotal
            do j = 1, size(boundaryDomain)
                if (boundaryDomain(j)%contain(point_t(x(:, i), i))) then
                    nbuffer = nbuffer + 1
                    Type(nbuffer) = 1
                    do d = 1, 3
                        x(d, nbuffer) = x(d, i)
                        v(d, nbuffer) = v(d, i)
                    end do
                    mass(nbuffer) = mass(i)
                    rho(nbuffer) = rho(i)
                    p(nbuffer) = p(i)
                    e(nbuffer) = e(i)
                    c(nbuffer) = c(i)
                    hsml(nbuffer) = hsml(i)
                end if
            end do
        end do

    end subroutine NRBC

end module boundary_condition_m