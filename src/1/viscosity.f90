module visc_m
    implicit none
    
contains
    subroutine viscosity(itype, eta)
        implicit none
        integer, intent(in) :: itype(:)
        real(8), intent(inout) :: eta(:)
        integer :: ntotal
        integer i

        ntotal = size(itype)

        !$OMP PARALLEL DO PRIVATE(i)
        do i = 1, ntotal
            select case(abs(itype(i)))
            case (1)
                eta(i) = 0
            case (2)
                eta(i) = 1e-3
            case (3)
                eta(i) = 0
            case (8)
                eta(i) = 8e10
            end select
        end do
        !$OMP END PARALLEL DO

    end subroutine viscosity
    
end module visc_m