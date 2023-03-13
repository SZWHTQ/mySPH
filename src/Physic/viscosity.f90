module visc_m
    implicit none
    
contains
    subroutine viscosity(Type, eta)
        implicit none
        integer, intent(in) :: Type(:)
        real(8), intent(inout) :: eta(:)
        integer :: ntotal
        integer i

        ntotal = size(Type)

        !$OMP PARALLEL DO PRIVATE(i)
        do i = 1, ntotal
            select case(abs(Type(i)))
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