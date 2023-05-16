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
            case (2)
                eta(i) = 1e-3
            case (8)
                eta(i) = 0.04585
            case (101)
                eta(i) = 8e10
            case (102)
                eta(i) = 4.27e6
            case (103)
                eta(i) = 0.5e6
            case default
                eta(i) = 0
            end select
        end do
        !$OMP END PARALLEL DO

    end subroutine viscosity
    
end module visc_m