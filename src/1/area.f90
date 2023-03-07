module area_m
    implicit none
    
contains
    subroutine calculate_area(ntotal, mass, rho, area)
        integer, intent(in)  :: ntotal
        real(8), intent(in)  :: mass(:)
        real(8), intent(in)  :: rho(:)
        real(8), intent(out) :: area

        integer i

        area = 0
        do i = 1, ntotal
            area = area + mass(i) / rho(i)
        end do
        
    end subroutine calculate_area
    
end module area_m