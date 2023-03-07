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

        do i = 1, ntotal
            if ( abs(itype(i)) == 1 ) then
                eta(i) = 0
            else if ( abs(itype(i)) == 2 ) then
                eta(i) = 1e-3
            else if ( abs(itype(i)) == 3 ) then
                eta(i) = 0
            end if
        end do
        
    end subroutine viscosity
    
end module visc_m