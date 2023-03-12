module boundary_condition_m
    implicit none

    private

contains
    subroutine NRBC(area, ntotal, nbuffer, itype, x, v, mass, rho, p, e, c, hsml)
        use geometry_m
        class(geometry_t), intent(in) :: area(:)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: nbuffer
        integer, intent(inout) :: itype(:)
        real(8), intent(inout) :: x(:, :)
        real(8), intent(inout) :: v(:, :)
        real(8), intent(inout) :: mass(:)
        real(8), intent(inout) :: rho(:)
        real(8), intent(inout) :: p(:)
        real(8), intent(inout) :: e(:)
        real(8), intent(inout) :: c(:)
        real(8), intent(inout) :: hsml(:)
        
    
        
    end subroutine NRBC

end module boundary_condition_m