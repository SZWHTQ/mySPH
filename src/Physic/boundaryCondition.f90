!! Temporarily discarded
module boundary_condition_m
    implicit none

    private

    public :: free_surface
contains
    subroutine free_surface(itype, x, rho, p)
        use parse_toml_m, only: nick
        integer, intent(in) :: itype(:)
        real(8), intent(in) :: x(:, :)
        real(8), intent(in) :: rho(:)
        real(8), intent(inout) :: p(:)
        real(8), parameter :: beta = 0.97
        real(8), parameter :: rho0 = 1000, pa = 101.325e3
        integer :: ntotal
        integer i

        ntotal = size(x, 2)

        do i = 1, ntotal
            if ( (itype(i) == 6 .or. itype(i) == 7) .and. rho(i) < 0.97 * rho0 ) then
                p(i) = pa
            end if
        end do

    end subroutine free_surface

end module boundary_condition_m