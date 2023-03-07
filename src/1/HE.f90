module he_m
    use parse_toml_m, only: nick
    implicit none
    private
    real(8), private :: D = 6930

    public detonation_wave
contains
    subroutine detonation_wave(ntotal, i_time_step, delta_t, x, itype)
        integer, intent(in)    :: ntotal
        integer, intent(in)    :: i_time_step
        real(8), intent(in)    :: delta_t
        real(8), intent(in)    :: x(:,:)
        integer, intent(inout) :: itype(:)

        select case(nick)
        case("tnt_bar")
            call tnt_bar_detonation(ntotal, i_time_step, delta_t, x, itype)
        case("undex_cylinder")
            call undex_cylinder_detonation(ntotal, i_time_step, delta_t, x, itype)
        end select

    end subroutine detonation_wave

    subroutine tnt_bar_detonation(ntotal, i_time_step, delta_t, x, itype)
        integer, intent(in)    :: ntotal
        integer, intent(in)    :: i_time_step
        real(8), intent(in)    :: delta_t
        real(8), intent(in)    :: x(:,:)
        integer, intent(inout) :: itype(:)
        integer i

        do i = 1, ntotal
            if ( x(1, i) < D*i_time_step*delta_t) then
                itype(i) = 5
            end if
        end do

    end subroutine tnt_bar_detonation

    subroutine undex_cylinder_detonation(ntotal, i_time_step, delta_t, x, itype)
        use geometry_m, only: point_t, circle_t
        integer, intent(in) :: ntotal
        integer, intent(in)    :: i_time_step
        real(8), intent(in)    :: delta_t
        real(8), intent(in)    :: x(:,:)
        integer, intent(inout) :: itype(:)
        type(circle_t) :: boundary
        integer i

        boundary = circle_t([0,0], D*i_time_step*delta_t, 0)

        do i = 1, ntotal
            if ( itype(i) == 0 .and. boundary%contain(point_t(x(:,i), i)) ) then
                itype(i) = 5
            end if
        end do
        
    end subroutine undex_cylinder_detonation

end module he_m