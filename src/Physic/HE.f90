module he_m
    use ctrl_dict, only: Project
    use sph,       only: Particle
    implicit none
    private
    real(8), private :: D = 6930

    public detonation_wave
contains
    subroutine detonation_wave(i_time_step, delta_t, Particles)
        integer, intent(in) :: i_time_step
        real(8), intent(in) :: delta_t
        type(Particle), intent(inout) :: Particles(:)
        integer :: ntotal

        ntotal = size(Particles)
        select case(Project%nick)
        case("tnt_bar")
            call tnt_bar_detonation(i_time_step, delta_t, Particles)
        case("undex_cylinder")
            call undex_cylinder_detonation(i_time_step, delta_t, Particles)
        end select

    end subroutine detonation_wave

    subroutine tnt_bar_detonation(i_time_step, delta_t, Particles)
        integer, intent(in) :: i_time_step
        real(8), intent(in) :: delta_t
        type(Particle), intent(inout) :: Particles(:)
        integer :: ntotal
        integer i

        ntotal = size(Particles)
        do i = 1, ntotal
            if ( Particles(i)%x(1) < D*i_time_step*delta_t) then
                Particles(i)%Type = 5
            end if
        end do

    end subroutine tnt_bar_detonation

    subroutine undex_cylinder_detonation(i_time_step, delta_t, Particles)
        use geometry_m, only: point_t, circle_t
        integer, intent(in) :: i_time_step
        real(8), intent(in) :: delta_t
        type(Particle), intent(inout) :: Particles(:)
        integer :: ntotal
        type(circle_t) :: boundary
        integer i

        ntotal = size(Particles)
        boundary = circle_t([0,0], D*i_time_step*delta_t, 0)

        do i = 1, ntotal
            if ( Particles(i)%Type == 0 .and. boundary%contain(point_t(Particles(i)%x(:), i)) ) then
                Particles(i)%Type = 5
            end if
        end do
        
    end subroutine undex_cylinder_detonation

end module he_m