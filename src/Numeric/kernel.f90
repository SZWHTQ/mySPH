module kernel_m
    use ctrl_dict, only: Field, Config
    use tools_m,   only: PI

    implicit none

contains
    pure subroutine kernel(r, dx, hsml, w, dwdx)
        implicit none
        real(8), intent(in)  :: r        !! Distance between particles i and j
        real(8), intent(in)  :: dx(:)    !! x-, y- and z-distance between i and j
        real(8), intent(in)  :: hsml     !! Smoothing length
        real(8), intent(inout) :: w        !! Kernel for all interaction pairs
        real(8), intent(inout) :: dwdx(:)  !! Derivative of kernel with respect to x, y and z

        select case (Config%skf)
        case (1)
            call cubic_spline_function(r, dx, hsml, w, dwdx)
        case (2)
            call gauss_kernel_function(r, dx, hsml, w, dwdx)
        case (3)
            call quintic_kernel_function(r, dx, hsml, w, dwdx)
        ! case default
        !     call print_error(skf, "Unsupported SKF", type="value")
        !     error stop "at subroutine kernel"
        end select

    end subroutine kernel

    pure subroutine cubic_spline_function(r, dx, hsml, w, dwdx)
        implicit none
        real(8), intent(in)  :: r        !! Distance between particles i and j
        real(8), intent(in)  :: dx(:)    !! x-, y- and z-distance between i and j
        real(8), intent(in)  :: hsml     !! Smoothing length
        real(8), intent(inout) :: w        !! Kernel for all interaction pairs
        real(8), intent(inout) :: dwdx(:)  !! Derivative of kernel with respect to x, y and z

        real(8) :: q       !! R = r/h = |x-x'|/h
        real(8) :: factor  !! alpha_d

        q = r / hsml
        w = 0
        dwdx = 0

        !!! Calculate factor alpha_d
        if (Field%dim == 1) then
            factor = 1.0_8 / (hsml)
        else if (Field%dim == 2) then
            factor = 15._8 / (7._8*PI*hsml**2)
        else if (Field%dim == 3) then
            factor = 3.0_8 / (2._8*PI*hsml**3)
        ! else
        !     call print_error(Field%dim, "Unsupported Dimension", type="value")
        !     error stop "at subroutine cubic_spline_function"
        end if

        !!! Main Function
        if (q >= 0 .and. q < 1) then
            w = factor * ( 2._8/3._8 - q*q + (q**3)/2._8 )
            dwdx = factor * ( -2._8 + 3._8/2._8*q ) / hsml**2 * dx
        else if ( q >= 1 .and. q <= 2 ) then
            w = factor * ( 1._8/6._8 ) * (2-q)**3
            dwdx = factor * -(1._8/2._8) * (2-q)**2 / hsml * (dx/r)
        else
            w = 0
            dwdx = 0
        end if


    end subroutine cubic_spline_function

    pure subroutine gauss_kernel_function(r, dx, hsml, w, dwdx)
        implicit none
        real(8), intent(in)  :: r        !! Distance between particles i and j
        real(8), intent(in)  :: dx(:)    !! x-, y- and z-distance between i and j
        real(8), intent(in)  :: hsml     !! Smoothing length
        real(8), intent(inout) :: w        !! Kernel for all interaction pairs
        real(8), intent(inout) :: dwdx(:)  !! Derivative of kernel with respect to x, y and z

        real(8) :: q       !! R = r/h = |x-x'|/h
        real(8) :: factor  !! alpha_d

        q = r/hsml
        w = 0
        dwdx = 0

        !!! Calculate factor alpha_d
        if (1 <= Field%dim .and. Field%dim <= 3) then ! NOTE: Var "dim" is an integer
            factor = 1.0_8 / (sqrt(PI)*hsml)**Field%dim
        ! else
        !     call print_error(Field%dim, "Unsupported Dimension", type="value")
        !     error stop "at subroutine gauss_kernel_function"
        end if

        !!! Main Function
        if (q >= 0 .and. q < 3) then
            w = factor * exp(-q*q)
            dwdx = w * ( -2._8 * dx/(hsml**2))
        else
            w = 0._8
            dwdx = 0._8
        end if


    end subroutine gauss_kernel_function

    pure subroutine quintic_kernel_function(r, dx, hsml, w, dwdx)
        implicit none
        real(8), intent(in)  :: r        !! Distance between particles i and j
        real(8), intent(in)  :: dx(:)    !! x-, y- and z-distance between i and j
        real(8), intent(in)  :: hsml     !! Smoothing length
        real(8), intent(inout) :: w        !! Kernel for all interaction pairs
        real(8), intent(inout) :: dwdx(:)  !! Derivative of kernel with respect to x, y and z

        real(8) :: q       !! R = r/h = |x-x'|/h
        real(8) :: factor  !! alpha_d

        q = r/hsml
        w = 0
        dwdx = 0

        !!! Calculate factor alpha_d
        if (Field%dim == 1) then
            factor = 1.0_8 / (120._8*hsml)
        else if (Field%dim == 2) then
            factor = 7.0_8 / (478._8*PI*hsml**2)
        else if (Field%dim == 3) then
            factor = 1.0_8 / (120._8*PI*hsml**3)
            ! factor = 3.0_8 / (359._8*PI*hsml**3) !! Chap.3 P90
        ! else
        !     call print_error(Field%dim, "Unsupported Dimension", type="value")
        !     error stop "at subroutine quintic_kernel_function"
        end if

        !!! Main Function
        if (q >= 0 .and. q < 1) then
            w = factor * ( (3-q)**5 - 6*(2-q)**5 + 15*(1-q)**5 )
            dwdx = factor * (-120 + 120*q - 50*q**2) / hsml**2 * dx
        else if ( q >= 1 .and. q < 2 ) then
            w = factor * ( (3-q)**5 - 6*(2-q)**5 )
            dwdx = factor * (-5*(3-q)**4 + 30*(2-q)**4) / hsml * (dx/r)
        else if ( q >= 2 .and. q <= 3) then
            w = factor * (3-q)**5
            dwdx = factor * (-5*(3-q)**4) / hsml * (dx/r)
        else
            w = 0
            dwdx = 0
        end if


    end subroutine quintic_kernel_function

end module kernel_m