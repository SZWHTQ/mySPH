module dummy_part_m
    use ctrl_dict, only: dim, skf
    use parse_toml_m, only: nick
    implicit none
    private

    public gen_dummy_particle
contains
    subroutine gen_dummy_particle(ntotal, ndummy, itype, x, v, mass, rho, p, e, c, hsml)
        ! use ctrl_dict, only: i_time_step
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy, itype(:)
        real(8), intent(inout) :: x(:,:), v(:,:), mass(:), rho(:), p(:), e(:), c(:), hsml(:)

        ndummy = 0

        select case(nick)
        case("shear_cavity")
            call shear_cavity_dp(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml)
        case("shock_tube")
            call shock_tube_dp_1(ntotal, ndummy, itype, x, v, mass, rho, p, e, c, hsml)
            call shock_tube_dp_2(ntotal, ndummy, itype, x, v, mass, rho, p, e, c, hsml)
        case("tnt_bar")
            call tnt_bar_dp_1(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml)
            call tnt_bar_dp_2(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml)
        case("undex_chamber")
            call undex_chamber_dp_1(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml)
            call undex_chamber_dp_2(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml)
        case("UNDEX")
            call undex_dp_1(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml)
            call undex_dp_2(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml)
        case("dam_break")
            call dam_break_dp(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml)
        case("taylor_rod")
            ! call taylor_rod_dp_1(ntotal, ndummy, itype, x, v, mass, rho, p, e, c, hsml)
            call taylor_rod_dp_2(ntotal, ndummy, itype, x, v, mass, rho, p, e, c, hsml)
        end select

    end subroutine gen_dummy_particle

    subroutine shear_cavity_dp(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml)
        ! use ctrl_dict, only: i_time_step
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy, itype(:)
        real(8), intent(inout) :: x(:,:), v(:,:), mass(:), rho(:), p(:), e(:), hsml(:)
        real(8) :: xl, dx, drive
        integer :: i, mp, l, layer

        ndummy = 0
        mp     = 40
        xl     = 1.0e-3
        dx     = xl / mp
        drive  = 1.0e-3
        layer = 1

        do l = 1, layer
            !!! Monaghan type dummy particle on the Upper side
            do i = 1, 2*(mp+l) - 1
                ndummy = ndummy + 1
                x(1, ntotal+ndummy) = dx * (i-l) / 2
                x(2, ntotal+ndummy) = xl + dx / 2 * (l-1)
                v(1, ntotal+ndummy) = drive
                v(2, ntotal+ndummy) = 0
            end do

            !!! Monaghan type dummy particle on the Lower side
            do i = 1, 2*(mp+l) - 1
                ndummy = ndummy + 1
                x(1, ntotal+ndummy) = dx * (i-l) / 2
                x(2, ntotal+ndummy) = 0 - dx / 2 * (l-1)
                v(1, ntotal+ndummy) = 0
                v(2, ntotal+ndummy) = 0
            end do

            !!! Monaghan type dummy particle on the Left side
            do i = 1, 2*(mp+l) - 3
                ndummy = ndummy + 1
                x(1, ntotal+ndummy) = 0 - dx / 2 * (l-1)
                x(2, ntotal+ndummy) = dx * (i-l+1) / 2
                v(1, ntotal+ndummy) = 0
                v(2, ntotal+ndummy) = 0
            end do

            do i = 1, 2*(mp+l) - 3
                ndummy = ndummy + 1
                x(1, ntotal+ndummy) = xl + dx / 2*(l-1)
                x(2, ntotal+ndummy) = dx * (i-l+1) / 2
                v(1, ntotal+ndummy) = 0
                v(2, ntotal+ndummy) = 0
            end do
        end do

        do i = 1, ndummy
            rho(ntotal+i)   = 1000
            mass(ntotal+i)  = rho(ntotal+i) * dx**2
            p(ntotal+i)     = 0
            e(ntotal+i)     = 357.1
            itype(ntotal+i) = -itype(1)
            hsml(ntotal+i)  = dx
        end do


    end subroutine shear_cavity_dp

    subroutine shock_tube_dp_1(ntotal, ndummy, itype, x, v, mass, rho, p, e, c, hsml)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy, itype(:)
        real(8), intent(inout) :: x(:,:), v(:,:), mass(:), rho(:), p(:), e(:), c(:), hsml(:)
        real(8) :: dx, gamma
        integer :: layer = 1
        integer i

        gamma = 1.4
        dx = 0.6 / ntotal * 5

        do i = 1, layer
            ndummy = ndummy + 1
            x(1, ntotal + ndummy) = -0.6 - dx / 4 * (i+1)
            v(1, ntotal + ndummy) =  0
            hsml(ntotal + ndummy) =  dx * 2
            itype(ntotal+ ndummy) = -1
            rho(ntotal  + ndummy) =  1
            mass(ntotal + ndummy) =  rho(ntotal +  ndummy) * dx / 4
            p(ntotal    + ndummy) =  1
            e(ntotal    + ndummy) =  2.5
            c(ntotal    + ndummy) =  sqrt(gamma*p(ntotal+ndummy)/rho(ntotal+ndummy))

            ndummy = ndummy + 1
            x(1, ntotal + ndummy) =  0.6 + dx * (i+1)
            v(1, ntotal + ndummy) =  0
            hsml(ntotal + ndummy) =  dx * 2
            itype(ntotal+ ndummy) = -1
            rho(ntotal  + ndummy) =  0.25
            mass(ntotal + ndummy) =  rho(ntotal +  ndummy) * dx
            p(ntotal    + ndummy) =  0.1795
            e(ntotal    + ndummy) =  1.795
            c(ntotal    + ndummy) =  sqrt(gamma*p(ntotal+ndummy)/rho(ntotal+ndummy))
        end do

    end subroutine shock_tube_dp_1

    subroutine shock_tube_dp_2(ntotal, ndummy, itype, x, v, mass, rho, p, e, c, hsml)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy, itype(:)
        real(8), intent(inout) :: x(:,:), v(:,:), mass(:), rho(:), p(:), e(:), c(:), hsml(:)
        real(8) :: gamma
        integer :: scale_k, n_dp_1
        integer i

        n_dp_1 = ndummy
        gamma = 1.4

        select case (skf)
        case (1)
            scale_k = 2
        case (2, 3)
            scale_k = 3
        end select

        do i = 1, ntotal
            if ( x(1, i) - x(1, ntotal+n_dp_1 - 1) < hsml(i) * scale_k ) then
                ndummy = ndummy + 1

                x(1, ntotal+ndummy) =  2*x(1, ntotal+n_dp_1 - 1) - x(1, i)
                v(1, ntotal+ndummy) = -v(1, i)
                mass(ntotal+ndummy) =  mass(i)
                hsml(ntotal+ndummy) =  hsml(i)
                itype(ntotal+ndummy)= -itype(i)
                e(ntotal+ndummy)    =  e(i)
                rho(ntotal+ndummy)  =  rho(i)
                p(ntotal+ndummy)    =  p(i)
                c(ntotal+ndummy)    =  c(i)
            end if
            if ( x(1, ntotal+n_dp_1) - x(1, i) < hsml(i) * scale_k ) then
                ndummy = ndummy + 1

                x(1, ntotal+ndummy) =  2*x(1, ntotal+n_dp_1) - x(1, i)
                v(1, ntotal+ndummy) = -v(1, i)
                mass(ntotal+ndummy) =  mass(i)
                hsml(ntotal+ndummy) =  hsml(i)
                itype(ntotal+ndummy)= -itype(i)
                e(ntotal+ndummy)    =  e(i)
                rho(ntotal+ndummy)  =  rho(i)
                p(ntotal+ndummy)    =  p(i)
                c(ntotal+ndummy)    =  c(i)
            end if
        end do

    end subroutine shock_tube_dp_2

    subroutine tnt_bar_dp_1(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy, itype(:)
        real(8), intent(inout) :: x(:,:), v(:,:), mass(:), rho(:), p(:), e(:), hsml(:)
        real(8) :: space_x
        real :: rho0 = 1630, E0 = 4.29e6
        integer :: layer = 2
        integer i

        space_x = 0.1 / ntotal

        do i = 1, layer
            ndummy = ndummy + 1
            x(1, ntotal + ndummy) =  0 -  space_x * i
            v(1, ntotal + ndummy) =  0
            mass(ntotal + ndummy) =  rho0 * space_x
            hsml(ntotal + ndummy) =  space_x * 1.5
            itype(ntotal+ ndummy) =  -itype(1)
            e(ntotal   +  ndummy) =  E0
            rho(ntotal +  ndummy) =  rho0
            p(ntotal   +  ndummy) =  0
        end do

    end subroutine tnt_bar_dp_1

    subroutine tnt_bar_dp_2(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy, itype(:)
        real(8), intent(inout) :: x(:,:), v(:,:), mass(:), rho(:), p(:), e(:), hsml(:)
        integer :: scale_k, n_dp_1
        integer i

        n_dp_1 = ndummy

        select case (skf)
        case (1)
            scale_k = 2
        case (2, 3)
            scale_k = 3
        end select

        do i = 1, ntotal
            if (  x(1, i) - x(1, ntotal+n_dp_1) < hsml(i) * scale_k ) then
                ndummy = ndummy + 1

                x(1, ntotal+ndummy) = 2*x(1, ntotal+n_dp_1) - x(1, i)
                v(1, ntotal+ndummy) = -v(1, i)
                mass(ntotal+ndummy) =  mass(i)
                hsml(ntotal+ndummy) =  hsml(i)
                itype(ntotal+ndummy)= -itype(i)
                e(ntotal+ndummy)    =  e(i)
                rho(ntotal+ndummy)  =  rho(i)
                p(ntotal+ndummy)    =  p(i)

            end if
        end do

    end subroutine tnt_bar_dp_2

    subroutine undex_chamber_dp_1(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy, itype(:)
        real(8), intent(inout) :: x(:,:), v(:,:), mass(:), rho(:), p(:), e(:), hsml(:)
        integer :: n(2) = [201, 201]
        real(8) :: length(2) = [1, 1], origin(2) = [-0.5, -0.5]
        real(8) :: delta(2)
        logical :: first_entry = .true.
        integer :: i, k
        integer :: l, layer

        ndummy = 0

        delta = length / (n - 1)

        layer = 6
        do l = 1, layer
            !!! Monaghan type dummy particle on the Upper side
            do i = 1, n(1) + 2*(l+1)
                ndummy = ndummy + 1
                x(1, ntotal+ndummy) = origin(1) + delta(1) * (i-l-2)
                x(2, ntotal+ndummy) = length(2)+origin(2) + delta(2) * (l+1)
            end do

            !!! Monaghan type dummy particle on the Lower side
            do i = 1, n(1) + 2*(l+1)
                ndummy = ndummy + 1
                x(1, ntotal+ndummy) = origin(1) + delta(1) * (i-l-2)
                x(2, ntotal+ndummy) = origin(2) - delta(2) * (l+1)
            end do

            !!! Monaghan type dummy particle on the Left side
            do i = 1, n(2) + 2*l
                ndummy = ndummy + 1
                x(1, ntotal+ndummy) = origin(1) - delta(1) * (l+1)
                x(2, ntotal+ndummy) = origin(2) + delta(2) * (i-l-1)
            end do

            !!! Monaghan type dummy particle on the Right side
            do i = 1, n(2) + 2*l
                ndummy = ndummy + 1
                x(1, ntotal+ndummy) = length(1) + origin(1) + delta(1) * (l+1)
                x(2, ntotal+ndummy) = origin(2) + delta(2) * (i-l-1)
            end do
        end do

        if (first_entry) then
            do i = 1, ndummy
                k = ntotal + i
                v(:, k) = 0
                rho(k)   = 1000
                mass(k)  = rho(k) * delta(1)*delta(2)
                p(k)     = 0
                e(k)     = 0
                itype(k) = -itype(1)
                hsml(k)  = 1.5 * sum(delta)/2
            end do
            first_entry = .false.
        end if

    end subroutine undex_chamber_dp_1

    subroutine undex_chamber_dp_2(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy, itype(:)
        real(8), intent(inout) :: x(:,:), v(:,:), mass(:), rho(:), p(:), e(:), hsml(:)
        real(8) :: boundary(2, 2), dx
        logical :: first_entry = .true.
        save boundary, dx, first_entry
        integer :: scale_k, n_dp_1
        integer i

        if ( first_entry ) then
            dx = x(2, 2) - x(2, 1)
            boundary = reshape([maxval(x(:,1:ntotal+ndummy), 2)+dx/2, &
                                minval(x(:,1:ntotal+ndummy), 2)-dx/2], &
                                shape(boundary))
            first_entry = .false.
        end if

        n_dp_1 = ndummy

        select case (skf)
        case (1)
            scale_k = 2
        case (2, 3)
            scale_k = 3
        end select

        do i = 1, ntotal
            !!! Dummy particle II on the Upper side
            if ( boundary(2, 1) - x(2, i) < hsml(i) * scale_k ) then
                ndummy = ndummy + 1
                x(1, ntotal+ndummy) =  x(1, i)
                x(2, ntotal+ndummy) =  2*boundary(2, 1) - x(2, i)
                v(:, ntotal+ndummy) =  v(:, i)
                mass(ntotal+ndummy) =  mass(i)
                hsml(ntotal+ndummy) =  hsml(i)
                itype(ntotal+ndummy)=  itype(i)
                e(ntotal+ndummy)    =  e(i)
                rho(ntotal+ndummy)  =  rho(i)
                p(ntotal+ndummy)    =  p(i)
            end if

            !!! Dummy particle II on the Lower side
            if ( x(2, i) - boundary(2, 2) < hsml(i) * scale_k ) then
                ndummy = ndummy + 1
                x(1, ntotal+ndummy) =  x(1, i)
                x(2, ntotal+ndummy) =  2*boundary(2, 2) - x(2, i)
                v(:, ntotal+ndummy) =  v(:, i)
                mass(ntotal+ndummy) =  mass(i)
                hsml(ntotal+ndummy) =  hsml(i)
                itype(ntotal+ndummy)=  itype(i)
                e(ntotal+ndummy)    =  e(i)
                rho(ntotal+ndummy)  =  rho(i)
                p(ntotal+ndummy)    =  p(i)
            end if

            !!! Monaghan type dummy particle on the Right side
            if ( boundary(1, 1) - x(1, i) < hsml(i) * scale_k ) then
                ndummy = ndummy + 1
                x(1, ntotal+ndummy) =  2*boundary(1, 1) - x(1, i)
                x(2, ntotal+ndummy) =  x(2, i)
                v(:, ntotal+ndummy) =  v(:, i)
                mass(ntotal+ndummy) =  mass(i)
                hsml(ntotal+ndummy) =  hsml(i)
                itype(ntotal+ndummy)=  itype(i)
                e(ntotal+ndummy)    =  e(i)
                rho(ntotal+ndummy)  =  rho(i)
                p(ntotal+ndummy)    =  p(i)
            end if

            !!! Monaghan type dummy particle on the Left side
            if ( x(1, i) - boundary(1, 2) < hsml(i) * scale_k ) then
                ndummy = ndummy + 1
                x(1, ntotal+ndummy) =  2*boundary(1, 2) - x(1, i)
                x(2, ntotal+ndummy) =  x(2, i)
                v(:, ntotal+ndummy) =  v(:, i)
                mass(ntotal+ndummy) =  mass(i)
                hsml(ntotal+ndummy) =  hsml(i)
                itype(ntotal+ndummy)=  itype(i)
                e(ntotal+ndummy)    =  e(i)
                rho(ntotal+ndummy)  =  rho(i)
                p(ntotal+ndummy)    =  p(i)
            end if
        end do

    end subroutine undex_chamber_dp_2

    subroutine undex_dp_1(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy, itype(:)
        real(8), intent(inout) :: x(:,:), v(:,:), mass(:), rho(:), p(:), e(:), hsml(:)
        real(8) :: wall_domain(4) = [0., 0., -0.05, 0.05]
        real(8) :: dx = 5e-4
        integer :: Ny
        integer :: layer

        integer :: index
        integer i, l

        ndummy = 0
        layer  = 4

        Ny = floor((wall_domain(4) - wall_domain(3)) / dx)
        do l = 1, layer
            do i = 1, Ny
                ndummy = ndummy + 1
                index  = ntotal + ndummy
                x(:, index)  = [wall_domain(1) - (l - 1) * dx, &
                                wall_domain(3) + (i - 1) * dx]
                v(:, index)  = 0
                rho(index)   = 1000
                mass(index)  = rho(index) * dx * dx
                p(index)     = 0
                e(index)     = 2e7
                itype(index) = -6
                hsml(index)  = dx
            end do
        end do

    end subroutine undex_dp_1

    subroutine undex_dp_2(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy, itype(:)
        real(8), intent(inout) :: x(:,:), v(:,:), mass(:), rho(:), p(:), e(:), hsml(:)
        real(8), save :: boundary(2, 2)
        integer :: scale_k, n_dp_1
        logical, save :: first_entry = .true.
        integer i

        if ( first_entry ) then
            boundary = reshape([maxval(x(:,1:ntotal), 2), &
                                minval(x(:,1:ntotal), 2)], &
                                shape(boundary))
            boundary = boundary * 1.01
            first_entry = .false.
        end if

        n_dp_1 = ndummy

        select case (skf)
        case (1)
            scale_k = 2
        case (2, 3)
            scale_k = 3
        end select

        do i = 1, ntotal
            !!! Dummy particle II Upside
            if ( boundary(2, 1) - x(2, i) < hsml(i) * scale_k ) then
                ndummy = ndummy + 1
                x(1, ntotal+ndummy) =  x(1, i)
                x(2, ntotal+ndummy) =  2*boundary(2, 1) - x(2, i)
                v(:, ntotal+ndummy) =  v(:, i)
                mass(ntotal+ndummy) =  mass(i)
                hsml(ntotal+ndummy) =  hsml(i)
                itype(ntotal+ndummy)=  itype(i)
                e(ntotal+ndummy)    =  e(i)
                rho(ntotal+ndummy)  =  rho(i)
                p(ntotal+ndummy)    =  p(i)
            end if

            !!! Dummy particle II Downside
            if ( x(2, i) - boundary(2, 2) < hsml(i) * scale_k ) then
                ndummy = ndummy + 1
                x(1, ntotal+ndummy) =  x(1, i)
                x(2, ntotal+ndummy) =  2*boundary(2, 2) - x(2, i)
                v(:, ntotal+ndummy) =  v(:, i)
                mass(ntotal+ndummy) =  mass(i)
                hsml(ntotal+ndummy) =  hsml(i)
                itype(ntotal+ndummy)=  itype(i)
                e(ntotal+ndummy)    =  e(i)
                rho(ntotal+ndummy)  =  rho(i)
                p(ntotal+ndummy)    =  p(i)
            end if

            !!! Monaghan type dummy particle on the Right side
            if ( boundary(1, 1) - x(1, i) < hsml(i) * scale_k ) then
                ndummy = ndummy + 1
                x(1, ntotal+ndummy) =  2*boundary(1, 1) - x(1, i)
                x(2, ntotal+ndummy) =  x(2, i)
                v(:, ntotal+ndummy) =  v(:, i)
                mass(ntotal+ndummy) =  mass(i)
                hsml(ntotal+ndummy) =  hsml(i)
                itype(ntotal+ndummy)=  itype(i)
                e(ntotal+ndummy)    =  e(i)
                rho(ntotal+ndummy)  =  rho(i)
                p(ntotal+ndummy)    =  p(i)
            end if

            ! !!! Monaghan type dummy particle on the Left side
            ! if ( x(1, i) - boundary(1, 2) < hsml(i) * scale_k ) then
            !     ndummy = ndummy + 1
            !     x(1, ntotal+ndummy) =  2*boundary(1, 2) - x(1, i)
            !     x(2, ntotal+ndummy) =  x(2, i)
            !     v(:, ntotal+ndummy) =  v(:, i)
            !     mass(ntotal+ndummy) =  mass(i)
            !     hsml(ntotal+ndummy) =  hsml(i)
            !     itype(ntotal+ndummy)=  itype(i)
            !     e(ntotal+ndummy)    =  e(i)
            !     rho(ntotal+ndummy)  =  rho(i)
            !     p(ntotal+ndummy)    =  p(i)
            ! end if
        end do

    end subroutine undex_dp_2

    subroutine dam_break_dp(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy, itype(:)
        real(8), intent(inout) :: x(:,:), v(:,:), mass(:), rho(:), p(:), e(:), hsml(:)
        real(8) :: dx
        save dx
        real(8) :: wall_domain(4)
        integer :: Nx, Ny
        integer :: layer
        logical :: first_entry = .true.
        save first_entry

        integer :: index
        integer i, l

        ndummy = 0

        if ( first_entry ) then
            dx = abs(x(1, 2) - x(2, 2))
            first_entry = .false.
        end if
        ! wall_domain = [0, 0, 0, 3]
        wall_domain = [0, 100, 0, 40]
        ! wall_domain(2) = floor(5.366/dx)*dx

        layer = 4
        do l = 1, layer
            !!! Dummy particle I on the Left side
            Ny = floor((wall_domain(4) - wall_domain(3))/dx) + layer
            do i = 1, Ny
                ndummy = ndummy + 1
                index = ntotal + ndummy
                x(:, index)  = [wall_domain(1) - (l-0.5)*dx, wall_domain(3) + (i-layer-0.5)*dx]
                v(:, index)  = 0
                rho(index)   = 1000
                mass(index)  = rho(index) * dx * dx
                p(index)     = 0
                e(index)     = 0
                itype(index) = -itype(1)
                hsml(index)  = dx
            end do

            !!! Dummy particle I on the Right side
            do i = 1, Ny
                ndummy = ndummy + 1
                index = ntotal + ndummy
                x(:, index)  = [wall_domain(2) + (l-0.5)*dx, &
                                wall_domain(3) + (i-layer-0.5)*dx]
                v(:, index)  = 0
                rho(index)   = 1000
                mass(index)  = rho(index) * dx * dx
                p(index)     = 0
                e(index)     = 0
                itype(index) = -itype(1)
                hsml(index)  = dx
            end do

            !!! Dummy particle I on the Bottom
            Nx = floor((wall_domain(2) - wall_domain(1))/dx) + 1
            do i = 1, Nx
                ndummy = ndummy + 1
                index = ntotal + ndummy
                x(:, index)  = [wall_domain(1) + (i-0.5)*dx, &
                                wall_domain(3) - (l-0.5)*dx]
                v(:, index)  = 0
                rho(index)   = 1000
                mass(index)  = rho(index) * dx * dx
                p(index)     = 0
                e(index)     = 0
                itype(index) = -itype(1)
                hsml(index)  = dx
            end do
        end do

        !!! Dummy particle I for baffle
        Nx = floor(3 / (dx / 2))
        do l = 1, layer
            do i = 1, Nx - l + 1
                ndummy = ndummy + 1
                index = ntotal + ndummy
                x(:, index)  = [sum(wall_domain(1:2))/2 + (i+l-1.5)*dx/2, &
                                wall_domain(3) + (i-0.5)*dx/2]
                v(:, index)  = 0
                rho(index)   = 1000
                mass(index)  = rho(index) * dx * dx / 4
                p(index)     = 0
                e(index)     = 0
                itype(index) = -itype(1)
                hsml(index)  = dx
            end do
        end do

    end subroutine dam_break_dp

    subroutine taylor_rod_dp_1(ntotal, ndummy, itype, x, v, mass, rho, p, e, c, hsml)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy, itype(:)
        real(8), intent(inout) :: x(:,:), v(:,:), mass(:), rho(:), p(:), e(:), c(:), hsml(:)
        real(8) :: dx
        save dx
        real(8) :: wall_domain(4)
        integer :: Nx, Ny
        integer :: layer
        logical :: first_entry = .true.
        save first_entry

        integer :: index
        integer i, j

        ndummy = 0

        if ( first_entry ) then
            dx = abs(x(1, 2) - x(2, 2)) / 2
            first_entry = .false.
        end if

        layer = 4
        wall_domain = [-1.14, 1.9, 0., 0.]  * 1e-2
        wall_domain(3) = - layer * dx

        Nx = floor((wall_domain(2) - wall_domain(1))/dx)
        Ny = floor((wall_domain(4) - wall_domain(3))/dx)
        do i = 1, Nx
            do j = 1, Ny
                ndummy = ndummy + 1
                index = ntotal + ndummy
                x(:, index)  = [wall_domain(2) - (i-0.5)*dx, &
                                wall_domain(4) - (j-0.5)*dx]
                v(:, index)  = 0
                rho(index)   = 7850
                mass(index)  = rho(index) * dx * dx
                p(index)     = 0
                e(index)     = 0
                c(index)     = 5000
                itype(index) = -itype(1)
                hsml(index)  = dx * 2
            end do
        end do

    end subroutine taylor_rod_dp_1

    subroutine taylor_rod_dp_2(ntotal, ndummy, itype, x, v, mass, rho, p, e, c, hsml)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy, itype(:)
        real(8), intent(inout) :: x(:,:), v(:,:), mass(:), rho(:), p(:), e(:), c(:), hsml(:)
        integer :: scale_k

        integer i

        select case (skf)
        case (1)
            scale_k = 2
        case (2, 3)
            scale_k = 3
        end select

        do i = 1, ntotal
            if ( x(2, i) - 0 < hsml(i) * scale_k ) then
                ndummy = ndummy + 1
                x(1, ntotal+ndummy)  =  x(1, i)
                x(2, ntotal+ndummy)  = -x(2, i)
                v(1, ntotal+ndummy)  =  v(1, i)
                v(2, ntotal+ndummy)  = -v(2, i)
                mass(ntotal+ndummy)  =  mass(i)
                hsml(ntotal+ndummy)  =  hsml(i)
                itype(ntotal+ndummy) = -itype(i)
                e(ntotal+ndummy)     =  e(i)
                rho(ntotal+ndummy)   =  rho(i)
                p(ntotal+ndummy)     =  p(i)
                c(ntotal+ndummy)     =  c(i)
            end if
        end do

    end subroutine taylor_rod_dp_2

end module dummy_part_m