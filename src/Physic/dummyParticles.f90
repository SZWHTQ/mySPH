module dummy_part_m
    use ctrl_dict, only: Field, Config, Project
    use sph
    implicit none
    private

    public gen_dummy_particle
contains
    subroutine gen_dummy_particle(ntotal, ndummy, Particles)
        ! use ctrl_dict, only: i_time_step
        integer, intent(in)    :: ntotal
        integer, intent(inout) :: ndummy
        type(Particle), intent(inout) :: Particles(:)

        ndummy = 0

        select case(Project%nick)
        case("shear_cavity")
            call shear_cavity_dp(ntotal, ndummy, Particles)
        case("shock_tube")
            call shock_tube_dp_1(ntotal, ndummy, Particles)
            call shock_tube_dp_2(ntotal, ndummy, Particles)
        case("tnt_bar")
            call tnt_bar_dp_1(ntotal, ndummy, Particles)
            call tnt_bar_dp_2(ntotal, ndummy, Particles)
        case("undex_chamber")
            call undex_chamber_dp_1(ntotal, ndummy, Particles)
            ! call undex_chamber_dp_2(ntotal, ndummy, Particles)
        case("UNDEX")
            call undex_dp_1(ntotal, ndummy, Particles)
            call undex_dp_2(ntotal, ndummy, Particles)
        case("dam_break")
            call dam_break_dp(ntotal, ndummy, Particles)
        case("taylor_rod")
            ! call taylor_rod_dp_1(ntotal, ndummy, Particles)
            call taylor_rod_dp_2(ntotal, ndummy, Particles)
        end select

    end subroutine gen_dummy_particle

    subroutine shear_cavity_dp(ntotal, ndummy, P)
        ! use ctrl_dict, only: i_time_step
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy
        type(Particle), intent(inout) :: P(:)
        real(8) :: xl, dx, drive
        integer :: i, mp, l, layer, index

        ndummy = 0
        mp     = 40
        xl     = 1.0e-3
        dx     = xl / mp
        if ( Config%i_time_step <= 10000 ) then
            drive = 0
        else
            drive = 1.0e-3
        end if
        layer = 1

        do l = 1, layer
            !!! Monaghan type dummy particle on the Upper side
            do i = 1, 2*(mp+l) - 1
                ndummy = ndummy + 1
                index = ntotal+ndummy
                P(index)%x(1) = dx * (i-l) / 2
                P(index)%x(2) = xl + dx / 2 * (l-1)
                P(index)%v(1) = drive
                P(index)%v(2) = 0
            end do

            !!! Monaghan type dummy particle on the Lower side
            do i = 1, 2*(mp+l) - 1
                ndummy = ndummy + 1
                index = ntotal+ndummy
                P(index)%x(1) = dx * (i-l) / 2
                P(index)%x(2) = 0 - dx / 2 * (l-1)
                P(index)%v(1) = 0
                P(index)%v(2) = 0
            end do

            !!! Monaghan type dummy particle on the Left side
            do i = 1, 2*(mp+l) - 3
                ndummy = ndummy + 1
                index = ntotal+ndummy
                P(index)%x(1) = 0 - dx / 2 * (l-1)
                P(index)%x(2) = dx * (i-l+1) / 2
                P(index)%v(1) = 0
                P(index)%v(2) = 0
            end do

            do i = 1, 2*(mp+l) - 3
                ndummy = ndummy + 1
                index = ntotal+ndummy
                P(index)%x(1) = xl + dx / 2*(l-1)
                P(index)%x(2) = dx * (i-l+1) / 2
                P(index)%v(1) = 0
                P(index)%v(2) = 0
            end do
        end do

        do i = 1, ndummy
            P(ntotal+i)%Density         = 1000
            P(ntotal+i)%Mass            = P(ntotal+i)%Density * dx**2
            P(ntotal+i)%Pressure        = 0
            P(ntotal+i)%InternalEnergy  = 357.1
            P(ntotal+i)%Type            = -P(1)%Type
            P(ntotal+i)%SmoothingLength = dx
        end do


    end subroutine shear_cavity_dp

    subroutine shock_tube_dp_1(ntotal, ndummy, P)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy
        type(Particle), intent(inout) :: P(:)
        real(8) :: dx, gamma
        integer :: layer = 1
        integer i

        gamma = 1.4
        dx = 0.6 / ntotal * 5

        do i = 1, layer
            ndummy = ndummy + 1
            P(ntotal+ndummy)%x(1)            = -0.6 - dx / 4 * (i+1)
            P(ntotal+ndummy)%v(1)            =  0
            P(ntotal+ndummy)%SmoothingLength =  dx * 1.5
            P(ntotal+ndummy)%Type            = -1
            P(ntotal+ndummy)%Density         =  1
            P(ntotal+ndummy)%Mass            =  P(ntotal+ndummy)%Density * dx / 4
            P(ntotal+ndummy)%Pressure        =  1
            P(ntotal+ndummy)%InternalEnergy  =  2.5
            P(ntotal+ndummy)%SoundSpeed      = sqrt(gamma                       &
                                                    * P(ntotal+ndummy)%Pressure &
                                                    / P(ntotal+ndummy)%Density)

            ndummy = ndummy + 1
            P(ntotal+ndummy)%x(1)            =  0.6 + dx * (i+1)
            P(ntotal+ndummy)%v(1)            =  0
            P(ntotal+ndummy)%SmoothingLength =  dx * 2
            P(ntotal+ndummy)%Type            = -1
            P(ntotal+ndummy)%Density         =  0.25
            P(ntotal+ndummy)%Mass            =  P(ntotal+ndummy)%Density * dx
            P(ntotal+ndummy)%Pressure        =  0.1795
            P(ntotal+ndummy)%InternalEnergy  =  1.795
            P(ntotal+ndummy)%SoundSpeed      =  sqrt(gamma                      &
                                                    * P(ntotal+ndummy)%Pressure &
                                                    / P(ntotal+ndummy)%Density)
        end do

    end subroutine shock_tube_dp_1

    subroutine shock_tube_dp_2(ntotal, ndummy, P)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy
        type(Particle), intent(inout) :: P(:)
        real(8) :: gamma
        integer :: scale_k, n_dp_1
        integer i, index

        n_dp_1 = ndummy
        gamma = 1.4
        scale_k = 0

        select case (Config%skf)
        case (1)
            scale_k = 2
        case (2, 3)
            scale_k = 3
        end select

        do i = 1, ntotal
            if ( P(i)%x(1) - P(ntotal+n_dp_1 - 1)%x(1) &
               < P(i)%SmoothingLength * scale_k ) then
                ndummy = ndummy + 1
                index = ntotal + ndummy
                P(index)%x(1) =  2*P(ntotal+n_dp_1 - 1)%x(1) - P(i)%x(1)
                P(index)%v(1) = -P(i)%v(1)
                P(index)%Mass = P(i)%Mass
                P(index)%SmoothingLength = P(i)%SmoothingLength
                P(index)%Type = -P(i)%Type
                P(index)%InternalEnergy = P(i)%InternalEnergy
                P(index)%Pressure = P(i)%Pressure
                P(index)%Density = P(i)%Density
                P(index)%SoundSpeed = P(i)%SoundSpeed
            end if
            if ( P(ntotal+n_dp_1)%x(1) - P(i)%x(1) &
               < P(i)%SmoothingLength * scale_k ) then
                ndummy = ndummy + 1
                index = ntotal + ndummy
                P(index)%x(1) =  2*P(ntotal+n_dp_1)%x(1) - P(i)%x(1)
                P(index)%v(1) = -P(i)%v(1)
                P(index)%Mass = P(i)%Mass
                P(index)%SmoothingLength = P(i)%SmoothingLength
                P(index)%Type = -P(i)%Type
                P(index)%InternalEnergy = P(i)%InternalEnergy
                P(index)%Pressure = P(i)%Pressure
                P(index)%Density = P(i)%Density
                P(index)%SoundSpeed = P(i)%SoundSpeed
            end if
        end do

    end subroutine shock_tube_dp_2

    subroutine tnt_bar_dp_1(ntotal, ndummy, P)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy
        type(Particle), intent(inout) :: P(:)
        real(8) :: space_x
        real :: rho0 = 1630, E0 = 4.29e6
        integer :: layer = 2
        integer i

        space_x = 0.1 / ntotal

        do i = 1, layer
            ndummy = ndummy + 1
            P(ntotal+ndummy)%x(1)            =  0 -  space_x * i
            P(ntotal+ndummy)%v(1)            =  0
            P(ntotal+ndummy)%Mass            =  rho0 * space_x
            P(ntotal+ndummy)%SmoothingLength =  space_x * 1.5
            P(ntotal+ndummy)%Type           =  -P(1)%Type
            P(ntotal+ndummy)%InternalEnergy  =  E0
            P(ntotal+ndummy)%Density         =  rho0
            P(ntotal+ndummy)%Pressure        =  0
        end do

    end subroutine tnt_bar_dp_1

    subroutine tnt_bar_dp_2(ntotal, ndummy, P)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy
        type(Particle), intent(inout) :: P(:)
        integer :: scale_k, n_dp_1
        integer i

        n_dp_1 = ndummy
        scale_k = 0

        select case (Config%skf)
        case (1)
            scale_k = 2
        case (2, 3)
            scale_k = 3
        end select

        do i = 1, ntotal
            if (  P(i)%x(1) - P(ntotal+n_dp_1)%x(1) &
                < P(i)%SmoothingLength * scale_k ) then
                ndummy = ndummy + 1
                P(ntotal+ndummy) = P(i)
                P(ntotal+ndummy)%x(1) = 2*P(ntotal+n_dp_1)%x(1) - P(i)%x(1)
                P(ntotal+ndummy)%v(1) = -P(i)%v(1)
                P(ntotal+ndummy)%Type= -P(i)%Type

            end if
        end do

    end subroutine tnt_bar_dp_2

    subroutine undex_chamber_dp_1(ntotal, ndummy, P)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy
        type(Particle), intent(inout) :: P(:)
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
                P(ntotal+ndummy)%x(1) = origin(1) + delta(1) * (i-l-2)
                P(ntotal+ndummy)%x(2) = length(2)+origin(2) + delta(2) * (l+1)
            end do

            !!! Monaghan type dummy particle on the Lower side
            do i = 1, n(1) + 2*(l+1)
                ndummy = ndummy + 1
                P(ntotal+ndummy)%x(1) = origin(1) + delta(1) * (i-l-2)
                P(ntotal+ndummy)%x(2) = origin(2) - delta(2) * (l+1)
            end do

            !!! Monaghan type dummy particle on the Left side
            do i = 1, n(2) + 2*l
                ndummy = ndummy + 1
                P(ntotal+ndummy)%x(1) = origin(1) - delta(1) * (l+1)
                P(ntotal+ndummy)%x(2) = origin(2) + delta(2) * (i-l-1)
            end do

            !!! Monaghan type dummy particle on the Right side
            do i = 1, n(2) + 2*l
                ndummy = ndummy + 1
                P(ntotal+ndummy)%x(1) = length(1) + origin(1) + delta(1) * (l+1)
                P(ntotal+ndummy)%x(2) = origin(2) + delta(2) * (i-l-1)
            end do
        end do

        if (first_entry) then
            do i = 1, ndummy
                k = ntotal + i
                P(k)%v(:)            = 0
                P(k)%Density         = 1000
                P(k)%Mass            = P(k)%Density * delta(1)*delta(2)
                P(k)%Pressure        = 0
                P(k)%InternalEnergy  = 0
                P(k)%Type            = -P(1)%Type
                P(k)%SmoothingLength = 1.5 * sum(delta)/2
            end do
            first_entry = .false.
        end if

    end subroutine undex_chamber_dp_1

    ! subroutine undex_chamber_dp_2(ntotal, ndummy, P)
    !     integer, intent(in) :: ntotal
    !     integer, intent(inout) :: ndummy
    !     type(Particle), intent(inout) :: P(:)
    !     real(8) :: boundary(2, 2), dx, thickness
    !     logical :: first_entry = .true.
    !     save boundary, dx, first_entry
    !     integer :: scale_k, n_dp_1
    !     integer i

    !     scale_k = 0
    !     boundary = reshape([-huge(0._8), -huge(0._8),  &
    !                          huge(0._8),  huge(0._8)], &
    !                          shape(boundary))
    !     if ( first_entry ) then
    !         dx = P(2)%x(2) - P(1)%x(2)
    !         do i = 1, ntotal+ndummy
    !             if ( P(i)%x(1) > boundary(1,1) ) then
    !                 boundary(1,1) = P(i)%x(1)
    !             end if
    !             if ( p(i)%x(2) > boundary(2,1) ) then
    !                 boundary(2,1) = P(i)%x(2)
    !             end if
    !             if ( p(i)%x(1) < boundary(1,2) ) then
    !                 boundary(1,2) = P(i)%x(1)
    !             end if
    !             if ( p(i)%x(2) < boundary(2,2) ) then
    !                 boundary(2,2) = P(i)%x(2)
    !             end if
    !         end do
    !         first_entry = .false.
    !     end if

    !     n_dp_1 = ndummy

    !     select case (Config%skf)
    !     case (1)
    !         scale_k = 2
    !     case (2, 3)
    !         scale_k = 3
    !     end select

    !     do i = 1, ntotal
    !         !!! Dummy particle II on the Upper side
    !         if ( boundary(2, 1) - P(i)%x(2) &
    !            < P(i)%SmoothingLength * scale_k ) then
    !             ndummy = ndummy + 1
    !             thickness = P(i)%SmoothingLength * scale_k
    !             P(ntotal+ndummy) = P(i)
    !             ! P(ntotal+ndummy)%x(2) =  2*boundary(2, 1) - P(i)%x(2)
    !             P(ntotal+ndummy)%x(2) =  P(i)%x(2) + thickness
    !             ! P(ntotal+ndummy)%v(:) =  P(i)%v(:)
    !             ! P(ntotal+ndummy)%Type =  P(i)%Type
    !         end if

    !         !!! Dummy particle II on the Lower side
    !         if ( P(i)%x(2) - boundary(2, 2) < P(i)%SmoothingLength * scale_k ) then
    !             ndummy = ndummy + 1
    !             thickness = P(i)%SmoothingLength * scale_k
    !             P(ntotal+ndummy) = P(i)
    !             ! P(ntotal+ndummy)%x(2) =  2*boundary(2, 2) - P(i)%x(2)
    !             P(ntotal+ndummy)%x(2) =  P(i)%x(2) - thickness
    !             ! P(ntotal+ndummy)%v(:) =  P(i)%v(:)
    !             ! P(ntotal+ndummy)%Type =  P(i)%Type
    !         end if
    !     end do

    !     n_dp_1 = ntotal+ndummy

    !     do i = 1, n_dp_1
    !         !!! Monaghan type dummy particle on the Right side
    !         if ( boundary(1, 1) - P(i)%x(1) < P(i)%SmoothingLength * scale_k ) then
    !             ndummy = ndummy + 1
    !             thickness = P(i)%SmoothingLength * scale_k
    !             P(ntotal+ndummy) = P(i)
    !             ! P(ntotal+ndummy)%x(1) =  2*boundary(1, 1) - P(i)%x(1)
    !             P(ntotal+ndummy)%x(1) =  P(i)%x(1) + thickness
    !             ! P(ntotal+ndummy)%v(:) =  P(i)%v(:)
    !             ! P(ntotal+ndummy)%Type =  P(i)%Type
    !         end if

    !         !!! Monaghan type dummy particle on the Left side
    !         if ( P(i)%x(1) - boundary(1, 2) < P(i)%SmoothingLength * scale_k ) then
    !             thickness = P(i)%SmoothingLength * scale_k
    !             ndummy = ndummy + 1
    !             P(ntotal+ndummy) = P(i)
    !             ! P(ntotal+ndummy)%x(1) =  2*boundary(1, 2) - P(i)%x(1)
    !             P(ntotal+ndummy)%x(1) =  P(i)%x(1) - thickness
    !             ! P(ntotal+ndummy)%v(:) =  P(i)%v(:)
    !             ! P(ntotal+ndummy)%Type =  P(i)%Type
    !         end if
    !     end do

    ! end subroutine undex_chamber_dp_2

    subroutine undex_dp_1(ntotal, ndummy, P)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy
        type(Particle), intent(inout) :: P(:)
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
                P(index)%x(:)  = [wall_domain(1) - (l - 1) * dx, &
                                wall_domain(3) + (i - 1) * dx]
                P(index)%v(:)            = 0
                P(index)%Density         = 1000
                P(index)%Mass            = P(index)%Density * dx * dx
                P(index)%Pressure        = 0
                P(index)%InternalEnergy  = 2e7
                P(index)%Type           = -6
                P(index)%SmoothingLength = dx
            end do
        end do

    end subroutine undex_dp_1

    subroutine undex_dp_2(ntotal, ndummy, P)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy
        type(Particle), intent(inout) :: P(:)
        real(8), save :: boundary(2, 2)
        integer :: scale_k, n_dp_1
        logical, save :: first_entry = .true.
        integer i

        scale_k = 0
        boundary = reshape([-huge(0._8), -huge(0._8),  &
                             huge(0._8),  huge(0._8)], &
                             shape(boundary))
        if ( first_entry ) then
            do i = 1, ntotal+ndummy
                if ( P(i)%x(1) > boundary(1,1) ) then
                    boundary(1,1) = P(i)%x(1)
                end if
                if ( p(i)%x(2) > boundary(2,1) ) then
                    boundary(2,1) = P(i)%x(2)
                end if
                if ( p(i)%x(1) < boundary(1,2) ) then
                    boundary(1,2) = P(i)%x(1)
                end if
                if ( p(i)%x(2) < boundary(2,2) ) then
                    boundary(2,2) = P(i)%x(2)
                end if
            end do
            first_entry = .false.
        end if

        n_dp_1 = ndummy

        select case (Config%skf)
        case (1)
            scale_k = 2
        case (2, 3)
            scale_k = 3
        end select

        do i = 1, ntotal
            !!! Dummy particle II Upside
            if ( boundary(2, 1) - P(i)%x(2) < P(i)%SmoothingLength * scale_k ) then
                ndummy = ndummy + 1
                P(ntotal+ndummy) = P(i)
                P(ntotal+ndummy)%x(2) =  2*boundary(2, 1) - P(i)%x(2)
            end if

            !!! Dummy particle II Downside
            if ( P(i)%x(2) - boundary(2, 2) < P(i)%SmoothingLength * scale_k ) then
                ndummy = ndummy + 1
                P(ntotal+ndummy) = P(i)
                P(ntotal+ndummy)%x(2) =  2*boundary(2, 2) - P(i)%x(2)
            end if

            !!! Monaghan type dummy particle on the Right side
            if ( boundary(1, 1) - P(i)%x(1) < P(i)%SmoothingLength * scale_k ) then
                ndummy = ndummy + 1
                P(ntotal+ndummy) = P(i)
                P(ntotal+ndummy)%x(1) =  2*boundary(1, 1) - P(i)%x(1)
            end if

            ! !!! Monaghan type dummy particle on the Left side
            ! if ( P(i)%x(1) - boundary(1, 2) < P(i)%SmoothingLength * scale_k ) then
            !     ndummy = ndummy + 1
            !     P(ntotal+ndummy) = P(i)
            !     P(ntotal+ndummy)%x(1) =  2*boundary(1, 2) - P(i)%x(1)
            ! end if
        end do

    end subroutine undex_dp_2

    subroutine dam_break_dp(ntotal, ndummy, P)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy
        type(Particle), intent(inout) :: P(:)
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
            dx = abs(P(2)%x(1) - P(2)%x(2))
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
                P(index)%x(:) = [wall_domain(1) - (l-0.5)*dx, &
                                 wall_domain(3) + (i-layer-0.5)*dx]
                P(index)%v(:) = 0
                P(index)%Density         = 1000
                P(index)%Mass            = P(index)%Density * dx * dx
                P(index)%Pressure        = 0
                P(index)%InternalEnergy  = 0
                P(index)%Type           = -P(1)%Type
                P(index)%SmoothingLength = dx
            end do

            !!! Dummy particle I on the Right side
            do i = 1, Ny
                ndummy = ndummy + 1
                index = ntotal + ndummy
                P(index)%x(:) = [wall_domain(2) + (l-0.5)*dx, &
                                 wall_domain(3) + (i-layer-0.5)*dx]
                P(index)%v(:) = 0
                P(index)%Density         = 1000
                P(index)%Mass            = P(index)%Density * dx * dx
                P(index)%Pressure        = 0
                P(index)%InternalEnergy  = 0
                P(index)%Type           = -P(1)%Type
                P(index)%SmoothingLength = dx
            end do

            !!! Dummy particle I on the Bottom
            Nx = floor((wall_domain(2) - wall_domain(1))/dx) + 1
            do i = 1, Nx
                ndummy = ndummy + 1
                index = ntotal + ndummy
                P(index)%x(:) = [wall_domain(1) + (i-0.5)*dx, &
                                 wall_domain(3) - (l-0.5)*dx]
                P(index)%v(:) = 0
                P(index)%Density         = 1000
                P(index)%Mass            = P(index)%Density * dx * dx
                P(index)%Pressure        = 0
                P(index)%InternalEnergy  = 0
                P(index)%Type            = -P(1)%Type
                P(index)%SmoothingLength = dx
            end do
        end do

        ! !!! Dummy particle I for baffle
        ! Nx = floor(3 / (dx / 2))
        ! do l = 1, layer
        !     do i = 1, Nx - l + 1
        !         ndummy = ndummy + 1
        !         index = ntotal + ndummy
        !         P(index)%x(:)  = [sum(wall_domain(1:2))/2 + (i+l-1.5)*dx/2, &
        !                           wall_domain(3) + (i-0.5)*dx/2]
        !         P(index)%v(:)  = 0
        !         P(index)%Density         = 1000
        !         P(index)%Mass            = P(index)%Density * dx * dx / 4
        !         P(index)%Pressure        = 0
        !         P(index)%InternalEnergy  = 0
        !         P(index)%Type            = -P(1)%Type
        !         P(index)%SmoothingLength = dx
        !     end do
        ! end do

        !!! Dummy particle for Lid
        Nx = floor((wall_domain(2) - wall_domain(1))/dx) + 1
        do l = 1, layer
            do i = 1, Nx
                ndummy = ndummy + 1
                index = ntotal + ndummy
                P(index)%x(:)  = [wall_domain(1) + (i-0.5)*dx, &
                                  wall_domain(4) - (l+0.5)*dx]
                P(index)%v(:)  = 0
                P(index)%Density         = 1000
                P(index)%Mass            = P(index)%Density * dx * dx
                P(index)%Pressure        = 0
                P(index)%InternalEnergy  = 0
                P(index)%Type            = -P(1)%Type
                P(index)%SmoothingLength = dx
            end do
        end do

    end subroutine dam_break_dp

    ! subroutine taylor_rod_dp_1(ntotal, ndummy, P)
    !     integer, intent(in) :: ntotal
    !     integer, intent(inout) :: ndummy
    !     type(Particle), intent(inout) :: P(:)
    !     real(8) :: dx
    !     save dx
    !     real(8) :: wall_domain(4)
    !     integer :: Nx, Ny
    !     integer :: layer
    !     logical :: first_entry = .true.
    !     save first_entry

    !     integer :: index
    !     integer i, j

    !     ndummy = 0

    !     if ( first_entry ) then
    !         dx = abs(P(2)%x(1) - P(2)%x(2)) / 2
    !         first_entry = .false.
    !     end if

    !     layer = 4
    !     wall_domain = [-1.14, 1.9, 0., 0.]  * 1e-2
    !     wall_domain(3) = - layer * dx

    !     Nx = floor((wall_domain(2) - wall_domain(1))/dx)
    !     Ny = floor((wall_domain(4) - wall_domain(3))/dx)
    !     do i = 1, Nx
    !         do j = 1, Ny
    !             ndummy = ndummy + 1
    !             index = ntotal + ndummy
    !             P(index)%x(:)  = [wall_domain(2) - (i-0.5)*dx, &
    !                             wall_domain(4) - (j-0.5)*dx]
    !             P(index)%v(:)  = 0
    !             P(index)%Density         = 7850
    !             P(index)%Mass            = P(index)%Density * dx * dx
    !             P(index)%Pressure        = 0
    !             P(index)%InternalEnergy  = 0
    !             P(index)%SoundSpeed      = 5000
    !             P(index)%Type           = -P(1)%Type
    !             P(index)%SmoothingLength = dx * 2
    !         end do
    !     end do

    ! end subroutine taylor_rod_dp_1

    subroutine taylor_rod_dp_2(ntotal, ndummy, P)
        integer, intent(in) :: ntotal
        integer, intent(inout) :: ndummy
        type(Particle), intent(inout) :: P(:)
        integer :: scale_k

        integer i

        scale_k = 0
        select case (Config%skf)
        case (1)
            scale_k = 2
        case (2, 3)
            scale_k = 3
        end select

        do i = 1, ntotal
            if ( P(i)%x(2) - 0 < P(i)%SmoothingLength * scale_k ) then
                ndummy = ndummy + 1
                P(ntotal+ndummy) = P(i)
                P(ntotal+ndummy)%x(2)  = -P(i)%x(2)
                P(ntotal+ndummy)%v(2)  = -P(i)%v(2)
                P(ntotal+ndummy)%Type = -P(i)%Type
            end if
        end do

    end subroutine taylor_rod_dp_2

end module dummy_part_m