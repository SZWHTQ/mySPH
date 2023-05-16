module input_m
    use ctrl_dict, only: Project, Config, Field
    use sph
    use output_m,  only: output
    use tools_m
    implicit none

contains
    subroutine input(ntotal, Particles)
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: Particles(:)

        if ( Config%read_ini_file_w ) then
            call read_initial_file(ntotal, Particles)
        else
            call write_initial_file(ntotal, Particles)
        end if

    end subroutine input

    !!! Load initial particle information from external disk file
    subroutine read_initial_file(ntotal, Particles)
        integer,        intent(inout) :: ntotal
        type(Particle), intent(inout) :: Particles(:)
        integer i, null


        open(11, file = Project%in_path // "ini.dat")

        write(*,*) "   Loading initial particle configuration..."
        read(11,*) ntotal
        write(*,*) "   Total number of particle: ", to_string(ntotal)
        write(*,*)

        do i = 1, ntotal
            read(11, "(I5, DT)", advance="no") null, Particles
            if ( null /= i ) then
                call print_warning(null, "Unconsistent particle index", BS)
            end if
        end do


    end subroutine read_initial_file

    subroutine write_initial_file(ntotal, Particles)
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: Particles(:)
        integer :: ini = 0

        select case(Project%nick)
        case("shear_cavity")
            call shear_cavity(ntotal, Particles)
        case("shock_tube")
            ntotal = Field%ntotal
            call shock_tube(ntotal, Particles)
        case("tnt_bar")
            ntotal = Field%ntotal
            call tnt_bar(ntotal, Particles)
        case("tnt_cylinder")
            call tnt_cylinder(ntotal, Particles)
        case("undex_cylinder")
            call undex_cylinder(ntotal, Particles)
        case("undex_chamber")
            call undex_chamber(ntotal, Particles)
        case("UNDEX")
            call undex(ntotal, Particles)
        case("dam_break")
            call damBreak(ntotal, Particles)
        case("db_gate")
            call damBreakwithElasticGate(ntotal, Particles)
        case("water_impact")
            call waterImpact(ntotal, Particles)
        case ("taylor_rod")
            call taylor_rod(ntotal, Particles)
        case ("can_beam")
            call cantilever_beam(ntotal, Particles)
        case ("beam_oil")
            write(*, "(A)") "!!! Unfinished Project !!!"
            call clammpedBeamwithOil(ntotal, Particles)
        end select

        call output(ini, Particles(1:ntotal))

        write(*,*) "Initial particle configuration generated"
        write(*,*) ">> Total number of particles: ", to_string(ntotal)
        write(*,*)

    end subroutine write_initial_file

    subroutine shock_tube(ntotal, P)
        integer, intent(in) :: ntotal
        type(Particle), intent(inout) :: P(:)
        real(8) :: dx
        real(8) :: gamma

        integer :: n_left
        integer i

        if ( Field%Dim /= 1 ) Field%Dim = 1

        gamma = 1.4
        dx = 0.6 / ntotal * 5
        n_left  = ntotal / 5 * 4


        do i = 1, n_left
            P(i)%x(1)            = -0.6 + dx / 4 * (i-1)
            P(i)%v(:)            = 0
            P(i)%Type            = 1
            P(i)%SmoothingLength = dx * 1.5
            P(i)%Density         = 1
            P(i)%Mass            = P(i)%Density * dx / 4
            P(i)%Pressure        = 1
            P(i)%InternalEnergy  = 2.5
            P(i)%SoundSpeed      = sqrt(gamma * P(i)%Pressure / P(i)%Density)
        end do

        do i = n_left+1, ntotal
            P(i)%x(1)            = 0 + dx * (i-n_left)
            P(i)%v(:)            = 0
            P(i)%Type            = 1
            P(i)%SmoothingLength = dx * 2
            P(i)%Density         = 0.25
            P(i)%Mass            = P(i)%Density * dx
            P(i)%Pressure        = 0.1795
            P(i)%InternalEnergy  = 1.795
            P(i)%SoundSpeed      = sqrt(gamma * P(i)%Pressure / P(i)%Density)
        end do

    end subroutine shock_tube

    subroutine shear_cavity(ntotal, P)
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: P(:)
        integer :: m, n, mp, np, k
        real(8) :: xl, yl, dx, dy

        integer i, j

        if ( Field%Dim /= 2 ) Field%Dim = 2

        m = 41
        n = 41
        mp = m-1
        np = n-1
        ntotal = mp * np
        xl = 1.e-3
        yl = 1.e-3
        dx = xl / mp
        dy = yl / np

        do i = 1, mp
            do j = 1, np
                k = j + (i-1)*np
                P(k)%x(1) = (i-1)*dx + dx/2._8
                P(k)%x(2) = (j-1)*dy + dy/2._8
            end do

        end do

        do i = 1, mp*np
            P(i)%v(1)            = 0
            P(i)%v(2)            = 0
            P(i)%Density         = 1000
            P(i)%Mass            = dx * dy * P(i)%Density
            P(i)%Pressure        = 0
            P(i)%InternalEnergy  = 357.1
            P(i)%Type           = 3
            P(i)%SmoothingLength = dx
        end do


    end subroutine shear_cavity

    subroutine tnt_bar(ntotal, P)
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: P(:)
        real(8) :: space_x
        real(8) :: rho0 = 1630, E0 = 4.29e6
        integer i

        if ( Field%Dim /= 1 ) Field%Dim = 1

        space_x = 0.1 / ntotal

        do i = 1, ntotal
            P(i)%x(:)            = 0 + space_x * i
            P(i)%v(:)            = 0
            P(i)%Density         = rho0
            P(i)%Mass            = P(i)%Density * space_x
            P(i)%Pressure        = 0
            P(i)%InternalEnergy  = E0
            P(i)%Type           = 0
            P(i)%SmoothingLength = space_x * 1.5
        end do

    end subroutine tnt_bar

    subroutine tnt_cylinder(ntotal, P)
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: P(:)
        real(8) :: rho0 = 1630, e0 = 4.29e6
        integer :: nr = 20, nt = 60
        real(8) :: dr
        real(8) :: theta, rot(2, 2)
        integer i, j, k

        if ( Field%Dim /= 2 ) Field%Dim = 2

        ntotal = nr * nt

        dr = 0.1 / nr
        theta = 2 * PI / nt
        rot = reshape( &
            [ cos(theta), sin(theta),   &
             -sin(theta), cos(theta) ], &
             shape(rot))

        do i = 1, nr
            P(i)%x(1)            = 0 + dr * i
            P(i)%x(2)            = 0
            P(i)%v(:)            = 0
            P(i)%Density         = rho0
            P(i)%Mass            = P(i)%Density * dr * ((2*PI*i*dr)/nt)
            P(i)%Pressure        = 0
            P(i)%InternalEnergy  = e0
            P(i)%Type           = 4
            P(i)%SmoothingLength = dr * 2
        end do

        do i = 1, nt - 1
            do j = 1, nr
                k = i * nr + j
                P(k)%x(:)            = matmul(rot, P(k-nr)%x)
                P(k)%v(:)            = 0
                P(k)%Density         = rho0
                P(k)%Mass            = P(k-nr)%Mass
                P(k)%Pressure        = 0
                P(k)%InternalEnergy  = e0
                P(k)%Type           = 4
                P(k)%SmoothingLength = dr * 2
            end do
        end do

    end subroutine tnt_cylinder

    subroutine undex_cylinder(ntotal, P)
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: P(:)
        real(8) :: rho0(2) = [1630., 1000.], e0(2) = [4.29e6, 0.]
        real(8) :: p0 = 101.325e3
        integer :: nr(2) = [10,  40], nt = 60
        real(8) :: r(2)  = [0.1, 0.4]
        real(8) :: dr
        real(8) :: theta, rot(2, 2)
        integer i, j, k, index

        if ( Field%Dim /= 2 ) Field%Dim = 2

        ntotal = sum(nr) * nt

        dr = sum(r) / sum(nr)
        theta = 2 * PI / nt
        rot = reshape( &
            [ cos(theta), sin(theta),   &
             -sin(theta), cos(theta) ], &
             shape(rot))

        do i = 1, nr(1)
            P(i)%x(1)            = 0 + dr * i
            P(i)%x(2)            = 0
            P(i)%v(:)            = 0
            P(i)%Density         = rho0(1)
            P(i)%mass            = P(i)%Density * dr * (theta*i*dr)
            P(i)%Pressure        = 0
            P(i)%InternalEnergy  = e0(1)
            P(i)%Type           = 0
            P(i)%SmoothingLength = dr * 2
        end do

        do i = nr(1) + 1, nr(1) + nr(2)
            P(i)%x(1)            = 0 + dr * i
            P(i)%x(2)            = 0
            P(i)%v(:)            = 0
            P(i)%Density         = rho0(2)
            P(i)%Mass            = P(i)%Density * dr * (theta*i*dr)
            P(i)%Pressure        = p0
            P(i)%InternalEnergy  = e0(2)
            P(i)%Type           = 6
            P(i)%SmoothingLength = dr * 2
        end do

        do i = 1, nt - 1
            do j = 1, sum(nr)
                k = i * sum(nr) + j
                index = k - sum(nr)
                P(k)%x(:)            = matmul(rot, P(index)%x(:))
                P(k)%v(:)            = 0
                P(k)%Density         = P(index)%Density
                P(k)%Mass            = P(index)%Mass
                P(k)%Pressure        = P(index)%Pressure
                P(k)%InternalEnergy  = P(index)%InternalEnergy
                P(k)%Type           = P(index)%Type
                P(k)%SmoothingLength = dr * 2
            end do
        end do

    end subroutine undex_cylinder

    subroutine undex_chamber(ntotal, P)
        use eos_m, only: mie_gruneisen_eos_of_water, jwl_eos
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: P(:)
        integer :: n(2),      tnt_n(2)
        real(8) :: origin(2), tnt_origin(2) !! Bottom left corner
        real(8) :: length(2), tnt_length(2)
        real(8) :: delta(2),  tnt_delta(2)
        ! integer :: layer, boundaryType
        ! integer :: l

        integer i, j, k


        if ( Field%Dim /= 2 )  Field%Dim = 2

        origin     = [real(8) :: -0.500, -0.500]
        length     = [real(8) ::  1.000,  1.000]
        tnt_origin = [real(8) :: -0.050, -0.050]
        tnt_length = [real(8) ::  0.100,  0.100]
        delta      = [real(8) ::  0.005,  0.005]
        tnt_delta  = delta / 3

        n     = int((length - tnt_length) / 2 / delta)
        tnt_n = int(tnt_length / tnt_delta) + 1

        k = 0
        ntotal = 0
        !!! Water to the left of the TNT
        do i = 1, n(1)
            do j = 1, n(2)*2 + int(tnt_length(2)/delta(2)+1)
                k = ntotal + (i-1) * (n(2)*2 + int(tnt_length(2)/delta(2)+1)) + j
                P(k)%x(:)            = origin + [i-1, j-1] * delta
                P(k)%Density         = 1000
                P(k)%Mass            = P(k)%Density * delta(1)*delta(2)
                P(k)%InternalEnergy  = 1e8
                P(k)%Type            = 6
                P(k)%SmoothingLength = 1.5 * sum(delta)/2
                call mie_gruneisen_eos_of_water(P(k)%Density, P(k)%InternalEnergy, &
                                                P(k)%Pressure, P(k)%SoundSpeed)
                if ( (P(k)%x(1) <= -0.45) .or. (abs(P(K)%x(2)) > 0.45) ) then
                    P(k)%Boundary = 2
                end if
            end do
        end do
        ntotal = k

        !!! Water below the TNT
        do i = 1, int(tnt_length(1)/delta(1)+1)
            do j = 1, n(2)
                k = ntotal + (i-1) * n(2) + j
                P(k)%x(:)            = [tnt_origin(1), origin(2)] + [i-1, j-1] * delta
                P(k)%Density         = 1000
                P(k)%Mass            = P(k)%Density * delta(1)*delta(2)
                P(k)%InternalEnergy  = 1e8
                P(k)%Type            = 6
                P(k)%SmoothingLength = 1.5 * sum(delta)/2
                call mie_gruneisen_eos_of_water(P(k)%Density, P(k)%InternalEnergy, &
                                                P(k)%Pressure, P(k)%SoundSpeed)
                if ( (P(K)%x(2) <= -0.45) ) then
                    P(k)%Boundary = 2
                end if
            end do
        end do
        ntotal = k

        !!! Water above the TNT
        do i = 1, int(tnt_length(1)/delta(1)+1)
            do j = 1, n(2)
                k = ntotal + (i-1) * n(2) + j
                P(k)%x(:)            = [tnt_origin(1), origin(2)] &
                                        + [0._8, length(2)] + [i-1, 1-j] * delta
                P(k)%Density         = 1000
                P(k)%Mass            = P(k)%Density * delta(1)*delta(2)
                P(k)%InternalEnergy  = 1e8
                P(k)%Type            = 6
                P(k)%SmoothingLength = 1.5 * sum(delta)/2
                call mie_gruneisen_eos_of_water(P(k)%Density, P(k)%InternalEnergy, &
                                                P(k)%Pressure, P(k)%SoundSpeed)
                if ( (P(K)%x(2) >= 0.45) ) then
                    P(k)%Boundary = 2
                end if
            end do
        end do
        ntotal = k

        !!! Water to the right of the TNT
        do i = 1, n(1)
            do j = 1, n(2)*2 + int(tnt_length(2)/delta(2)+1)
                k = ntotal + (i-1) * (n(2)*2 + int(tnt_length(2)/delta(2)+1)) + j
                P(k)%x(:)            = [tnt_origin(1)+tnt_length(1)+delta(1), origin(2)] &
                                        + [i-1, j-1] * delta
                P(k)%Density         = 1000
                P(k)%Mass            = P(k)%Density * delta(1)*delta(2)
                P(k)%InternalEnergy  = 1e8
                P(k)%Type            = 6
                P(k)%SmoothingLength = 1.5 * sum(delta)/2
                call mie_gruneisen_eos_of_water(P(k)%Density, P(k)%InternalEnergy, &
                                                P(k)%Pressure, P(k)%SoundSpeed)
                if ( (P(k)%x(1) >= 0.45) .or. (abs(P(K)%x(2)) > 0.45) ) then
                    P(k)%Boundary = 2
                end if
            end do
        end do
        ntotal = k

        !!! TNT
        k = 0
        do i = 1, tnt_n(1)
            do j = 1, tnt_n(2)
                k = ntotal + (i-1) * tnt_n(2) + j
                P(k)%x(:)             = tnt_origin + [i-1, j-1] * tnt_delta
                P(k)%Density          = 1630
                P(k)%Mass             = P(k)%Density * tnt_delta(1)*tnt_delta(2)
                P(k)%InternalEnergy   = 4.29e6
                P(k)%Type             = 5
                P(k)%SmoothingLength  = 1.5 * sum(tnt_delta)/2
                call jwl_eos(P(k)%Density, P(k)%InternalEnergy, P(k)%Pressure)
            end do
        end do
        ntotal = k

        n = int(length / delta) + 1

        ! !!! Sponge Layer
        ! layer = 20
        ! boundaryType = 2
        ! do l = 1, layer
        !     ! if ( l == layer ) then
        !     !     boundaryType = 1
        !     ! else
        !     !     boundaryType = 2
        !     ! end if
        !     !!! Upper side
        !     do i = 1, n(1) + 2*l
        !         k = ntotal + i
        !         P(k)%x(1) = origin(1) + delta(1) * (i-l-1)
        !         P(k)%x(2) = length(2) + origin(2) + delta(2) * l
        !         P(k)%v(:)            = 0
        !         P(k)%Density         = 1000
        !         P(k)%Mass            = P(k)%Density * delta(1)*delta(2)
        !         P(k)%InternalEnergy  = 1e8
        !         P(k)%Type            = 6
        !         P(k)%SmoothingLength = 1.5 * sum(delta)/2
        !         P(k)%Boundary        = boundaryType
        !         call mie_gruneisen_eos_of_water(P(k)%Density, P(k)%InternalEnergy, &
        !                                         P(k)%Pressure, P(k)%SoundSpeed)
        !     end do
        !     ntotal = k

        !     !!! Lower side
        !     do i = 1, n(1) + 2*l
        !         k = ntotal + i
        !         P(k)%x(1) = origin(1) + delta(1) * (i-l-1)
        !         P(k)%x(2) = origin(2) - delta(2) * l
        !         P(k)%v(:)            = 0
        !         P(k)%Density         = 1000
        !         P(k)%Mass            = P(k)%Density * delta(1)*delta(2)
        !         P(k)%InternalEnergy  = 1e8
        !         P(k)%Type            = 6
        !         P(k)%SmoothingLength = 1.5 * sum(delta)/2
        !         P(k)%Boundary        = boundaryType
        !         call mie_gruneisen_eos_of_water(P(k)%Density, P(k)%InternalEnergy, &
        !                                         P(k)%Pressure, P(k)%SoundSpeed)
        !     end do
        !     ntotal = k

        !     !!! Left side
        !     do i = 1, n(2) + 2*l
        !         k = ntotal + i
        !         P(k)%x(1) = origin(1) - delta(1) * l
        !         P(k)%x(2) = origin(2) + delta(2) * (i-l-1)
        !         P(k)%v(:)            = 0
        !         P(k)%Density         = 1000
        !         P(k)%Mass            = P(k)%Density * delta(1)*delta(2)
        !         P(k)%InternalEnergy  = 1e8
        !         P(k)%Type            = 6
        !         P(k)%SmoothingLength = 1.5 * sum(delta)/2
        !         P(k)%Boundary        = boundaryType
        !         call mie_gruneisen_eos_of_water(P(k)%Density, P(k)%InternalEnergy, &
        !                                         P(k)%Pressure, P(k)%SoundSpeed)
        !     end do
        !     ntotal = k

        !     !!! Right side
        !     do i = 1, n(2) + 2*l
        !         k = ntotal + i
        !         P(k)%x(1) = length(1) + origin(1) + delta(1) * l
        !         P(k)%x(2) = origin(2) + delta(2) * (i-l-1)
        !         P(k)%v(:)            = 0
        !         P(k)%Density         = 1000
        !         P(k)%Mass            = P(k)%Density * delta(1)*delta(2)
        !         P(k)%InternalEnergy  = 1e8
        !         P(k)%Type            = 6
        !         P(k)%SmoothingLength = 1.5 * sum(delta)/2
        !         P(k)%Boundary        = boundaryType
        !         call mie_gruneisen_eos_of_water(P(k)%Density, P(k)%InternalEnergy, &
        !                                         P(k)%Pressure, P(k)%SoundSpeed)
        !     end do
        !     ntotal = k
        ! end do

        ! !!! Solid
        ! layer = 6
        ! do l = 1, layer
        !     !!! Upper side
        !     do i = 1, n(1) + 2*(l+1)
        !         k = ntotal + i
        !         P(k)%x(1) = origin(1) + delta(1) * (i-l-2)
        !         P(k)%x(2) = length(2)+origin(2) + delta(2) * (l+1)
        !         P(k)%v(:)            = 0
        !         P(k)%Density         = 7850
        !         P(k)%Mass            = P(k)%Density * delta(1)*delta(2)
        !         P(k)%Pressure        = 0
        !         P(k)%InternalEnergy  = 0
        !         P(k)%SoundSpeed      = 5000
        !         P(k)%Type            = 8
        !         P(k)%SmoothingLength = 1.5 * sum(delta)/2
        !     end do
        !     ntotal = k

        !     !!! Lower side
        !     do i = 1, n(1) + 2*(l+1)
        !         k = ntotal + i
        !         P(k)%x(1) = origin(1) + delta(1) * (i-l-2)
        !         P(k)%x(2) = origin(2) - delta(2) * (l+1)
        !         P(k)%v(:)            = 0
        !         P(k)%Density         = 7850
        !         P(k)%Mass            = P(k)%Density * delta(1)*delta(2)
        !         P(k)%Pressure        = 0
        !         P(k)%InternalEnergy  = 0
        !         P(k)%SoundSpeed      = 5000
        !         P(k)%Type            = 8
        !         P(k)%SmoothingLength = 1.5 * sum(delta)/2
        !     end do
        !     ntotal = k

        !     !!! Left side
        !     do i = 1, n(2) + 2*l
        !         k = ntotal + i
        !         P(k)%x(1) = origin(1) - delta(1) * (l+1)
        !         P(k)%x(2) = origin(2) + delta(2) * (i-l-1)
        !         P(k)%v(:)            = 0
        !         P(k)%Density         = 7850
        !         P(k)%Mass            = P(k)%Density * delta(1)*delta(2)
        !         P(k)%Pressure        = 0
        !         P(k)%InternalEnergy  = 0
        !         P(k)%SoundSpeed      = 5000
        !         P(k)%Type            = 8
        !         P(k)%SmoothingLength = 1.5 * sum(delta)/2
        !     end do
        !     ntotal = k

        !     !!! Right side
        !     do i = 1, n(2) + 2*l
        !         k = ntotal + i
        !         P(k)%x(1) = length(1) + origin(1) + delta(1) * (l+1)
        !         P(k)%x(2) = origin(2) + delta(2) * (i-l-1)
        !         P(k)%v(:)            = 0
        !         P(k)%Density         = 7850
        !         P(k)%Mass            = P(k)%Density * delta(1)*delta(2)
        !         P(k)%Pressure        = 0
        !         P(k)%InternalEnergy  = 0
        !         P(k)%SoundSpeed      = 5000
        !         P(k)%Type            = 8
        !         P(k)%SmoothingLength = 1.5 * sum(delta)/2
        !     end do
        !     ntotal = k
        ! end do

    end subroutine undex_chamber

    subroutine undex(ntotal, P)
        use geometry_m
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: P(:)
        real(8) :: fluid_domian(4) = [0.0, 0.1, -0.045, 0.045]
        real(8) :: confined = 1
        type(rectangle_t) :: confined_domian
        real(8) :: dx = 5e-4
        real(8) :: tnt_center(2) = [0.05, 0.0]
        real(8) :: tnt_radius = 0.001397433858052
        type(circle_t) :: tnt_domain
        integer :: Nx, Ny

        integer i, j

        ntotal = 0
        confined_domian = rectangle_t(tnt_center, [4, 4]*tnt_radius, 0)
        tnt_domain = circle_t(tnt_center, tnt_radius, 0)
        Nx = floor((fluid_domian(2) - fluid_domian(1)) / dx) + 1
        Ny = floor((fluid_domian(4) - fluid_domian(3)) / dx) + 1
        do i = 1, Nx
            do j = 1, Ny
                ntotal = ntotal + 1
                P(ntotal)%x(:)  = [fluid_domian(1) + i * dx, &
                                 fluid_domian(3) + j * dx]
                if ( confined_domian%contain(point_t(P(ntotal)%x(:), 0)) ) then
                    ntotal = ntotal - 1
                    cycle
                end if
                P(ntotal)%Density   = 1000
                P(ntotal)%Pressure     = 0
                P(ntotal)%InternalEnergy     = 2e7
                P(ntotal)%Type = 6
                P(ntotal)%Mass  = P(ntotal)%Density * dx*dx
                P(ntotal)%SmoothingLength  = dx
            end do
        end do

        Nx = floor(confined_domian%length(1) / (dx/confined))
        Ny = floor(confined_domian%length(2) / (dx/confined))
        do i = 1, Nx
            do j = 1, Ny
                ntotal = ntotal + 1
                P(ntotal)%x(:) = confined_domian%center   &
                             - confined_domian%length/2 &
                             + [i-0.3, j-0.3] * (dx/confined)
                ! if ( tnt_domain%contain(point_t(x(:,ntotal), 0)) ) then
                    !!! TNT
                    P(ntotal)%Density   = 1630
                    P(ntotal)%Pressure     = 0
                    P(ntotal)%InternalEnergy     = 4.29e6
                    P(ntotal)%Type = 5
                ! else
                !     ! ntotal = ntotal - 1
                !     !!! Water
                !     P(ntotal)%Density   = 1000
                !     P(ntotal)%Pressure     = 0
                !     P(ntotal)%InternalEnergy     = 0
                !     P(ntotal)%Type = 6
                ! end if
                P(ntotal)%Mass = P(ntotal)%Density * dx*dx / (confined**2)
                P(ntotal)%SmoothingLength = dx / (confined**2)
            end do
        end do

        ! write(*,*) ntotal

    end subroutine undex

    subroutine damBreak(ntotal, P)
        use eos_m, only: arti_water_eos_1
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: P(:)
        ! real(8) :: fluid_domian(4) = [0, 2, 0, 1]
        real(8) :: fluid_domian(4) = [0, 25, 0, 25]
        ! real(8) :: dx = 0.02
        real(8) :: dx = 0.2!, bulk=2e7, rho0=1100, c
        integer :: Nx, Ny

        integer i, j, k

        Nx = int((fluid_domian(2) - fluid_domian(1)) / dx)
        Ny = int((fluid_domian(4) - fluid_domian(3)) / dx)
        do i = 1, Nx
            do j = 1, Ny
                k = (i-1) * Ny + j
                P(k)%x(:)            = [fluid_domian(1) + (i-0.5) * dx, &
                                        fluid_domian(3) + (j-0.5) * dx]
                P(K)%v(:)            = 0
                P(k)%Density         = 1000
                P(k)%Mass            = P(k)%Density * dx * dx
                P(k)%Pressure        = P(k)%Density * 9.81 * (fluid_domian(4) - P(k)%x(2))
                P(k)%InternalEnergy  = 0
                P(k)%SoundSpeed      = 50
                P(k)%Type            = 2
                P(k)%SmoothingLength = dx
                call arti_water_eos_1(P(k)%Density, P(k)%Pressure, Density=.true.)
            end do
        end do

        ntotal = Nx * Ny

        ! c = sqrt(bulk / rho0)
        ! Nx = int(0.6 / (dx))
        ! Ny = int(10  / (dx))
        ! do i = 1, Nx
        !     do j = 1, Ny
        !         k = ntotal + (i-1) * Ny + j
        !         P(k)%x(:)            = [50 + (i-(Nx/2.)-0.5) * (dx), &
        !                                 fluid_domian(3) + (j-0.5) * (dx)]
        !         P(K)%v(:)            = 0
        !         P(k)%Density         = rho0
        !         P(k)%Mass            = P(k)%Density * dx * dx / 4
        !         P(k)%Pressure        = 0
        !         P(k)%InternalEnergy  = 0
        !         P(k)%SoundSpeed      = c
        !         P(k)%Type            = 102
        !         P(k)%SmoothingLength = dx
        !         if ( j == 1 ) then
        !             P(k)%Boundary = 1
        !         end if
        !     end do
        ! end do

        ! ntotal = ntotal + Nx * Ny

    end subroutine damBreak

    subroutine damBreakwithElasticGate(ntotal, P)
        use eos_m, only: arti_water_eos_1
        use geometry_m, only: rectangle_t
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: P(:)
        type(rectangle_t) :: domain
        real(8) :: dx, h
        integer :: nx, ny

        integer i, j, k

        ntotal = 0

        dx = 0.001
        h = 0.14

        !!! Fluid
        domain = rectangle_t([0.05, 0.07], [0.10, 0.14], 0)
        nx = int((domain%length(1)) / dx)
        ny = int((domain%length(2)) / dx) + 1
        do i = 1, nx
            do j = 1, ny
                k = ntotal + (i-1) * ny + j
                P(k)%x(:)            = domain%center        &
                                        - domain%length / 2 &
                                        + [i-0.5, j-0.5] * dx
                P(K)%v(:)            = 0
                P(k)%Density         = 1000
                P(k)%Mass            = P(k)%Density * dx * dx
                P(k)%Pressure        = P(k)%Density * 9.81 * max((h - P(k)%x(2)), 0._8)
                P(k)%InternalEnergy  = 0
                P(k)%SoundSpeed      = 50
                P(k)%Type            = 2
                P(k)%SmoothingLength = dx
                call arti_water_eos_1(P(k)%Density, P(k)%Pressure, Density=.true.)
            end do
        end do
        ntotal = ntotal + nx * ny

        !!! Solid
        domain = rectangle_t([0.1025, 0.0395], [0.005, 0.079], 0)
        nx = int((domain%length(1)) / dx) + 1
        ny = int((domain%length(2)) / dx) + 1
        do i = 1, nx
            do j = 1, ny
                k = ntotal + (i-1) * ny + j
                P(k)%x(:)            = domain%center        &
                                        - domain%length / 2 &
                                        + [i-0.5, j+0.5] * dx
                P(K)%v(:)            = 0
                P(k)%Density         = 1100
                P(k)%Mass            = P(k)%Density * dx * dx
                P(k)%Pressure        = 0
                P(k)%InternalEnergy  = 0
                P(k)%SoundSpeed      = 134.84
                P(k)%Type            = 102
                P(k)%SmoothingLength = dx
                if ( j == ny ) then
                    P(k)%Boundary = 1
                end if
            end do
        end do
        ntotal = ntotal + nx * ny

    end subroutine damBreakwithElasticGate

    subroutine waterImpact(ntotal, P)
        use eos_m, only: arti_water_eos_1
        use geometry_m, only: rectangle_t
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: P(:)
        type(rectangle_t) :: domain
        real(8) :: dx, h
        integer :: nx, ny

        integer i, j, k

        ntotal = 0

        dx = 0.002
        h = 0.292

        !!! Fluid
        domain = rectangle_t([0.073, 0.146], [0.146, 0.292], 0)
        nx = int((domain%length(1)) / dx) + 1
        ny = int((domain%length(2)) / dx) + 1
        do i = 1, nx
            do j = 1, ny
                k = ntotal + (i-1) * ny + j
                P(k)%x(:)            = domain%center        &
                                        - domain%length / 2 &
                                        + [i-0.5, j-0.5] * dx
                P(K)%v(:)            = 0
                P(k)%Density         = 1000
                P(k)%Mass            = P(k)%Density * dx * dx
                P(k)%Pressure        = P(k)%Density * 9.81 * max((h - P(k)%x(2)), 0._8)
                P(k)%InternalEnergy  = 0
                P(k)%SoundSpeed      = 50
                P(k)%Type            = 2
                P(k)%SmoothingLength = dx
                call arti_water_eos_1(P(k)%Density, P(k)%Pressure, Density=.true.)
            end do
        end do
        ntotal = ntotal + nx * ny

        !!! Solid
        domain = rectangle_t([0.292, 0.04], [0.012, 0.08], 0)
        nx = int((domain%length(1)) / dx) + 1
        ny = int((domain%length(2)) / dx) + 1
        do i = 1, nx
            do j = 1, ny
                k = ntotal + (i-1) * ny + j
                P(k)%x(:)            = domain%center        &
                                        - domain%length / 2 &
                                        + [i-0.5, j+0.5] * dx
                P(K)%v(:)            = 0
                P(k)%Density         = 2500
                P(k)%Mass            = P(k)%Density * dx * dx
                P(k)%Pressure        = 0
                P(k)%InternalEnergy  = 0
                P(k)%SoundSpeed      = 11.489
                P(k)%Type            = 103
                P(k)%SmoothingLength = dx
                if ( j == 1 ) then
                    P(k)%Boundary = 1
                end if
            end do
        end do
        ntotal = ntotal + nx * ny

    end subroutine waterImpact

    subroutine clammpedBeamwithOil(ntotal, P)
        use eos_m, only: oil_eos, arti_eos_of_102
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: P(:)
        real(8) :: domain(4)
        real(8) :: dx = 0.001, bulk=2e7, rho0=1100, c
        integer :: Nx, Ny

        integer i, j, k

        ntotal = 0
        c = sqrt(bulk / rho0)

        !!! Left side
        domain = [real(8) :: -304.5, -2, 0, 115] * 1e-3
        Nx = int((domain(2) - domain(1)) / dx) - 2
        Ny = int((domain(4) - domain(3)) / dx) + 1
        do i = 1, Nx
            do j = 1, Ny
                k = ntotal + (i-1) * Ny + j
                P(k)%x(:)            = [domain(1) + (i-0.5) * dx, &
                                        domain(3) + (j-0.5) * dx]
                P(K)%v(:)            = 0
                P(k)%Density         = 917
                P(k)%Mass            = P(k)%Density * dx * dx
                P(k)%InternalEnergy  = 0
                P(k)%SoundSpeed      = 50
                P(k)%Type            = 8
                P(k)%SmoothingLength = dx
                call oil_eos(P(k)%Density, P(k)%Pressure)
            end do
        end do
        ntotal = ntotal + Nx * Ny

        !!! Right side
        domain = [real(8) :: 2, 304.5, 0, 115] * 1e-3
        Nx = int((domain(2) - domain(1)) / dx) - 2
        Ny = int((domain(4) - domain(3)) / dx) + 1
        do i = 1, Nx
            do j = 1, Ny
                k = ntotal + (i-1) * Ny + j
                P(k)%x(:)            = [domain(1) + (i+2.0) * dx, &
                                        domain(3) + (j-0.5) * dx]
                P(K)%v(:)            = 0
                P(k)%Density         = 917
                P(k)%Mass            = P(k)%Density * dx * dx
                P(k)%InternalEnergy  = 0
                P(k)%Type            = 8
                P(k)%SmoothingLength = dx
                call oil_eos(P(k)%Density, P(k)%Pressure)
            end do
        end do
        ntotal = ntotal + Nx * Ny

        !!! Solid
        domain = [real(8) :: -2, 2, 0, 115] * 1e-3
        Nx = int((domain(2) - domain(1)) / dx)
        Ny = int((domain(4) - domain(3)) / dx)
        do i = 1, Nx
            do j = 1, Ny
                k = ntotal + (i-1) * Ny + j
                P(k)%x(:)            = [domain(1) + (i-0.5) * dx, &
                                        domain(3) + (j-0.5) * dx]
                P(K)%v(:)            = 0
                P(k)%Density         = rho0
                P(k)%Mass            = P(k)%Density * dx * dx
                P(k)%Pressure        = 0
                P(k)%InternalEnergy  = 0
                P(k)%SoundSpeed      = c
                P(k)%Type            = 102
                P(k)%SmoothingLength = dx
                if ( j == 1 ) then
                    P(k)%Boundary = 1
                end if
            end do
        end do
        ntotal = ntotal + Nx * Ny

    end subroutine clammpedBeamwithOil

    subroutine taylor_rod(ntotal, P)
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: P(:)
        real(8) :: solid_domain(4) = [0., 0.760, &
                                      0., 2.546] * 1e-2
        real(8) :: dx = 0.038 * 1e-2
        integer :: Nx, Ny

        integer i, j, k

        Nx = floor((solid_domain(2) - solid_domain(1)) / dx)
        Ny = floor((solid_domain(4) - solid_domain(3)) / dx)
        do j = 1, Ny
            do i = 1, Nx
                ! k = (i-1) * Ny + j
                k = (j-1) * Nx + i
                P(k)%x(:)            = [solid_domain(1) + (i-0.5) * dx, &
                                        solid_domain(3) + (j-0.5) * dx]
                P(k)%v(1)            = 0
                P(k)%v(2)            = -221
                P(k)%Density         = 7850
                P(k)%Mass            = P(k)%Density * dx * dx
                P(k)%Pressure        = 0
                P(k)%InternalEnergy  = 0
                P(k)%SoundSpeed      = 5000
                P(k)%Type            = 101
                P(k)%SmoothingLength = dx * 2
            end do
        end do

        ntotal = Nx * Ny

    end subroutine taylor_rod

    subroutine cantilever_beam(ntotal, P)
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: P(:)
        real(8) :: dx = 1e-3, bulk=2e7, rho0=1100, c
        integer :: nx = 101, ny = 5

        integer i, j, k

        c = sqrt(bulk / rho0)

        do i = 1, nx
            do j = 1, ny
                k = (i-1) * ny + j
                P(k)%x(:)            = [i-1, j-1] * dx
                P(k)%v(:)            = 0
                P(k)%Density         = rho0
                P(k)%Mass            = P(k)%Density * dx * dx
                P(k)%Pressure        = 0
                P(k)%InternalEnergy  = 0
                P(k)%SoundSpeed      = c
                P(k)%Type            = 102
                P(k)%SmoothingLength = dx * 2
                if ( i == 1 ) then
                    P(k)%Boundary = 1
                end if
            end do
        end do

        ntotal = nx * ny

    end subroutine cantilever_beam

end module input_m