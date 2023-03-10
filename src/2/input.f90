module input_m
    use ctrl_dict, only: ntotal, dim, &
                         read_ini_file_W
    use initial_m, only: itype, x, v, mass, rho, p, e, c, hsml
    use parse_toml_m, only: in_path, out_path, vtk_path, nick
    use output_m,  only: output
    use tools_m
    implicit none

contains
    subroutine input()

        if ( read_ini_file_w ) then
            call read_initial_file()
        else
            call write_initial_file()
        end if

    end subroutine input

    !!! Load initial particle information from external disk file
    subroutine read_initial_file()
        integer i, null


        open(11, file = in_path // "ini.dat")

        write(*,*) "   Loading initial particle configuration..."
        read(11,*) ntotal
        write(*,*) "   Total number of particle: ", to_string(ntotal)
        write(*,*)

        do i = 1, ntotal
            read(11, *) null    , x(:, i), v(:, i), &
                        mass(i) , rho(i) , p(i)   , e(i), &
                        itype(i), hsml(i)
            if ( null /= i ) then
                call print_warning(null, "Unconsistent particle index", BS)
            end if
        end do


    end subroutine read_initial_file

    subroutine write_initial_file()
        use ctrl_dict, only: maxn
        real(8) :: div_r(maxn)
        real(8) :: Stress(dim,dim,maxn)
        integer :: ini = 0

        select case(nick)
        case("shear_cavity")
            call shear_cavity()
        case("shock_tube")
            call shock_tube()
        case("tnt_bar")
            call tnt_bar()
        case("tnt_cylinder")
            call tnt_cylinder()
        case("undex_cylinder")
            call undex_cylinder()
        case("undex_chamber")
            call undex_chamber()
        case("UNDEX")
            call undex()
        case("dam_break")
            call dam_break()
        case ("taylor_rod")
            call taylor_rod()
        end select

        call output(ini, ntotal, itype, x, v, mass, rho, p, e, c, hsml, div_r, Stress)

        write(*,*) "Initial particle configuration generated"
        write(*,*) ">> Total number of particles: ", to_string(ntotal)
        write(*,*)

    end subroutine write_initial_file

    subroutine shock_tube()
        real(8) :: dx
        real(8) :: gamma

        integer :: n_left
        integer i

        if ( dim /= 1 ) dim = 1

        gamma = 1.4
        dx = 0.6 / ntotal * 5
        n_left  = ntotal / 5 * 4


        do i = 1, n_left
            x(1, i)  = -0.6 + dx / 4 * (i-1)
            v(:, i)  = 0
            itype(i) = 1
            hsml(i)  = dx * 2
            rho(i)   = 1
            mass(i)  = rho(i) * dx / 4
            ! mass(i) = 0.75/ ntotal
            p(i)     = 1
            e(i)     = 2.5
            c(i)     = sqrt(gamma*p(i)/rho(i))
        end do

        do i = n_left+1, ntotal
            x(1, i)  = 0 + dx * (i-n_left)
            v(:, i)  = 0
            itype(i) = 1
            hsml(i)  = dx * 2
            rho(i)   = 0.25
            mass(i)  = rho(i) * dx
            ! mass(i) = 0.75/ ntotal
            p(i)     = 0.1795
            e(i)     = 1.795
            c(i)     = sqrt(gamma*p(i)/rho(i))
        end do

    end subroutine shock_tube

    subroutine shear_cavity()
        integer :: m, n, mp, np, k
        real(8) :: xl, yl, dx, dy

        integer i, j

        if ( dim /= 2 ) dim = 2

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
                x(1, k) = (i-1)*dx + dx/2._8
                x(2, k) = (j-1)*dy + dy/2._8
            end do

        end do

        do i = 1, mp*np
            v(1, i)  = 0
            v(2, i)  = 0
            rho(i)   = 1000
            mass(i)  = dx*dy*rho(i)
            p(i)     = 0
            e(i)     = 357.1
            itype(i) = 3
            hsml(i)  = dx
        end do


    end subroutine shear_cavity

    subroutine tnt_bar()
        real(8) :: space_x
        real(8) :: rho0 = 1630, E0 = 4.29e6
        integer i

        if ( dim /= 1 ) dim = 1

        space_x = 0.1 / ntotal

        do i = 1, ntotal
            x(:, i)  = 0 + space_x * i
            v(:, i)  = 0
            rho(i)   = rho0
            mass(i)  = rho(i) * space_x
            p(i)     = 0
            e(i)     = E0
            itype(i) = 0
            hsml(i)  = space_x * 1.5
        end do

    end subroutine tnt_bar

    subroutine tnt_cylinder()
        real(8) :: rho0 = 1630, e0 = 4.29e6
        integer :: nr = 20, nt = 60
        real(8) :: dr
        real(8) :: theta, rot(2, 2)
        integer i, j, k

        if ( dim /= 2 ) dim = 2

        ntotal = nr * nt

        dr = 0.1 / nr
        theta = 2 * PI / nt
        rot = reshape( &
            [ cos(theta), sin(theta),   &
             -sin(theta), cos(theta) ], &
             shape(rot))

        do i = 1, nr
            x(1, i)  = 0 + dr * i
            x(2, i)  = 0
            v(:, i)  = 0
            rho( i)  = rho0
            mass(i)  = rho(i) * dr * ((2*PI*i*dr)/nt)
            p(i)     = 0
            e(i)     = e0
            itype(i) = 4
            hsml( i) = dr * 2
        end do

        do i = 1, nt - 1
            do j = 1, nr
                k = i * nr + j
                x(:, k)  = matmul(rot, x(:, k-nr))
                v(:, k)  = 0
                rho( k)  = rho0
                mass(k)  = mass(k-nr)
                p(k)     = 0
                e(k)     = e0
                itype(k) = 4
                hsml( k) = dr * 2
            end do
        end do

    end subroutine tnt_cylinder

    subroutine undex_cylinder()
        real(8) :: rho0(2) = [1630., 1000.], e0(2) = [4.29e6, 0.]
        real(8) :: p0 = 101.325e3
        integer :: nr(2) = [10,  40], nt = 60
        real(8) :: r(2)  = [0.1, 0.4]
        real(8) :: dr
        real(8) :: theta, rot(2, 2)
        integer i, j, k

        if ( dim /= 2 ) dim = 2

        ntotal = sum(nr) * nt

        dr = sum(r) / sum(nr)
        theta = 2 * PI / nt
        rot = reshape( &
            [ cos(theta), sin(theta),   &
             -sin(theta), cos(theta) ], &
             shape(rot))

        do i = 1, nr(1)
            x(1, i)  = 0 + dr * i
            x(2, i)  = 0
            v(:, i)  = 0
            rho( i)  = rho0(1)
            mass(i)  = rho(i) * dr * (theta*i*dr)
            p(i)     = 0
            e(i)     = e0(1)
            itype(i) = 0
            hsml( i) = dr * 2
        end do

        do i = nr(1) + 1, nr(1) + nr(2)
            x(1, i)  = 0 + dr * i
            x(2, i)  = 0
            v(:, i)  = 0
            rho( i)  = rho0(2)
            mass(i)  = rho(i) * dr * (theta*i*dr)
            p(i)     = p0
            e(i)     = e0(2)
            itype(i) = 6
            hsml( i) = dr * 2
        end do

        do i = 1, nt - 1
            do j = 1, sum(nr)
                k = i * sum(nr) + j
                x(:, k)  = matmul(rot, x(:, k-sum(nr)))
                v(:, k)  = 0
                rho( k)  = rho(k-sum(nr))
                mass(k)  = mass(k-sum(nr))
                p(k)     = p(k-sum(nr))
                e(k)     = e(k-sum(nr))
                itype(k) = itype(k-sum(nr))
                hsml( k) = dr * 2
            end do
        end do

    end subroutine undex_cylinder

    subroutine undex_chamber()
        integer :: water_n(2) = [90, 90], &
                     tnt_n(2) = [61, 61]
        real(8) :: water_origin(2) = [-0.5, -0.5], &
                     tnt_origin(2) = [-0.05, -0.05] !! Bottom left corner
        real(8) :: length(2) = [1, 1], tnt_length(2) = [0.1, 0.1]
        real(8) :: delta(2), tnt_delta(2)
        integer i, j, k

        if ( dim /= 2 )  dim = 2

        delta = (length - tnt_length) / (water_n*2)
        tnt_delta = tnt_length / (tnt_n - 1)

        k = 0
        ntotal = 0
        !!! Water to the left of the TNT
        do i = 1, water_n(1)
            do j = 1, water_n(2)*2 + int(tnt_length(2)/delta(2)+1)
                k = ntotal + (i-1) * (water_n(2)*2 + int(tnt_length(2)/delta(2)+1)) + j
                rho(k)   = 1000
                mass(k)  = rho(k) * delta(1)*delta(2)
                p(k)     = 0
                e(k)     = 1e8
                itype(k) = 6
                hsml(k)  = 1.5 * sum(delta)/2
                x(:, k) = water_origin + [i-1, j-1] * delta
            end do
        end do
        ntotal = k

        !!! Water below the TNT
        do i = 1, int(tnt_length(1)/delta(1)+1)
            do j = 1, water_n(2)
                k = ntotal + (i-1) * water_n(2) + j
                x(:, k) = [tnt_origin(1), water_origin(2)] + [i-1, j-1] * delta
                rho(k)   = 1000
                mass(k)  = rho(k) * delta(1)*delta(2)
                p(k)     = 0
                e(k)     = 1e8
                itype(k) = 6
                hsml(k)  = 1.5 * sum(delta)/2
            end do
        end do
        ntotal = k

        !!! Water above the TNT
        do i = 1, int(tnt_length(1)/delta(1)+1)
            do j = 1, water_n(2)
                k = ntotal + (i-1) * water_n(2) + j
                x(:, k) = [tnt_origin(1), water_origin(2)] + [0._8, length(2)] + [i-1, 1-j] * delta
                rho(k)   = 1000
                mass(k)  = rho(k) * delta(1)*delta(2)
                p(k)     = 0
                e(k)     = 1e8
                itype(k) = 6
                hsml(k)  = 1.5 * sum(delta)/2
            end do
        end do
        ntotal = k

        !!! Water to the right of the TNT
        do i = 1, water_n(1)
            do j = 1, water_n(2)*2 + int(tnt_length(2)/delta(2)+1)
                k = ntotal + (i-1) * (water_n(2)*2 + int(tnt_length(2)/delta(2)+1)) + j
                x(:, k) = [tnt_origin(1)+tnt_length(1)+delta(1), water_origin(2)] + [i-1, j-1] * delta
                rho(k)   = 1000
                mass(k)  = rho(k) * delta(1)*delta(2)
                p(k)     = 0
                e(k)     = 1e8
                itype(k) = 6
                hsml(k)  = 1.5 * sum(delta)/2
            end do
        end do
        ntotal = k

        !!! TNT
        k = 0
        do i = 1, tnt_n(1)
            do j = 1, tnt_n(2)
                k = ntotal + (i-1) * tnt_n(2) + j
                x(:, k) = tnt_origin + [i-1, j-1] * tnt_delta
                rho(k)   = 1630
                mass(k)  = rho(k) * tnt_delta(1)*tnt_delta(2)
                p(k)     = 0
                e(k)     = 4.29e6
                itype(k) = 5
                hsml(k)  = 1.5 * sum(tnt_delta)/2
            end do
        end do
        ntotal = k

    end subroutine undex_chamber

    subroutine undex()
        use geometry_m
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
                x(:, ntotal)  = [fluid_domian(1) + i * dx, &
                                 fluid_domian(3) + j * dx]
                if ( confined_domian%contain(point_t(x(:,ntotal), 0)) ) then
                    ntotal = ntotal - 1
                    cycle
                end if
                rho(ntotal)   = 1000
                p(ntotal)     = 0
                e(ntotal)     = 2e7
                itype(ntotal) = 6
                mass(ntotal)  = rho(ntotal) * dx*dx 
                hsml(ntotal)  = dx
            end do
        end do

        Nx = floor(confined_domian%length(1) / (dx/confined))
        Ny = floor(confined_domian%length(2) / (dx/confined))
        do i = 1, Nx
            do j = 1, Ny
                ntotal = ntotal + 1
                x(:, ntotal) = confined_domian%center   &
                             - confined_domian%length/2 &
                             + [i-0.3, j-0.3] * (dx/confined)
                ! if ( tnt_domain%contain(point_t(x(:,ntotal), 0)) ) then
                    !!! TNT
                    rho(ntotal)   = 1630
                    p(ntotal)     = 0
                    e(ntotal)     = 4.29e6
                    itype(ntotal) = 5
                ! else
                !     ! ntotal = ntotal - 1
                !     !!! Water
                !     rho(ntotal)   = 1000
                !     p(ntotal)     = 0
                !     e(ntotal)     = 0
                !     itype(ntotal) = 6
                ! end if
                mass(ntotal) = rho(ntotal) * dx*dx / (confined**2)
                hsml(ntotal) = dx / (confined**2)
            end do
        end do

        ! write(*,*) ntotal

    end subroutine undex

    subroutine dam_break()
        ! real(8) :: fluid_domian(4) = [0, 2, 0, 1]
        real(8) :: fluid_domian(4) = [0, 25, 0, 25]
        ! real(8) :: dx = 0.02
        real(8) :: dx = 0.2
        integer :: Nx, Ny

        integer i, j, k

        Nx = int((fluid_domian(2) - fluid_domian(1)) / dx)
        Ny = int((fluid_domian(4) - fluid_domian(3)) / dx)
        do i = 1, Nx
            do j = 1, Ny
                k = (i-1) * Ny + j
                x(:, k)  = [fluid_domian(1) + (i-0.5) * dx, &
                            fluid_domian(3) + (j-0.5) * dx]
                v(:, k)  = 0
                rho(k)   = 1000
                mass(k)  = rho(k) * dx * dx
                p(k)     = 0
                e(k)     = 0
                itype(k) = 2
                hsml(k)  = dx
            end do
        end do

        ntotal = Nx * Ny

    end subroutine dam_break

    subroutine taylor_rod()
        real(8) :: solid_domain(4) = [0., 0.760, &
                                      0., 2.546] * 1e-2
        real(8) :: dx = 0.038 * 1e-2
        integer :: Nx, Ny

        integer i, j, k

        Nx = floor((solid_domain(2) - solid_domain(1)) / dx)
        Ny = floor((solid_domain(4) - solid_domain(3)) / dx)
        do i = 1, Nx
            do j = 1, Ny
                k = (i-1) * Ny + j
                x(:, k)  = [solid_domain(1) + (i-0.5) * dx, &
                            solid_domain(3) + (j-0.5) * dx]
                v(1, k)  = 0
                v(2, k)  = -221
                rho(k)   = 7850
                mass(k)  = rho(k) * dx * dx
                p(k)     = 0
                e(k)     = 0
                c(k)     = 5000
                itype(k) = 8
                hsml(k)  = dx * 2
            end do
        end do

        ntotal = Nx * Ny

    end subroutine taylor_rod

    subroutine config_input_3()
        integer i, null


        open(11, file = in_path // "/f_xv.dat")
        open(12, file = in_path // "/f_state.dat")
        open(13, file = in_path // "/f_other.dat")

        write(*,*) repeat("#", 60)
        write(*,*) "Loading initial particle configuration..."
        read(11,*) ntotal
        write(*,*) "Total number of particle: ", to_string(ntotal)
        write(*,*) repeat("#", 60)

        do i = 1, ntotal
            read(11,*) null, x(:, i), v(:, i)
            read(12,*) null, mass(i), rho(i), p(i), e(i)
            read(13,*) null, itype(i), hsml(i)
        end do

        close(11)
        close(12)
        close(13)

    end subroutine config_input_3

    subroutine input_3()
        integer i

        open(11, file = in_path // "/ini_xv.dat")
        open(12, file = in_path // "/ini_state.dat")
        open(13, file = in_path // "/ini_other.dat")

        select case(nick)
        case("shear_cavity")
            call shear_cavity()
        case("shock_tube")
            call shock_tube()
        case("tnt_bar")
            call tnt_bar()
        case("tnt_cylinder")
            call tnt_cylinder()
        end select

        do i = 1, ntotal
            write(11, 1001) i, x(:, i), v(:, i)
            write(12, 1002) i, mass(i), rho(i), p(i), e(i)
            write(13, 1003) i, itype(i), hsml(i)
        end do

        1001 format(I5, 6(2X, E15.8))
        1002 format(I5, 4(2X, E15.8))
        1003 format(I5, 2X, I2, 2X, E15.8)

        write(*,*) repeat("#", 60)
        write(*,*) "Initial particle configuration generated"
        write(*,*) "Total number of particles: ", to_string(ntotal)
        write(*,*) repeat("#", 60)

        close(11)
        close(12)
        close(13)

    end subroutine input_3

end module input_m