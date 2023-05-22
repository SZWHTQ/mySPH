module ex_force_m
    use ctrl_dict, only: Config, Field, Project
    use sph,       only: Particle
    implicit none

contains
    subroutine external_force(ntotal, P, dvdt)
        integer, intent(in) :: ntotal
        type(Particle), intent(in) :: P(:)
        real(8), intent(inout) :: dvdt(:, :)
        real(8) :: dx(Field%Dim), dr, r
        real(8), save :: factor_s, r0, p1, p2
        real(8), save :: factor_p, pe, n1, n2
        real(8), save :: delta, eta, chi, factor_c
        logical, save :: first_entry = .true.
        real(8) :: alpha, ax, ay

        integer i, j, k, d

        do i=1, Field%Dim
            do j=1, ntotal
                dvdt(i, j) = 0
            end do
        end do

        !!! Consider gravity or not
        if ( Config%gravity_w ) then
            select case (Project%nick)
            case ("beam_oil")
                alpha = slosh(Config%i_time_step * Config%delta_t)
                ax = -9.81 * sin(alpha)
                ay = -9.81 * cos(alpha)
                do i = 1, ntotal
                    if ( P(i)%Type > 100 ) then
                        cycle
                    end if
                    dvdt(1, i) = ax
                    dvdt(2, i) = ay
                end do
            case ("db_gate")
                do i = 1, ntotal
                    if ( P(i)%Type > 100 ) then
                        cycle
                    end if
                    dvdt(Field%Dim, i) = -9.81
                end do
            case default
                do i = 1, ntotal
                    dvdt(Field%Dim, i) = -9.81
                end do
            end select
        end if

        !!! Boundary particle force and penalty anti-penetration force
        if ( first_entry ) then
            factor_s = 1e-2
            r0 = 1.25e-5
            p1 = 12
            p2 = 4
    
            factor_p = 1e5
            n1 = 6
            n2 = 4

            delta = 0.001

            select case (Project%nick)
            case ("shock_tube")
                factor_s = 10
                r0 = abs(P(2)%x(1) - P(1)%x(1)) * 0.5
            case ("tnt_bar")
                factor_s = 1e4
                r0 = abs(P(2)%x(1) - P(1)%x(1))
            case ("undex_chamber")
                factor_s = 1e4
                r0 = abs(P(2)%x(1) - P(1)%x(1)) * 0.5
            case ("dam_break")
                factor_s = 10
                r0 = abs(P(2)%x(1) - P(2)%x(2)) * 1
                delta = 0.2
            case ("taylor_rod")
                factor_s = 1e-2
                r0 = abs(P(2)%x(1) - P(2)%x(2)) * 1
            case ("db_gate")
                ! factor_p = 1e4
                delta = 0.001
            case ("water_impact")
                delta = 0.002
            case ("undex_plate")
                delta = 0.0025
            end select

            first_entry = .false.
        end if

        !$OMP PARALLEL DO PRIVATE(i, j, k, d, dx, dr, r, pe)
        do i = 1, ntotal !! All particles
            do k = 1, P(i)%neighborNum !! All neighbors
                j = P(i)%neighborList(k)

                !!! Interaction between real particle and dummy particle
                if ( P(i)%Type > 0 .and. P(j)%Type < 0 ) then
                    !!! Calculate the distance 'r' between particle i and j
                    dx(1) = P(i)%x(1) - P(j)%x(1)
                    dr = dx(1)*dx(1)
                    do d = 2, Field%Dim
                        dx(d) = P(i)%x(d) - P(j)%x(d)
                        dr = dr + dx(d)*dx(d)
                    end do
                    r = sqrt(dr)
                    !!! Calculate the force between particle i and j,
                    !!! if the distance 'r' is smaller than r0
                    if ( r < r0 ) then
                        dvdt(:, i) = dvdt(:, i) &
                            + factor_s * ((r0/r)**p1 - (r0/r)**p2) * dx(:) / r**2 !/ mass(i)
                    end if
                end if
                !!! Interaction between different phase particles
                if ( (P(i)%Type /= P(j)%Type) .and. (P(i)%Type > 0 .and. P(j)%Type > 0) ) then
                    if ( (P(i)%Type < 100 .and. P(j)%Type < 100) ) then
                        !!! Calculate the distance 'r' between particle i and j
                        dx(1) = P(i)%x(1) - P(j)%x(1)
                        dr = dx(1)*dx(1)
                        do d = 2, Field%Dim
                            dx(d) = P(i)%x(d) - P(j)%x(d)
                            dr = dr + dx(d)*dx(d)
                        end do
                        r = sqrt(dr)
                        !!! Calculate the force between particle i and j
                        pe = (P(i)%SmoothingLength + P(j)%SmoothingLength) / (2 * r)
                        if ( pe >= 1 ) then
                            dvdt(:, i) = dvdt(:, i) &
                                + factor_p * (pe**n1 - pe**n2) * dx(:) / r**2
                        end if
                    else if ( .not. (P(i)%Type > 100 .and. P(j)%Type > 100) ) then
                        !!! Calculate the distance 'r' between particle i and j
                        dx(1) = P(i)%x(1) - P(j)%x(1)
                        dr = dx(1)*dx(1)
                        do d = 2, Field%Dim
                            dx(d) = P(i)%x(d) - P(j)%x(d)
                            dr = dr + dx(d)*dx(d)
                        end do
                        r = sqrt(dr)
                        eta = r / 0.75 / ((P(i)%SmoothingLength + P(j)%SmoothingLength)*0.5)
                        chi = 1 - r / delta
                        if ( chi < 0 .or. chi > 1 ) then
                            chi = 0
                        end if

                        if ( eta <= 0 .or. eta >= 2 ) then
                            factor_c = 0
                        else if ( eta <= 2._8 / 3 ) then
                            factor_c = 2._8 / 3
                        else if ( eta <= 1 ) then
                            factor_c = 2 * eta - 1.5 * eta*eta
                        else if ( eta < 2 ) then
                            factor_c = 0.5*(2-eta)**2
                        end if

                        dvdt(:, i) = dvdt(:, i)               &
                            + 0.01 * P(i)%SoundSpeed**2 * chi &
                                * factor_c * dx(:) / r**2
                    end if
                end if
            end do !! k
        end do !! i
        !$OMP END PARALLEL DO

    contains
        function slosh(t) result(rad)
            use tools_m, only: PI
            real(8), intent(in) :: t
            real(8) :: rad
        
            if ( t < 0.56 ) then
                rad = -143.42*(t**4) + 110.16*(t**3) - 5.866*(t**2) + 0.9767*t + 0.0077
            else
                rad = 4 * sin((2*PI/1.21)*t - 1.3)
            end if

            rad = rad / 180 * PI
            
        end function slosh
    end subroutine external_force

end module ex_force_m