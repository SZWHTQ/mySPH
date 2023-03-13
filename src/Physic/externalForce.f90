module ex_force_m
    use ctrl_dict, only: dim, i_time_step, gravity_w
    use sph
    use parse_toml_m, only: nick
    implicit none

contains
    subroutine ex_force(P, dvdt)
        type(Particle), intent(in) :: P(:)
        real(8), intent(inout) :: dvdt(:, :)
        integer :: ntotal
        real(8) :: dx(dim), dr, r
        real(8), save :: factor_s, r0, p1, p2
        real(8), save :: factor_p, pe, n1, n2
        logical, save :: first_entry = .true.

        integer i, j, k, d

        ntotal = size(P)
        forall (i=1:dim, j=1:ntotal) dvdt(i, j) = 0

        !!! Consider gravity or not
        if ( gravity_w ) then
            forall (i=1:ntotal) dvdt(dim, i) = -9.8
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

            select case (nick)
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
            case ("taylor_rod")
                factor_s = 1e-2
                r0 = abs(P(2)%x(1) - P(2)%x(2)) * 1
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
                    do d = 2, dim
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
                if ( (P(i)%Type > 0 .and. P(j)%Type > 0) .and. P(i)%Type /= P(j)%Type ) then
                    !!! Calculate the distance 'r' between particle i and j
                    dx(1) = P(i)%x(1) - P(j)%x(1)
                    dr = dx(1)*dx(1)
                    do d = 2, dim
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
                end if
            end do !! k
        end do !! i
        !$OMP END PARALLEL DO


    end subroutine ex_force

end module ex_force_m