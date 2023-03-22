module density_m
    use ctrl_dict, only: Config, Field
    use sph,       only: Particle
    use kernel_m,  only: kernel

    implicit none

contains
    !!! Subroutine to calculate the density with SPH summation algorithm
    subroutine sum_density(P)
        type(Particle), intent(inout) :: P(:)
        integer :: ntotal
        real(8) :: self
        real(8) :: hv(Field%dim)
        real(8), allocatable :: wi(:)  !! Integration of the kernel itself
        integer i, j, k

        ntotal = size(P)
        allocate(wi(ntotal), source=0._8)
        hv = 0

        !!! Self density of each particle: Wii (Kernel for distance 0)
        !!! and take contribution of particle itself

        !!! Firstly, calculate the integration of the kernel over the space
        if ( Config%norm_dens_w ) then
            !$OMP PARALLEL DO PRIVATE(i, self, hv)
            do i = 1, ntotal
                if ( P(i)%State /= 0 ) cycle
                if ( P(i)%divergencePosition < 1.5 ) then
                    call kernel(dble(0), 1*hv, P(i)%SmoothingLength, self, hv)
                    wi(i) = P(i)%Mass / P(i)%Density * self
                end if
            end do
            !$OMP END PARALLEL DO
        end if

        if ( Config%norm_dens_w ) then
            !$OMP PARALLEL DO PRIVATE(i, j, k) REDUCTION(+:wi)
            do i = 1, ntotal
                if ( P(i)%State /= 0 ) cycle
                if ( P(i)%divergencePosition < 1.5 ) then
                    do k = 1, P(i)%neighborNum
                        j = P(i)%neighborList(k)
                        wi(i) = wi(i) + P(j)%Mass/P(j)%Density * P(i)%w(k)
                    end do
                end if
            end do
            !$OMP END PARALLEL DO
        end if

        !!! Secondly, calculate the rho integration over the space
        !$OMP PARALLEL DO PRIVATE(i, self, hv)
        do i = 1, ntotal
            if ( P(i)%State /= 0 ) cycle
            call kernel(dble(0), 1*hv, P(i)%SmoothingLength, self, hv)
            P(i)%Density = P(i)%Mass * self
        end do
        !$OMP END PARALLEL DO

        !!! Calculate SPH sum for rho:
        ! !$OMP PARALLEL DO PRIVATE(i, j, k) REDUCTION(+:rho)
        do i = 1, ntotal
            if ( P(i)%State /= 0 ) cycle
            do k = 1, P(i)%neighborNum
                j = P(i)%neighborList(k)
                P(i)%Density = P(i)%Density + P(j)%Mass * P(i)%w(k)
            end do
        end do
        ! !$OMP END PARALLEL DO

        !!! Thirdly, calculate the normalized rho, rho = Σrho / Σw
        if ( Config%norm_dens_w ) then
            !$OMP PARALLEL DO PRIVATE(i)
            do i = 1, ntotal
                if ( P(i)%State /= 0 ) cycle
                if ( P(i)%divergencePosition < 1.5 ) then
                    P(i)%Density = P(i)%Density / wi(i)
                end if
            end do
            !$OMP END PARALLEL DO
        end if

        deallocate(wi)

    end subroutine sum_density

    !!! Subroutine to calculate the density with SPH continuity approach
    subroutine con_density(P, drhodt)
        type(Particle), intent(in) :: P(:)
        integer :: ntotal
        real(8), intent(inout) :: drhodt(:)    !! Density change rate of each particle

        integer i, j, k

        ntotal = size(P)
        do i = 1, ntotal
            drhodt(i) = 0
        end do

        !$OMP PARALLEL DO PRIVATE(i, j, k) REDUCTION(+:drhodt)
        do i = 1, ntotal
            if ( P(i)%State /= 0 ) cycle
            do k = 1, P(i)%neighborNum
                j = P(i)%neighborList(k)
                drhodt(i) = drhodt(i) &
                    + P(j)%Mass         &
                    * sum((P(i)%v(:) - P(j)%v(:)) * P(i)%dwdx(:, k))
            end do
        end do
        !$OMP END PARALLEL DO


    end subroutine con_density

    !!! PVRS Riemann Solver
    subroutine con_density_riemann(P, drhodt)
        type(Particle), intent(in) :: P(:)
        integer :: ntotal
        real(8), intent(inout) :: drhodt(:)    !! Density change rate of each particle
        real(8) :: Z_l, Z_r, v_l, v_r
        real(8) :: v_ij
        real(8) :: v_star(Field%dim), e_ij(Field%dim)
        integer i, j, k

        ntotal = size(P)
        do i = 1, ntotal
            drhodt(i) = 0
        end do
        v_ij = 0
        e_ij = 0

        do i = 1, ntotal
            if ( P(i)%State /= 0 ) cycle
            do k = 1, P(i)%neighborNum
                j = P(i)%neighborList(k)

                Z_l = P(i)%Density * P(j)%SoundSpeed
                Z_r = P(j)%Density * P(j)%SoundSpeed

                e_ij = (P(j)%x(:) - P(i)%x(:)) / norm2(P(j)%x(:) - P(i)%x(:))

                v_l = dot_product(P(i)%v(:), e_ij)
                v_r = dot_product(P(j)%v(:), e_ij)

                v_ij = ( Z_l*v_l + Z_r*v_r + (P(i)%Pressure-P(j)%Pressure) ) &
                     / ( Z_l + Z_r )
                v_star = v_ij*e_ij + ((P(i)%v(:)+P(j)%v(:))/2 - ((v_l+v_r)/2)*e_ij)

                drhodt(i) = drhodt(i)                           &
                    + 2 * P(i)%Density * P(j)%Mass/P(j)%Density &
                        * dot_product((P(i)%v(:) - v_star), P(i)%dwdx(:, K))
            end do
        end do

    end subroutine con_density_riemann

    subroutine sum_density_dsph(P)
        type(Particle), intent(inout) :: P(:)
        integer :: ntotal
        real(8) :: self
        real(8) :: hv(Field%dim)
        real(8) :: rho_max, rho_min, criteria, ratio
        real(8), allocatable :: wi(:)  !! Integration of the kernel itself
        integer, allocatable :: dc_point(:)
        integer i, j, k, s

        ntotal = size(P)
        allocate(wi(ntotal), source=0._8)
        allocate(dc_point(ntotal), source=0)
        hv = 0

        !!! Self density of each particle: Wii (Kernel for distance 0)
        !!! and take contribution of particle itself

        !!! Firstly, calculate the integration of the kernel over the space
        do i = 1, ntotal
            if ( P(i)%State /= 0 ) cycle
            call kernel(dble(0), 1*hv, P(i)%SmoothingLength, self, hv)
            wi(i) = P(i)%mass/P(i)%Density * self
        end do

        criteria = 0.2
        ratio    = 0
        rho_max  = maxval(P%Density)
        rho_min  = minval(P%Density)
        dc_point = ntotal + 1
        s = 0

        do i = 1, ntotal
            if ( P(i)%State /= 0 ) cycle
            do k = 1, P(i)%neighborNum
                j = P(i)%neighborList(k)
                wi(i) = wi(i) + P(j)%Mass/P(j)%Density * P(i)%w(k)
                if ( abs((P(j)%Density-P(i)%Density)/(rho_max-rho_min)) >= criteria &
               .and. abs((P(j)%Density-P(i)%Density)/(rho_max-rho_min)) >= ratio ) then
                    ratio = abs((P(j)%Density-P(i)%Density)/(rho_max-rho_min))
                    if ( s /= i ) then
                          dc_point(i) = j
                          s = i
                    end if
                 end if
            end do
        end do

        !!! Secondly, calculate the rho integration over the space
        do i = 1, ntotal
            if ( P(i)%State /= 0 ) cycle
            call kernel(dble(0), 1*hv, P(i)%SmoothingLength, self, hv)
            P(i)%Density = P(i)%mass * self
        end do

        do i = 1, ntotal
            if ( P(i)%State /= 0 ) cycle
            do k = 1, P(i)%neighborNum
                j = P(i)%neighborList(k)
                s = dc_point(i)
                P(i)%Density = P(i)%Density + P(j)%Mass*P(i)%w(k)
                if ( j >= s ) then
                    P(i)%Density = P(i)%Density - (P(s)%Density-P(i)%Density) * P(j)%Mass/P(j)%Density * P(i)%w(k)
                end if
            end do
        end do


        !!! Thirdly, calculate the normalized rho, rho = Σrho / Σw
        do i = 1, ntotal
            if ( P(i)%State /= 0 ) cycle
            P(i)%Density = P(i)%Density / wi(i)
        end do

        deallocate(wi, dc_point)

    end subroutine sum_density_dsph

end module density_m